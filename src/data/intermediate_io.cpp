#include "nsmsc/intermediate_io.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace nsmsc {
namespace {

constexpr std::array<char, 8> kCoverageMagic = {
    'N', 'S', 'M', 'C', 'O', 'V', '1', '\0'};
constexpr std::array<char, 8> kComplementarityMagic = {
    'N', 'S', 'M', 'C', 'O', 'M', '1', '\0'};
constexpr std::array<char, 8> kComplementarityMagicV2 = {
    'N', 'S', 'M', 'C', 'O', 'M', '2', '\0'};
constexpr std::uint32_t kFormatVersion = 1;
constexpr std::uint32_t kFormatVersionV2 = 2;

static_assert(sizeof(double) == sizeof(std::uint64_t));
static_assert(std::numeric_limits<double>::is_iec559);

template <typename T>
void writeUnsigned(std::ostream& output, T value) {
    static_assert(std::is_unsigned_v<T>);
    std::array<char, sizeof(T)> bytes{};
    for (std::size_t i = 0; i < sizeof(T); ++i) {
        bytes[i] = static_cast<char>(value & 0xffU);
        value >>= 8U;
    }
    output.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
}

template <typename T>
T readUnsigned(std::istream& input) {
    static_assert(std::is_unsigned_v<T>);
    std::array<unsigned char, sizeof(T)> bytes{};
    input.read(reinterpret_cast<char*>(bytes.data()),
               static_cast<std::streamsize>(bytes.size()));
    if (!input) {
        throw std::runtime_error("Unexpected end of intermediate file");
    }
    T value = 0;
    for (std::size_t i = 0; i < sizeof(T); ++i) {
        value |= static_cast<T>(bytes[i]) << (8U * i);
    }
    return value;
}

void writeDouble(std::ostream& output, double value) {
    writeUnsigned(output, std::bit_cast<std::uint64_t>(value));
}

double readDouble(std::istream& input) {
    return std::bit_cast<double>(readUnsigned<std::uint64_t>(input));
}

void writeMagic(std::ostream& output, const std::array<char, 8>& magic) {
    output.write(magic.data(), static_cast<std::streamsize>(magic.size()));
}

void expectMagic(std::istream& input, const std::array<char, 8>& expected) {
    std::array<char, 8> actual{};
    input.read(actual.data(), static_cast<std::streamsize>(actual.size()));
    if (!input || actual != expected) {
        throw std::runtime_error("Invalid intermediate-file magic header");
    }
}

void checkVersion(std::istream& input) {
    const auto version = readUnsigned<std::uint32_t>(input);
    static_cast<void>(readUnsigned<std::uint32_t>(input));
    if (version != kFormatVersion) {
        throw std::runtime_error("Unsupported intermediate-file version");
    }
}

std::size_t checkedSize(std::uint64_t value, const char* field) {
    if (value > std::numeric_limits<std::size_t>::max()) {
        throw std::runtime_error(std::string(field) + " does not fit in memory");
    }
    return static_cast<std::size_t>(value);
}

std::uint64_t checkedAdd(std::uint64_t lhs, std::uint64_t rhs) {
    if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
        throw std::runtime_error("Intermediate-file size overflow");
    }
    return lhs + rhs;
}

std::uint64_t checkedMultiply(std::uint64_t lhs, std::uint64_t rhs) {
    if (lhs != 0 && rhs > std::numeric_limits<std::uint64_t>::max() / lhs) {
        throw std::runtime_error("Intermediate-file size overflow");
    }
    return lhs * rhs;
}

void ensureOutput(std::ofstream& output, const std::filesystem::path& path) {
    output.flush();
    if (!output) {
        throw std::runtime_error("Failed to write " + path.string());
    }
}

void ensureNoTrailingBytes(std::istream& input) {
    if (input.peek() != std::char_traits<char>::eof()) {
        throw std::runtime_error("Intermediate file has trailing bytes");
    }
}

}  // namespace

void writeCoverageBinary(
    const std::filesystem::path& path,
    const CoverageInstance& instance) {
    instance.validate();
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("Cannot open output file " + path.string());
    }

    writeMagic(output, kCoverageMagic);
    writeUnsigned(output, kFormatVersion);
    writeUnsigned(output, std::uint32_t{0});
    writeUnsigned(output, static_cast<std::uint64_t>(instance.targetCount()));
    writeUnsigned(output, static_cast<std::uint64_t>(instance.elementCount()));
    writeUnsigned(output, static_cast<std::uint64_t>(instance.incidenceCount()));
    writeUnsigned(output, instanceFingerprint(instance));
    for (const double weight : instance.targetWeights) {
        writeDouble(output, weight);
    }
    for (const std::uint32_t cost : instance.elementCosts) {
        writeUnsigned(output, cost);
    }
    for (const std::uint64_t offset : instance.offsets) {
        writeUnsigned(output, offset);
    }
    for (const TargetId target : instance.targetIds) {
        writeUnsigned(output, target);
    }
    ensureOutput(output, path);
}

CoverageInstance readCoverageBinary(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Cannot open coverage file " + path.string());
    }
    expectMagic(input, kCoverageMagic);
    checkVersion(input);
    const std::uint64_t targetCount64 = readUnsigned<std::uint64_t>(input);
    const std::uint64_t elementCount64 = readUnsigned<std::uint64_t>(input);
    const std::uint64_t incidenceCount64 = readUnsigned<std::uint64_t>(input);
    const std::uint64_t storedFingerprint = readUnsigned<std::uint64_t>(input);

    std::uint64_t expectedSize = 48;
    expectedSize = checkedAdd(
        expectedSize, checkedMultiply(targetCount64, sizeof(double)));
    expectedSize = checkedAdd(
        expectedSize, checkedMultiply(elementCount64, sizeof(std::uint32_t)));
    expectedSize = checkedAdd(
        expectedSize, checkedMultiply(checkedAdd(elementCount64, 1),
                                      sizeof(std::uint64_t)));
    expectedSize = checkedAdd(
        expectedSize, checkedMultiply(incidenceCount64, sizeof(TargetId)));
    if (std::filesystem::file_size(path) != expectedSize) {
        throw std::runtime_error("Coverage file size does not match its header");
    }

    const std::size_t targetCount = checkedSize(targetCount64, "target count");
    const std::size_t elementCount = checkedSize(elementCount64, "element count");
    const std::size_t incidenceCount =
        checkedSize(incidenceCount64, "incidence count");
    CoverageInstance instance;
    instance.targetWeights.resize(targetCount);
    instance.elementCosts.resize(elementCount);
    instance.offsets.resize(elementCount + 1);
    instance.targetIds.resize(incidenceCount);
    for (double& weight : instance.targetWeights) {
        weight = readDouble(input);
    }
    for (std::uint32_t& cost : instance.elementCosts) {
        cost = readUnsigned<std::uint32_t>(input);
    }
    for (std::uint64_t& offset : instance.offsets) {
        offset = readUnsigned<std::uint64_t>(input);
    }
    for (TargetId& target : instance.targetIds) {
        target = readUnsigned<TargetId>(input);
    }
    ensureNoTrailingBytes(input);
    instance.validate();
    if (instanceFingerprint(instance) != storedFingerprint) {
        throw std::runtime_error("Coverage-file fingerprint mismatch");
    }
    return instance;
}

void writeComplementarityBinary(
    const std::filesystem::path& path,
    const ComplementarityConfig& config,
    const CoverageInstance& instance) {
    config.validate(instance);
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("Cannot open output file " + path.string());
    }

    const bool v2 = std::any_of(
        config.components.begin(), config.components.end(), [](const auto& component) {
            return component.kind == ComponentKind::AlternativeWitnesses;
        });
    writeMagic(output, v2 ? kComplementarityMagicV2 : kComplementarityMagic);
    writeUnsigned(output, v2 ? kFormatVersionV2 : kFormatVersion);
    writeUnsigned(output, std::uint32_t{0});
    writeUnsigned(output, config.seed);
    writeUnsigned(output, static_cast<std::uint64_t>(config.tau));
    writeUnsigned(output, static_cast<std::uint64_t>(config.components.size()));
    writeUnsigned(output, static_cast<std::uint64_t>(instance.targetCount()));
    writeUnsigned(output, static_cast<std::uint64_t>(instance.elementCount()));
    writeUnsigned(output, config.coverageFingerprint);
    for (const auto& component : config.components) {
        if (v2) {
            writeUnsigned(output, static_cast<std::uint32_t>(component.kind));
            writeUnsigned(output, std::uint32_t{0});
        }
        writeUnsigned(output, component.target);
        writeUnsigned(output, component.threshold);
        writeDouble(output, component.value);
        writeUnsigned(
            output, static_cast<std::uint64_t>(component.witness.size()));
        for (const ElementId element : component.witness) {
            writeUnsigned(output, element);
        }
        if (v2) {
            writeUnsigned(output,
                static_cast<std::uint64_t>(component.alternativeCount()));
            for (std::size_t i = 0; i < component.alternativeCount(); ++i) {
                const auto alternative = component.alternative(i);
                writeUnsigned(output,
                    static_cast<std::uint64_t>(alternative.size()));
                for (ElementId element : alternative) writeUnsigned(output, element);
            }
        }
    }
    ensureOutput(output, path);
}

ComplementarityConfig readComplementarityBinary(
    const std::filesystem::path& path,
    const CoverageInstance& instance) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error(
            "Cannot open complementarity file " + path.string());
    }
    std::array<char, 8> magic{};
    input.read(magic.data(), static_cast<std::streamsize>(magic.size()));
    const bool v2 = magic == kComplementarityMagicV2;
    if (!input || (!v2 && magic != kComplementarityMagic)) {
        throw std::runtime_error("Invalid complementarity-file magic header");
    }
    const auto version = readUnsigned<std::uint32_t>(input);
    static_cast<void>(readUnsigned<std::uint32_t>(input));
    if ((!v2 && version != kFormatVersion) ||
        (v2 && version != kFormatVersionV2)) {
        throw std::runtime_error("Unsupported complementarity-file version");
    }

    ComplementarityConfig config;
    config.seed = readUnsigned<std::uint64_t>(input);
    config.tau = checkedSize(readUnsigned<std::uint64_t>(input), "tau");
    const std::size_t componentCount = checkedSize(
        readUnsigned<std::uint64_t>(input), "component count");
    const std::uint64_t targetCount = readUnsigned<std::uint64_t>(input);
    const std::uint64_t elementCount = readUnsigned<std::uint64_t>(input);
    config.coverageFingerprint = readUnsigned<std::uint64_t>(input);
    if (targetCount != instance.targetCount() ||
        elementCount != instance.elementCount()) {
        throw std::runtime_error(
            "Complementarity dimensions do not match coverage data");
    }

    const std::uint64_t fileSize = std::filesystem::file_size(path);
    constexpr std::uint64_t kConfigHeaderSize = 64;
    const std::uint64_t minimumComponentSize = v2 ? 40 : 24;
    if (fileSize < kConfigHeaderSize ||
        componentCount >
            (fileSize - kConfigHeaderSize) / minimumComponentSize) {
        throw std::runtime_error(
            "Complementarity component count is inconsistent with file size");
    }

    config.components.reserve(componentCount);
    for (std::size_t i = 0; i < componentCount; ++i) {
        ThresholdCoverageComponent component;
        if (v2) {
            const auto kind = readUnsigned<std::uint32_t>(input);
            static_cast<void>(readUnsigned<std::uint32_t>(input));
            if (kind > static_cast<std::uint32_t>(
                           ComponentKind::AlternativeWitnesses)) {
                throw std::runtime_error("Unknown component kind in V2 file");
            }
            component.kind = static_cast<ComponentKind>(kind);
        }
        component.target = readUnsigned<TargetId>(input);
        component.threshold = readUnsigned<std::uint32_t>(input);
        component.value = readDouble(input);
        const std::size_t witnessCount = checkedSize(
            readUnsigned<std::uint64_t>(input), "witness count");
        const auto position = input.tellg();
        if (position < 0 ||
            witnessCount >
                (fileSize - static_cast<std::uint64_t>(position)) /
                    sizeof(ElementId)) {
            throw std::runtime_error(
                "Complementarity witness exceeds the remaining file size");
        }
        component.witness.resize(witnessCount);
        for (ElementId& element : component.witness) {
            element = readUnsigned<ElementId>(input);
        }
        if (v2) {
            const std::size_t alternativeCount = checkedSize(
                readUnsigned<std::uint64_t>(input), "alternative count");
            component.alternativeOffsets.reserve(alternativeCount + 1);
            component.alternativeOffsets.push_back(0);
            for (std::size_t j = 0; j < alternativeCount; ++j) {
                const std::size_t size = checkedSize(
                    readUnsigned<std::uint64_t>(input), "alternative size");
                const auto alternativePosition = input.tellg();
                if (alternativePosition < 0 ||
                    size > (fileSize - static_cast<std::uint64_t>(alternativePosition)) /
                               sizeof(ElementId)) {
                    throw std::runtime_error(
                        "Alternative witness exceeds remaining file size");
                }
                for (std::size_t k = 0; k < size; ++k) {
                    component.alternativeElements.push_back(
                        readUnsigned<ElementId>(input));
                }
                component.alternativeOffsets.push_back(
                    component.alternativeElements.size());
            }
        }
        config.components.push_back(std::move(component));
    }
    ensureNoTrailingBytes(input);
    config.validate(instance);
    return config;
}

}  // namespace nsmsc
