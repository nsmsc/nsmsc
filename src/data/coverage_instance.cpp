#include "nsmsc/coverage_instance.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <unordered_set>

namespace nsmsc {
namespace {

constexpr std::uint64_t kFnvOffset = 14695981039346656037ULL;
constexpr std::uint64_t kFnvPrime = 1099511628211ULL;

void hashByte(std::uint64_t& hash, std::uint8_t value) {
    hash ^= value;
    hash *= kFnvPrime;
}

template <typename T>
void hashInteger(std::uint64_t& hash, T value) {
    static_assert(std::is_unsigned_v<T>);
    for (std::size_t i = 0; i < sizeof(T); ++i) {
        hashByte(hash, static_cast<std::uint8_t>(value & 0xffU));
        value >>= 8U;
    }
}

}  // namespace

std::size_t ThresholdCoverageComponent::alternativeCount() const {
    return alternativeOffsets.empty() ? 0 : alternativeOffsets.size() - 1;
}

std::span<const ElementId> ThresholdCoverageComponent::alternative(
    std::size_t index) const {
    if (index >= alternativeCount()) {
        throw std::out_of_range("Alternative witness is out of range");
    }
    const auto begin = static_cast<std::size_t>(alternativeOffsets[index]);
    const auto size = static_cast<std::size_t>(
        alternativeOffsets[index + 1] - alternativeOffsets[index]);
    return {alternativeElements.data() + begin, size};
}

std::size_t CoverageInstance::targetCount() const {
    return targetWeights.size();
}

std::size_t CoverageInstance::elementCount() const {
    return elementCosts.size();
}

std::size_t CoverageInstance::incidenceCount() const {
    return targetIds.size();
}

std::span<const TargetId> CoverageInstance::targets(ElementId element) const {
    const auto index = static_cast<std::size_t>(element);
    if (index >= elementCount()) {
        throw std::out_of_range("ElementId is out of range");
    }
    const auto begin = static_cast<std::size_t>(offsets[index]);
    const auto length = static_cast<std::size_t>(offsets[index + 1] - offsets[index]);
    return std::span<const TargetId>(targetIds.data() + begin, length);
}

void CoverageInstance::validate() const {
    if (targetCount() >
            static_cast<std::size_t>(std::numeric_limits<TargetId>::max()) + 1ULL ||
        elementCount() >
            static_cast<std::size_t>(std::numeric_limits<ElementId>::max()) + 1ULL) {
        throw std::invalid_argument("Coverage instance exceeds 32-bit IDs");
    }
    if (offsets.size() != elementCount() + 1 || offsets.empty() ||
        offsets.front() != 0 || offsets.back() != targetIds.size()) {
        throw std::invalid_argument("Invalid CSR offsets");
    }
    for (const double weight : targetWeights) {
        if (!std::isfinite(weight) || weight < 0.0) {
            throw std::invalid_argument(
                "Target weights must be finite and nonnegative");
        }
    }
    for (const std::uint32_t cost : elementCosts) {
        if (cost == 0) {
            throw std::invalid_argument("Element costs must be positive");
        }
    }
    for (std::size_t element = 0; element < elementCount(); ++element) {
        if (offsets[element] > offsets[element + 1]) {
            throw std::invalid_argument("CSR offsets must be nondecreasing");
        }
        bool first = true;
        TargetId previous = 0;
        for (const TargetId target : targets(static_cast<ElementId>(element))) {
            if (static_cast<std::size_t>(target) >= targetCount()) {
                throw std::invalid_argument("TargetId is out of range");
            }
            if (!first && target <= previous) {
                throw std::invalid_argument(
                    "Targets for each element must be sorted and unique");
            }
            first = false;
            previous = target;
        }
    }
}

void ComplementarityConfig::validate(const CoverageInstance& instance) const {
    instance.validate();
    if (tau == 0) {
        throw std::invalid_argument("Configuration tau must be positive");
    }
    if (coverageFingerprint != instanceFingerprint(instance)) {
        throw std::invalid_argument(
            "Complementarity configuration does not match coverage data");
    }
    if (components.size() >
        static_cast<std::size_t>(std::numeric_limits<ComponentId>::max()) + 1ULL) {
        throw std::invalid_argument("Too many complementary components");
    }

    std::unordered_set<TargetId> componentTargets;
    for (const auto& component : components) {
        if (!std::isfinite(component.value) || component.value < 0.0) {
            throw std::invalid_argument(
                "Component values must be finite and nonnegative");
        }

        if (component.witness.empty() || component.witness.size() > tau) {
            throw std::invalid_argument("Canonical witness size is invalid");
        }
        std::unordered_set<ElementId> witnessElements;
        for (const ElementId element : component.witness) {
            if (static_cast<std::size_t>(element) >= instance.elementCount()) {
                throw std::invalid_argument("Witness element is out of range");
            }
            if (!witnessElements.insert(element).second) {
                throw std::invalid_argument("A witness contains duplicates");
            }
        }
        if (component.kind == ComponentKind::ThresholdCoverage) {
            if (static_cast<std::size_t>(component.target) >= instance.targetCount()) {
                throw std::invalid_argument("Component target is out of range");
            }
            if (!componentTargets.insert(component.target).second) {
                throw std::invalid_argument(
                    "A target has more than one complementary component");
            }
            if (component.threshold == 0 || component.threshold > tau ||
                component.witness.size() != component.threshold ||
                !component.alternativeOffsets.empty() ||
                !component.alternativeElements.empty()) {
                throw std::invalid_argument("Invalid threshold component");
            }
            for (const ElementId element : component.witness) {
                const auto covered = instance.targets(element);
                if (!std::binary_search(covered.begin(), covered.end(),
                                        component.target)) {
                    throw std::invalid_argument(
                        "A witness element does not cover its component target");
                }
            }
        } else if (component.kind == ComponentKind::AlternativeWitnesses) {
            if (component.alternativeOffsets.size() < 2 ||
                component.alternativeOffsets.front() != 0 ||
                component.alternativeOffsets.back() !=
                    component.alternativeElements.size()) {
                throw std::invalid_argument("Invalid alternative-witness offsets");
            }
            std::vector<ElementId> previous;
            bool firstAlternative = true;
            bool canonicalFound = false;
            for (std::size_t i = 0; i < component.alternativeCount(); ++i) {
                const auto alternative = component.alternative(i);
                if (alternative.empty() || alternative.size() > tau ||
                    !std::is_sorted(alternative.begin(), alternative.end()) ||
                    std::adjacent_find(alternative.begin(), alternative.end()) !=
                        alternative.end()) {
                    throw std::invalid_argument("Invalid alternative witness");
                }
                for (ElementId element : alternative) {
                    if (static_cast<std::size_t>(element) >= instance.elementCount()) {
                        throw std::invalid_argument(
                            "Alternative witness element is out of range");
                    }
                }
                const std::vector<ElementId> current(alternative.begin(), alternative.end());
                if (!firstAlternative && !(previous < current)) {
                    throw std::invalid_argument(
                        "Alternative witnesses must be sorted and unique");
                }
                canonicalFound = canonicalFound || current == component.witness;
                previous = current;
                firstAlternative = false;
            }
            if (!canonicalFound) {
                throw std::invalid_argument(
                    "Canonical witness is not an alternative");
            }
        } else {
            throw std::invalid_argument("Unknown component kind");
        }
    }
}

std::uint64_t instanceFingerprint(const CoverageInstance& instance) {
    instance.validate();
    std::uint64_t hash = kFnvOffset;
    hashInteger(hash, static_cast<std::uint64_t>(instance.targetCount()));
    hashInteger(hash, static_cast<std::uint64_t>(instance.elementCount()));
    hashInteger(hash, static_cast<std::uint64_t>(instance.incidenceCount()));
    for (const double weight : instance.targetWeights) {
        hashInteger(hash, std::bit_cast<std::uint64_t>(weight));
    }
    for (const std::uint32_t cost : instance.elementCosts) {
        hashInteger(hash, cost);
    }
    for (const std::uint64_t offset : instance.offsets) {
        hashInteger(hash, offset);
    }
    for (const TargetId target : instance.targetIds) {
        hashInteger(hash, target);
    }
    return hash;
}

}  // namespace nsmsc
