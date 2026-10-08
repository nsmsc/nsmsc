#include "nsmsc/intermediate_io.hpp"
#include "nsmsc/mip.hpp"

#include "Highs.h"

#include <charconv>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

std::size_t parseSize(std::string_view text, const char* name) {
    std::uint64_t value = 0;
    const auto result = std::from_chars(
        text.data(), text.data() + text.size(), value);
    if (result.ec != std::errc{} || result.ptr != text.data() + text.size() ||
        value > std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument(std::string("Invalid ") + name);
    }
    return static_cast<std::size_t>(value);
}

double parsePositiveDouble(std::string_view text, const char* name) {
    std::size_t consumed = 0;
    const std::string owned(text);
    const double value = std::stod(owned, &consumed);
    if (consumed != owned.size() || !std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(std::string("Invalid ") + name);
    }
    return value;
}

std::string jsonEscape(std::string_view text) {
    std::string result;
    for (const char character : text) {
        switch (character) {
        case '\\': result += "\\\\"; break;
        case '"': result += "\\\""; break;
        case '\n': result += "\\n"; break;
        case '\r': result += "\\r"; break;
        case '\t': result += "\\t"; break;
        default: result += character; break;
        }
    }
    return result;
}

void writeNumber(std::ostream& output, const char* key, double value) {
    output << ",\"" << key << "\":";
    if (std::isfinite(value)) output << value;
    else output << "null";
}

std::string inferDatasetName(const std::filesystem::path& inputPath) {
    bool followsPreparedData = false;
    for (const auto& part : inputPath) {
        if (followsPreparedData) return part.string();
        followsPreparedData = part == "prepared_data";
    }
    const auto instanceDirectory = inputPath.parent_path();
    return instanceDirectory.has_parent_path()
        ? instanceDirectory.parent_path().filename().string()
        : "unknown";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: run_mip <coverage.bin> <config.bin> "
                     "<total-budget> <time-limit-seconds> <output.jsonl>\n";
        return 2;
    }

    try {
        using Clock = std::chrono::steady_clock;
        const std::filesystem::path coveragePath(argv[1]);
        const std::filesystem::path configPath(argv[2]);
        const std::size_t budget = parseSize(argv[3], "budget");
        const double timeLimit = parsePositiveDouble(argv[4], "time limit");

        const auto loadStart = Clock::now();
        auto coverage = nsmsc::readCoverageBinary(coveragePath);
        auto config = nsmsc::readComplementarityBinary(configPath, coverage);
        const auto loadEnd = Clock::now();
        const double loadSeconds =
            std::chrono::duration<double>(loadEnd - loadStart).count();

        const auto result = nsmsc::solveMip(
            coverage, config, budget, timeLimit);

        const std::filesystem::path outputPath(argv[5]);
        if (outputPath.has_parent_path()) {
            std::filesystem::create_directories(outputPath.parent_path());
        }
        std::ofstream output(outputPath, std::ios::trunc);
        if (!output) {
            throw std::runtime_error(
                "Cannot open output file " + outputPath.string());
        }

        const auto instanceDirectory = coveragePath.parent_path();
        const std::string instanceName =
            instanceDirectory.filename().string();
        const std::string datasetName = inferDatasetName(coveragePath);
        const double missing = std::numeric_limits<double>::quiet_NaN();

        output << std::setprecision(17)
               << "{\"dataset\":\"" << jsonEscape(datasetName) << "\""
               << ",\"instance\":\"" << jsonEscape(instanceName) << "\""
               << ",\"method\":\"MIP-HiGHS\""
               << ",\"comparison_group\":\"same_budget\""
               << ",\"solver_version\":\"" << highsVersion() << "\""
               << ",\"status\":\"" << jsonEscape(result.status) << "\""
               << ",\"optimal\":" << (result.optimal ? "true" : "false")
               << ",\"has_incumbent\":"
               << (result.feasible ? "true" : "false")
               << ",\"allowed_budget\":" << budget
               << ",\"time_limit\":" << timeLimit
               << ",\"threads\":1"
               << ",\"algorithm_seed\":1"
               << ",\"coverage_fingerprint\":"
               << config.coverageFingerprint
               << ",\"selected_size\":" << result.elements.size()
               << ",\"selected_elements\":[";
        for (std::size_t i = 0; i < result.elements.size(); ++i) {
            if (i != 0) output << ',';
            output << result.elements[i];
        }
        output << ']';
        writeNumber(output, "f_value", result.feasible ? result.fValue : missing);
        writeNumber(
            output, "complementary_value",
            result.feasible ? result.totalValue - result.fValue : missing);
        writeNumber(
            output, "total_value",
            result.feasible ? result.totalValue : missing);
        writeNumber(output, "upper_bound", result.upperBound);
        writeNumber(output, "mip_gap", result.gap);
        writeNumber(output, "load_seconds", loadSeconds);
        writeNumber(output, "model_build_seconds", result.buildSeconds);
        writeNumber(output, "solve_seconds", result.solveSeconds);
        writeNumber(
            output, "optimization_seconds",
            result.buildSeconds + result.solveSeconds);
        output << "}\n";
        output.flush();
        if (!output) {
            throw std::runtime_error(
                "Failed to write output file " + outputPath.string());
        }

        std::cout << result.status;
        if (result.feasible) std::cout << " F=" << result.totalValue;
        else std::cout << " no-incumbent";
        std::cout << " bound=" << result.upperBound
                  << " gap=" << result.gap << '\n';
    } catch (const std::exception& error) {
        std::cerr << "run_mip: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
