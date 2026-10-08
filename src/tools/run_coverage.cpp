#include "nsmsc/complementary_completion.hpp"
#include "nsmsc/coverage_oracle.hpp"
#include "nsmsc/intermediate_io.hpp"

#include <charconv>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
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

nsmsc::GreedyStrategy parseStrategy(std::string_view text) {
    if (text == "lazy") {
        return nsmsc::GreedyStrategy::Lazy;
    }
    if (text == "standard") {
        return nsmsc::GreedyStrategy::Standard;
    }
    throw std::invalid_argument("Strategy must be lazy or standard");
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: run_coverage <coverage-input> <config-input> "
                     "<k> <augmentation> <lazy|standard>\n";
        return 2;
    }
    try {
        const auto start = std::chrono::steady_clock::now();
        auto instance = std::make_shared<nsmsc::CoverageInstance>(
            nsmsc::readCoverageBinary(argv[1]));
        auto config = nsmsc::readComplementarityBinary(argv[2], *instance);
        const auto loaded = std::chrono::steady_clock::now();
        const std::size_t k = parseSize(argv[3], "k");
        const std::size_t augmentation = parseSize(argv[4], "augmentation");
        const auto strategy = parseStrategy(argv[5]);
        const std::size_t tau = config.tau;
        const std::size_t h = config.components.size();

        nsmsc::CoverageOracle oracle(instance, std::move(config));
        const auto solution = nsmsc::greedyCompletion(
            oracle, k, augmentation, tau, strategy);
        const double loadSeconds =
            std::chrono::duration<double>(loaded - start).count();

        std::cout << std::fixed << std::setprecision(9)
                  << "{\"targets\":" << instance->targetCount()
                  << ",\"elements\":" << instance->elementCount()
                  << ",\"h\":" << h
                  << ",\"tau\":" << tau
                  << ",\"k\":" << k
                  << ",\"augmentation\":" << augmentation
                  << ",\"selected_size\":" << solution.elements.size()
                  << ",\"completion_count\":" << solution.completionCount
                  << ",\"f_value\":" << solution.fValue
                  << ",\"complementary_value\":"
                  << solution.complementaryValue
                  << ",\"total_value\":" << solution.totalValue
                  << ",\"marginal_calls\":" << solution.marginalCalls
                  << ",\"value_calls\":" << solution.valueCalls
                  << ",\"load_seconds\":" << loadSeconds
                  << ",\"optimization_seconds\":"
                  << solution.runtimeSeconds << "}\n";
    } catch (const std::exception& error) {
        std::cerr << "run_coverage: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
