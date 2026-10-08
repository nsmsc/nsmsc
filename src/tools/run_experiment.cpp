#include "nsmsc/baselines.hpp"
#include "nsmsc/complementary_completion.hpp"
#include "nsmsc/coverage_oracle.hpp"
#include "nsmsc/experiment_protocol.hpp"
#include "nsmsc/greedy.hpp"
#include "nsmsc/intermediate_io.hpp"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

std::uint64_t parseUnsigned(std::string_view text, const char* name) {
    std::uint64_t value = 0;
    const auto result = std::from_chars(
        text.data(), text.data() + text.size(), value);
    if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
        throw std::invalid_argument(std::string("Invalid ") + name);
    }
    return value;
}

std::size_t parseSize(std::string_view text, const char* name) {
    const std::uint64_t value = parseUnsigned(text, name);
    if (value > std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument(std::string(name) + " is too large");
    }
    return static_cast<std::size_t>(value);
}

double parseDouble(std::string_view text, const char* name) {
    std::size_t consumed = 0;
    const std::string owned(text);
    const double value = std::stod(owned, &consumed);
    if (consumed != owned.size()) {
        throw std::invalid_argument(std::string("Invalid ") + name);
    }
    return value;
}

std::uint64_t resolveSeed(std::string_view text) {
    if (text == "time") {
        return static_cast<std::uint64_t>(
            std::chrono::system_clock::now().time_since_epoch().count());
    }
    return parseUnsigned(text, "seed");
}

std::string jsonEscape(std::string_view text) {
    std::string result;
    result.reserve(text.size());
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

struct RunResult {
    std::string method;
    std::string comparisonGroup;
    std::size_t allowedBudget;
    std::size_t augmentation;
    nsmsc::Solution solution;
    double epsilon = 0.0;
    double lsEpsilon = 0.0;
    double ftgEpsilon = 0.0;
    double ftgAlpha = 0.0;
};

void writeResult(
    std::ostream& output,
    const RunResult& result,
    const std::string& datasetName,
    const std::string& instanceName,
    std::size_t k,
    std::size_t h,
    std::size_t tau,
    std::size_t theoreticalAugmentation,
    std::uint64_t configSeed,
    std::uint64_t algorithmSeed,
    double rewardMinimum,
    double rewardMaximum,
    double rewardTotal,
    double loadSeconds,
    bool includeSelectedElements) {
    output << std::fixed << std::setprecision(9)
           << "{\"dataset\":\"" << jsonEscape(datasetName) << "\""
           << ",\"instance\":\"" << jsonEscape(instanceName) << "\""
           << ",\"method\":\"" << result.method << "\""
           << ",\"comparison_group\":\"" << result.comparisonGroup << "\""
           << ",\"k\":" << k
           << ",\"h\":" << h
           << ",\"tau\":" << tau
           << ",\"theoretical_augmentation\":" << theoreticalAugmentation
           << ",\"augmentation\":" << result.augmentation
           << ",\"allowed_budget\":" << result.allowedBudget
           << ",\"config_seed\":" << configSeed
           << ",\"algorithm_seed\":" << algorithmSeed
           << ",\"component_reward_min\":" << rewardMinimum
           << ",\"component_reward_max\":" << rewardMaximum
           << ",\"component_reward_total\":" << rewardTotal
           << ",\"threshold_epsilon\":" << result.epsilon
           << ",\"ls_epsilon\":" << result.lsEpsilon
           << ",\"ftg_epsilon\":" << result.ftgEpsilon
           << ",\"ftg_alpha\":" << result.ftgAlpha
           << ",\"ls_swaps\":" << result.solution.localSearchSwaps
           << ",\"ls_exchange_evaluations\":" << result.solution.localSearchEvaluations
           << ",\"selected_size\":" << result.solution.elements.size()
           << ",\"completion_count\":" << result.solution.completionCount
           << ",\"winning_branch\":\""
           << jsonEscape(result.solution.winningBranch) << "\""
           << ",\"f_value\":" << result.solution.fValue
           << ",\"complementary_value\":"
           << result.solution.complementaryValue
           << ",\"total_value\":" << result.solution.totalValue
           << ",\"marginal_calls\":" << result.solution.marginalCalls
           << ",\"value_calls\":" << result.solution.valueCalls
           << ",\"load_seconds\":" << loadSeconds
           << ",\"optimization_seconds\":"
           << result.solution.runtimeSeconds;
    if (includeSelectedElements) {
        output << ",\"selected_elements\":[";
        for (std::size_t i = 0; i < result.solution.elements.size(); ++i) {
            if (i != 0) {
                output << ',';
            }
            output << result.solution.elements[i];
        }
        output << ']';
    }
    output << "}\n";
}

}  // namespace

int main(int argc, char** argv) {
    bool omitSelectedElements = false;
    std::string augmentationOverrideText;
    int optionArgc = argc;
    while (optionArgc >= 2) {
        const std::string_view option(argv[optionArgc - 1]);
        if (option == "--omit-selected-elements") {
            if (omitSelectedElements) {
                std::cerr << "Duplicate --omit-selected-elements option\n";
                return 2;
            }
            omitSelectedElements = true;
            --optionArgc;
        } else if (option.starts_with("--augmentation=")) {
            if (!augmentationOverrideText.empty()) {
                std::cerr << "Duplicate --augmentation option\n";
                return 2;
            }
            augmentationOverrideText = option.substr(15);
            --optionArgc;
        } else {
            break;
        }
    }
    if (optionArgc < 6 || optionArgc > 10) {
        std::cerr << "Usage: run_experiment <coverage-input> <config-input> "
                     "<k> <result-output> <seed|time> [threshold-epsilon] "
                     "[ls-epsilon] [ftg-epsilon] [ftg-alpha] "
                     "[--augmentation=<b>] [--omit-selected-elements]\n";
        return 2;
    }
    try {
        double epsilon = 0.1;
        if (optionArgc >= 7) {
            epsilon = parseDouble(argv[6], "threshold epsilon");
        }
        if (!(epsilon > 0 && epsilon < 1))
            throw std::invalid_argument("epsilon must be in (0,1)");
        double lsEpsilon = 0.1;
        if (optionArgc >= 8) {
            lsEpsilon = parseDouble(argv[7], "LS epsilon");
        }
        if (!(lsEpsilon > 0 && lsEpsilon < 1))
            throw std::invalid_argument("LS epsilon must be in (0,1)");
        double ftgEpsilon = 0.1;
        if (optionArgc >= 9) {
            ftgEpsilon = parseDouble(argv[8], "FTG epsilon");
        }
        if (!(ftgEpsilon > 0 && ftgEpsilon < 1))
            throw std::invalid_argument("FTG epsilon must be in (0,1)");
        double ftgAlpha = 1.0;
        if (optionArgc == 10) {
            ftgAlpha = parseDouble(argv[9], "FTG alpha");
        }
        if (!(ftgAlpha >= 1.0))
            throw std::invalid_argument("FTG alpha must be at least 1");
        const auto loadStart = std::chrono::steady_clock::now();
        auto instance = std::make_shared<nsmsc::CoverageInstance>(
            nsmsc::readCoverageBinary(argv[1]));
        auto config = nsmsc::readComplementarityBinary(argv[2], *instance);
        const auto loadEnd = std::chrono::steady_clock::now();

        const std::size_t k = parseSize(argv[3], "k");
        const std::filesystem::path outputPath(argv[4]);
        const std::uint64_t algorithmSeed = resolveSeed(argv[5]);
        const std::size_t h = config.components.size();
        const std::size_t tau = config.tau;
        const std::uint64_t configSeed = config.seed;
        double rewardMinimum = 0.0;
        double rewardMaximum = 0.0;
        double rewardTotal = 0.0;
        if (!config.components.empty()) {
            rewardMinimum = config.components.front().value;
            rewardMaximum = config.components.front().value;
            for (const auto& component : config.components) {
                rewardMinimum = std::min(rewardMinimum, component.value);
                rewardMaximum = std::max(rewardMaximum, component.value);
                rewardTotal += component.value;
            }
        }
        const std::size_t theoreticalAugmentation =
            nsmsc::theoreticalAugmentation(h, tau);
        const std::size_t augmentation = augmentationOverrideText.empty()
            ? theoreticalAugmentation
            : parseSize(augmentationOverrideText, "augmentation");
        if (k > instance->elementCount() ||
            augmentation > instance->elementCount() - k) {
            throw std::invalid_argument(
                "k plus theoretical augmentation exceeds the ground set");
        }
        const std::size_t augmentedBudget = k + augmentation;
        const double loadSeconds = std::chrono::duration<double>(
            loadEnd - loadStart).count();
        nsmsc::CoverageOracle oracle(instance, std::move(config));
        const char* methodFilterValue = std::getenv("METHOD_FILTER");
        const std::string methodFilter =
            methodFilterValue == nullptr ? "" : methodFilterValue;
        const auto enabled = [&](std::string_view method) {
            return methodFilter.empty() || methodFilter == method;
        };

        std::vector<RunResult> results;
        if (enabled("FTG-2022")) {
        results.push_back({
            "FTG-2022", "augmented", augmentedBudget, augmentation,
            nsmsc::fastThresholdGreedyOnTotal(
                oracle, augmentedBudget, ftgEpsilon, ftgAlpha),
            0, 0, ftgEpsilon, ftgAlpha});
        results.push_back({
            "FTG-2022", "original", k, 0,
            nsmsc::fastThresholdGreedyOnTotal(
                oracle, k, ftgEpsilon, ftgAlpha),
            0, 0, ftgEpsilon, ftgAlpha});
        }
        if (enabled("F-Greedy-LS")) {
        results.push_back({"F-Greedy-LS", "augmented", augmentedBudget, augmentation,
            nsmsc::fGreedyLocalSearch(oracle, augmentedBudget, lsEpsilon), 0, lsEpsilon});
        results.push_back({"F-Greedy-LS", "original", k, 0,
            nsmsc::fGreedyLocalSearch(oracle, k, lsEpsilon), 0, lsEpsilon});
        }
        if (enabled("RRG")) {
        results.push_back({
            "RRG", "augmented", augmentedBudget, augmentation,
            nsmsc::residualRandomGreedy(
                oracle, augmentedBudget, algorithmSeed)});
        results.push_back({
            "RRG", "original", k, 0,
            nsmsc::residualRandomGreedy(oracle, k, algorithmSeed)});
        }
        if (enabled("threshold_greedy_2024")) {
        results.push_back({
            "threshold_greedy_2024", "augmented", augmentedBudget, augmentation,
            nsmsc::thresholdGreedyOnTotal(oracle, augmentedBudget, epsilon), epsilon});
        results.push_back({
            "threshold_greedy_2024", "original", k, 0,
            nsmsc::thresholdGreedyOnTotal(oracle, k, epsilon), epsilon});
        }
        if (enabled("proposed")) {
        results.push_back({
            "proposed", "augmented", augmentedBudget, augmentation,
            nsmsc::greedyCompletion(
                oracle, k, augmentation, tau, nsmsc::GreedyStrategy::Lazy)});
        results.push_back({
            "proposed", "original", k, 0,
            nsmsc::greedyCompletion(
                oracle, k, 0, tau, nsmsc::GreedyStrategy::Lazy)});
        }
        if (enabled("proposed-original")) {
        results.push_back({
            "proposed-original", "augmented", augmentedBudget, augmentation,
            nsmsc::greedyCompletionOriginal(
                oracle, k, augmentation, tau, nsmsc::GreedyStrategy::Lazy)});
        results.push_back({
            "proposed-original", "original", k, 0,
            nsmsc::greedyCompletionOriginal(
                oracle, k, 0, tau, nsmsc::GreedyStrategy::Lazy)});
        }
        if (enabled("proposed-no-enumeration")) {
        results.push_back({
            "proposed-no-enumeration", "augmented", augmentedBudget, augmentation,
            nsmsc::greedyCompletionFixedQ(
                oracle, k, augmentation, tau, nsmsc::GreedyStrategy::Lazy)});
        results.push_back({
            "proposed-no-enumeration", "original", k, 0,
            nsmsc::greedyCompletionFixedQ(
                oracle, k, 0, tau, nsmsc::GreedyStrategy::Lazy)});
        }
        if (enabled("proposed-no-ranking")) {
        results.push_back({
            "proposed-no-ranking", "augmented", augmentedBudget, augmentation,
            nsmsc::greedyCompletionRandomRanking(
                oracle, k, augmentation, tau, algorithmSeed,
                nsmsc::GreedyStrategy::Lazy)});
        results.push_back({
            "proposed-no-ranking", "original", k, 0,
            nsmsc::greedyCompletionRandomRanking(
                oracle, k, 0, tau, algorithmSeed,
                nsmsc::GreedyStrategy::Lazy)});
        }
        if (enabled("proposed-threshold")) {
        results.push_back({
            "proposed-threshold", "augmented", augmentedBudget, augmentation,
            nsmsc::thresholdTargetedCompletion(
                oracle, k, augmentation, tau, nsmsc::GreedyStrategy::Lazy)});
        results.push_back({
            "proposed-threshold", "original", k, 0,
            nsmsc::thresholdTargetedCompletion(
                oracle, k, 0, tau, nsmsc::GreedyStrategy::Lazy)});
        }
        if (enabled("proposed-cost-aware")) {
        results.push_back({
            "proposed-cost-aware", "augmented", augmentedBudget, augmentation,
            nsmsc::costAwareCompletion(
                oracle, k, augmentation, tau, nsmsc::GreedyStrategy::Lazy)});
        results.push_back({
            "proposed-cost-aware", "original", k, 0,
            nsmsc::costAwareCompletion(
                oracle, k, 0, tau, nsmsc::GreedyStrategy::Lazy)});
        }
        if (enabled("proposed-adaptive")) {
        auto adaptive = nsmsc::adaptiveCompletion(
            oracle, k, augmentation, tau, nsmsc::GreedyStrategy::Lazy);
        const std::size_t adaptiveAugmentation = adaptive.elements.size() - k;
        results.push_back({
            "proposed-adaptive", "augmented", k + adaptiveAugmentation,
            adaptiveAugmentation, std::move(adaptive)});
        auto adaptiveOriginal = nsmsc::adaptiveCompletion(
            oracle, k, 0, tau, nsmsc::GreedyStrategy::Lazy);
        results.push_back({
            "proposed-adaptive", "original", k, 0,
            std::move(adaptiveOriginal)});
        }
        if (enabled("f_greedy")) {
        results.push_back({
            "f_greedy", "augmented", augmentedBudget, augmentation,
            nsmsc::greedyOnF(
                oracle, augmentedBudget, nsmsc::GreedyStrategy::Lazy)});
        results.push_back({
            "f_greedy", "original", k, 0,
            nsmsc::greedyOnF(oracle, k, nsmsc::GreedyStrategy::Lazy)});
        }
        if (enabled("F_greedy")) {
        results.push_back({
            "F_greedy", "augmented", augmentedBudget, augmentation,
            nsmsc::greedyOnTotal(oracle, augmentedBudget)});
        results.push_back({
            "F_greedy", "original", k, 0,
            nsmsc::greedyOnTotal(oracle, k)});
        }
        if (enabled("random")) {
        results.push_back({
            "random", "augmented", augmentedBudget, augmentation,
            nsmsc::randomSelection(oracle, augmentedBudget, algorithmSeed)});
        results.push_back({
            "random", "original", k, 0,
            nsmsc::randomSelection(oracle, k, algorithmSeed)});
        }
        if (results.empty()) {
            throw std::invalid_argument("Unknown METHOD_FILTER: " + methodFilter);
        }

        if (outputPath.has_parent_path()) {
            std::filesystem::create_directories(outputPath.parent_path());
        }
        std::ofstream output(outputPath, std::ios::trunc);
        if (!output) {
            throw std::runtime_error(
                "Cannot open result file " + outputPath.string());
        }
        const auto instanceDirectory =
            std::filesystem::path(argv[1]).parent_path();
        const std::string instanceName = instanceDirectory.filename().string();
        const std::string datasetName =
            instanceDirectory.has_parent_path()
                ? instanceDirectory.parent_path().filename().string()
                : "unknown";
        for (const auto& result : results) {
            writeResult(
                output, result, datasetName, instanceName, k, h, tau,
                theoreticalAugmentation, configSeed, algorithmSeed,
                rewardMinimum, rewardMaximum, rewardTotal,
                loadSeconds, !omitSelectedElements);
        }
        output.flush();
        if (!output) {
            throw std::runtime_error(
                "Failed to write result file " + outputPath.string());
        }

        std::cout << "output=" << outputPath.string() << '\n'
                  << "k=" << k << '\n'
                  << "h=" << h << '\n'
                  << "tau=" << tau << '\n'
                  << "theoretical_augmentation=" << theoreticalAugmentation << '\n'
                  << "augmentation=" << augmentation << '\n'
                  << "algorithm_seed=" << algorithmSeed << '\n'
                  << "runs=" << results.size() << '\n';
    } catch (const std::exception& error) {
        std::cerr << "run_experiment: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
