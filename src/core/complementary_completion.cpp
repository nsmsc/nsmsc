#include "nsmsc/complementary_completion.hpp"

#include "greedy_internal.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <numbers>
#include <optional>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

namespace nsmsc {
namespace {

struct Residual {
    double value;
    ComponentId component;
};

enum class ComponentRanking {
    ResidualValue,
    Random,
};

bool residualGreater(const Residual& lhs, const Residual& rhs) {
    if (lhs.value != rhs.value) {
        return lhs.value > rhs.value;
    }
    return lhs.component < rhs.component;
}

std::vector<ElementId> sortedElements(const std::vector<ElementId>& elements) {
    std::vector<ElementId> result = elements;
    std::sort(result.begin(), result.end());
    return result;
}

bool isBetter(const Solution& candidate, const Solution& incumbent) {
    if (candidate.totalValue != incumbent.totalValue) {
        return candidate.totalValue > incumbent.totalValue;
    }
    if (candidate.elements.size() != incumbent.elements.size()) {
        return candidate.elements.size() < incumbent.elements.size();
    }
    return sortedElements(candidate.elements) < sortedElements(incumbent.elements);
}

void validateCompletionInput(
    const ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau) {
    if (tau == 0) {
        throw std::invalid_argument("tau must be positive");
    }
    if (k > oracle.groundSize() ||
        augmentation > oracle.groundSize() - k) {
        throw std::invalid_argument("k + augmentation exceeds ground-set size");
    }
    if (oracle.componentCount() >
        static_cast<std::size_t>(std::numeric_limits<ComponentId>::max()) + 1ULL) {
        throw std::invalid_argument("Too many complementary components");
    }
    for (std::size_t i = 0; i < oracle.componentCount(); ++i) {
        const auto component = static_cast<ComponentId>(i);
        const auto componentWitness = oracle.witness(component);
        if (componentWitness.size() > tau) {
            throw std::invalid_argument("A witness is larger than tau");
        }
        std::unordered_set<ElementId> uniqueElements;
        for (const ElementId element : componentWitness) {
            if (static_cast<std::size_t>(element) >= oracle.groundSize()) {
                throw std::invalid_argument("A witness element is out of range");
            }
            if (!uniqueElements.insert(element).second) {
                throw std::invalid_argument("A witness contains duplicates");
            }
        }
    }
}

detail::GreedyRun initializeFromOrder(
    ObjectiveOracle& oracle,
    std::size_t budget,
    std::span<const ElementId> order) {
    if (order.size() < budget) {
        throw std::invalid_argument("External f-order is shorter than the budget");
    }
    oracle.reset();
    detail::GreedyRun run;
    run.selectionOrder.reserve(budget);
    for (std::size_t i = 0; i < budget; ++i) {
        const ElementId element = order[i];
        if (static_cast<std::size_t>(element) >= oracle.groundSize() ||
            oracle.selected(element)) {
            throw std::invalid_argument("External f-order is invalid");
        }
        oracle.add(element);
        run.selectionOrder.push_back(element);
    }
    return run;
}

}  // namespace

Solution greedyCompletionImpl(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    GreedyStrategy strategy,
    bool fillRemainingBudget,
    std::span<const ElementId> externalOrder = {},
    std::optional<std::size_t> fixedQ = std::nullopt,
    ComponentRanking ranking = ComponentRanking::ResidualValue,
    std::uint64_t rankingSeed = 0) {
    validateCompletionInput(oracle, k, augmentation, tau);
    const auto start = std::chrono::steady_clock::now();
    const std::size_t totalBudget = k + augmentation;
    const std::size_t maxQ = std::min(
        oracle.componentCount(), totalBudget / tau);

    Solution best;
    best.totalValue = -std::numeric_limits<double>::infinity();
    std::uint64_t totalMarginalCalls = 0;
    std::uint64_t totalValueCalls = 0;

    const std::size_t firstQ = fixedQ.value_or(0);
    const std::size_t lastQ = fixedQ.value_or(maxQ);
    if (firstQ > maxQ || lastQ > maxQ) {
        throw std::logic_error("Fixed completion count is infeasible");
    }
    std::mt19937_64 random(rankingSeed);

    for (std::size_t q = firstQ; q <= lastQ; ++q) {
        const std::size_t greedyBudget = totalBudget - q * tau;
        auto run = externalOrder.empty()
            ? (oracle.reset(), detail::runGreedy(oracle, greedyBudget, strategy))
            : initializeFromOrder(oracle, greedyBudget, externalOrder);
        totalMarginalCalls += run.marginalCalls;

        std::vector<ComponentId> rankedComponents;
        rankedComponents.reserve(oracle.componentCount());
        if (ranking == ComponentRanking::ResidualValue) {
            std::vector<Residual> residuals;
            residuals.reserve(oracle.componentCount());
            for (std::size_t i = 0; i < oracle.componentCount(); ++i) {
                const auto component = static_cast<ComponentId>(i);
                const double maximum = oracle.gMaximum(component);
                const double current = oracle.gValue(component);
                totalValueCalls += 2;
                if (!std::isfinite(maximum) || !std::isfinite(current)) {
                    throw std::runtime_error("Oracle returned a non-finite g value");
                }
                if (current > maximum + 1e-12) {
                    throw std::runtime_error("A g value exceeds its declared maximum");
                }
                residuals.push_back({maximum - current, component});
            }
            std::sort(residuals.begin(), residuals.end(), residualGreater);
            for (const auto& residual : residuals) {
                rankedComponents.push_back(residual.component);
            }
        } else {
            for (std::size_t i = 0; i < oracle.componentCount(); ++i) {
                rankedComponents.push_back(static_cast<ComponentId>(i));
            }
            std::shuffle(rankedComponents.begin(), rankedComponents.end(), random);
        }

        std::vector<ElementId> selectionOrder = std::move(run.selectionOrder);
        for (std::size_t index = 0; index < q; ++index) {
            for (const ElementId element :
                 oracle.witness(rankedComponents[index])) {
                if (!oracle.selected(element)) {
                    oracle.add(element);
                    selectionOrder.push_back(element);
                }
            }
        }
        if (oracle.selectedSize() > totalBudget) {
            throw std::logic_error("Completion exceeded the augmented budget");
        }
        if (fillRemainingBudget && oracle.selectedSize() < totalBudget) {
            auto fill = detail::continueGreedy(
                oracle, totalBudget - oracle.selectedSize(), strategy);
            totalMarginalCalls += fill.marginalCalls;
            selectionOrder.insert(
                selectionOrder.end(), fill.selectionOrder.begin(),
                fill.selectionOrder.end());
        }

        Solution candidate;
        candidate.elements = std::move(selectionOrder);
        candidate.completionCount = q;
        candidate.fValue = oracle.fValue();
        ++totalValueCalls;
        for (std::size_t i = 0; i < oracle.componentCount(); ++i) {
            candidate.complementaryValue +=
                oracle.gValue(static_cast<ComponentId>(i));
            ++totalValueCalls;
        }
        candidate.totalValue = candidate.fValue + candidate.complementaryValue;

        if (isBetter(candidate, best)) {
            best = std::move(candidate);
        }
    }

    // Keep the mutable oracle consistent with the returned solution.
    oracle.reset();
    for (const ElementId element : best.elements) {
        if (!oracle.selected(element)) {
            oracle.add(element);
        }
    }

    const auto end = std::chrono::steady_clock::now();
    best.marginalCalls = totalMarginalCalls;
    best.valueCalls = totalValueCalls;
    best.runtimeSeconds = std::chrono::duration<double>(end - start).count();
    return best;
}

Solution greedyCompletion(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    GreedyStrategy strategy) {
    return greedyCompletionImpl(
        oracle, k, augmentation, tau, strategy, true);
}

Solution greedyCompletionOriginal(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    GreedyStrategy strategy) {
    return greedyCompletionImpl(
        oracle, k, augmentation, tau, strategy, false);
}

Solution greedyCompletionFixedQ(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    GreedyStrategy strategy) {
    validateCompletionInput(oracle, k, augmentation, tau);
    const double fraction = 1.0 - 1.0 / std::numbers::e;
    const std::size_t theoreticalQ = static_cast<std::size_t>(
        std::ceil(fraction * static_cast<double>(oracle.componentCount())));
    const std::size_t maxQ = std::min(
        oracle.componentCount(), (k + augmentation) / tau);
    return greedyCompletionImpl(
        oracle, k, augmentation, tau, strategy, true, {},
        std::min(theoreticalQ, maxQ));
}

Solution greedyCompletionRandomRanking(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    std::uint64_t seed,
    GreedyStrategy strategy) {
    return greedyCompletionImpl(
        oracle, k, augmentation, tau, strategy, true, {}, std::nullopt,
        ComponentRanking::Random, seed);
}

Solution greedyCompletionFromFOrder(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    std::span<const ElementId> fOrder) {
    return greedyCompletionImpl(
        oracle, k, augmentation, tau, GreedyStrategy::Lazy, true, fOrder);
}

Solution thresholdTargetedCompletion(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    GreedyStrategy strategy) {
    validateCompletionInput(oracle, k, augmentation, tau);
    const auto start = std::chrono::steady_clock::now();

    oracle.reset();
    auto greedy = detail::runGreedy(oracle, k, strategy);

    std::vector<Residual> residuals;
    residuals.reserve(oracle.componentCount());
    double complementaryMaximum = 0.0;
    std::uint64_t valueCalls = 0;
    for (std::size_t i = 0; i < oracle.componentCount(); ++i) {
        const auto component = static_cast<ComponentId>(i);
        const double maximum = oracle.gMaximum(component);
        const double current = oracle.gValue(component);
        valueCalls += 2;
        if (!std::isfinite(maximum) || !std::isfinite(current)) {
            throw std::runtime_error("Oracle returned a non-finite g value");
        }
        if (current > maximum + 1e-12) {
            throw std::runtime_error("A g value exceeds its declared maximum");
        }
        complementaryMaximum += maximum;
        residuals.push_back({maximum - current, component});
    }
    std::sort(residuals.begin(), residuals.end(), residualGreater);

    const double target = (1.0 - 1.0 / std::numbers::e) * complementaryMaximum;
    const auto complementaryValue = [&]() {
        double value = 0.0;
        for (std::size_t i = 0; i < oracle.componentCount(); ++i) {
            value += oracle.gValue(static_cast<ComponentId>(i));
            ++valueCalls;
        }
        return value;
    };

    double currentComplementary = complementaryValue();
    std::size_t completionCount = 0;
    for (const auto& residual : residuals) {
        if (currentComplementary + 1e-12 >= target) {
            break;
        }
        std::size_t newElementCount = 0;
        for (const ElementId element : oracle.witness(residual.component)) {
            newElementCount += !oracle.selected(element);
        }
        if (newElementCount == 0 ||
            newElementCount > k + augmentation - oracle.selectedSize()) {
            continue;
        }
        for (const ElementId element : oracle.witness(residual.component)) {
            if (!oracle.selected(element)) {
                oracle.add(element);
                greedy.selectionOrder.push_back(element);
            }
        }
        ++completionCount;
        currentComplementary = complementaryValue();
    }

    Solution solution;
    solution.elements = std::move(greedy.selectionOrder);
    solution.completionCount = completionCount;
    solution.fValue = oracle.fValue();
    ++valueCalls;
    solution.complementaryValue = currentComplementary;
    solution.totalValue = solution.fValue + solution.complementaryValue;
    solution.marginalCalls = greedy.marginalCalls;
    solution.valueCalls = valueCalls;
    solution.runtimeSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    return solution;
}

Solution thresholdTargetedCompletionFromOrder(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    GreedyStrategy /*strategy*/,
    std::span<const ElementId> order) {
    validateCompletionInput(oracle, k, augmentation, tau);
    const auto start = std::chrono::steady_clock::now();
    auto greedy = initializeFromOrder(oracle, k, order);
    std::vector<Residual> residuals;
    double maximumTotal = 0.0;
    std::uint64_t valueCalls = 0;
    for (std::size_t i = 0; i < oracle.componentCount(); ++i) {
        const auto component = static_cast<ComponentId>(i);
        const double maximum = oracle.gMaximum(component);
        const double current = oracle.gValue(component);
        maximumTotal += maximum;
        residuals.push_back({maximum - current, component});
        valueCalls += 2;
    }
    std::sort(residuals.begin(), residuals.end(), residualGreater);
    const double target = (1.0 - 1.0 / std::numbers::e) * maximumTotal;
    const auto currentG = [&]() {
        double value = 0.0;
        for (std::size_t i = 0; i < oracle.componentCount(); ++i) {
            value += oracle.gValue(static_cast<ComponentId>(i));
            ++valueCalls;
        }
        return value;
    };
    double complementary = currentG();
    std::size_t completed = 0;
    for (const auto& residual : residuals) {
        if (complementary + 1e-12 >= target) break;
        std::size_t cost = 0;
        for (ElementId element : oracle.witness(residual.component)) {
            cost += !oracle.selected(element);
        }
        if (cost == 0 || cost > k + augmentation - oracle.selectedSize()) continue;
        for (ElementId element : oracle.witness(residual.component)) {
            if (!oracle.selected(element)) {
                oracle.add(element);
                greedy.selectionOrder.push_back(element);
            }
        }
        ++completed;
        complementary = currentG();
    }
    Solution solution;
    solution.elements = std::move(greedy.selectionOrder);
    solution.completionCount = completed;
    solution.fValue = oracle.fValue();
    solution.complementaryValue = complementary;
    solution.totalValue = solution.fValue + complementary;
    solution.valueCalls = valueCalls + 1;
    solution.runtimeSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    return solution;
}

Solution costAwareCompletionImpl(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    GreedyStrategy strategy,
    std::span<const ElementId> externalOrder = {}) {
    validateCompletionInput(oracle, k, augmentation, tau);
    const auto start = std::chrono::steady_clock::now();
    const std::size_t totalBudget = k + augmentation;

    auto greedy = externalOrder.empty()
        ? (oracle.reset(), detail::runGreedy(oracle, k, strategy))
        : initializeFromOrder(oracle, k, externalOrder);

    std::vector<double> maxima(oracle.componentCount());
    double complementaryMaximum = 0.0;
    double currentComplementary = 0.0;
    std::uint64_t valueCalls = 0;
    for (std::size_t i = 0; i < oracle.componentCount(); ++i) {
        const auto component = static_cast<ComponentId>(i);
        maxima[i] = oracle.gMaximum(component);
        const double current = oracle.gValue(component);
        valueCalls += 2;
        if (!std::isfinite(maxima[i]) || !std::isfinite(current)) {
            throw std::runtime_error("Oracle returned a non-finite g value");
        }
        if (current > maxima[i] + 1e-12) {
            throw std::runtime_error("A g value exceeds its declared maximum");
        }
        complementaryMaximum += maxima[i];
        currentComplementary += current;
    }

    const double target =
        (1.0 - 1.0 / std::numbers::e) * complementaryMaximum;
    std::size_t completionCount = 0;
    while (currentComplementary + 1e-12 < target &&
           oracle.selectedSize() < totalBudget) {
        const double currentTotal = oracle.totalValue();
        ++valueCalls;
        bool found = false;
        double bestScore = -std::numeric_limits<double>::infinity();
        double bestGain = -std::numeric_limits<double>::infinity();
        std::size_t bestCost = 0;
        ComponentId bestComponent = 0;
        std::vector<ElementId> bestMissing;

        for (std::size_t i = 0; i < oracle.componentCount(); ++i) {
            const auto component = static_cast<ComponentId>(i);
            const double current = oracle.gValue(component);
            ++valueCalls;
            if (current + 1e-12 >= maxima[i]) {
                continue;
            }

            std::vector<ElementId> missing;
            for (const ElementId element : oracle.witness(component)) {
                if (!oracle.selected(element)) {
                    missing.push_back(element);
                }
            }
            const std::size_t cost = missing.size();
            if (cost == 0 || cost > totalBudget - oracle.selectedSize()) {
                continue;
            }

            for (const ElementId element : missing) {
                oracle.add(element);
            }
            const double gain = oracle.totalValue() - currentTotal;
            ++valueCalls;
            for (auto iterator = missing.rbegin(); iterator != missing.rend();
                 ++iterator) {
                oracle.remove(*iterator);
            }
            if (!std::isfinite(gain)) {
                throw std::runtime_error("Oracle returned a non-finite F gain");
            }
            const double score = gain / static_cast<double>(cost);
            const bool better = !found || score > bestScore + 1e-12 ||
                (std::abs(score - bestScore) <= 1e-12 &&
                 (gain > bestGain + 1e-12 ||
                  (std::abs(gain - bestGain) <= 1e-12 &&
                   (cost < bestCost ||
                    (cost == bestCost && component < bestComponent)))));
            if (better) {
                found = true;
                bestScore = score;
                bestGain = gain;
                bestCost = cost;
                bestComponent = component;
                bestMissing = std::move(missing);
            }
        }

        if (!found) {
            break;
        }
        for (const ElementId element : bestMissing) {
            oracle.add(element);
            greedy.selectionOrder.push_back(element);
        }
        ++completionCount;
        currentComplementary = 0.0;
        for (std::size_t i = 0; i < oracle.componentCount(); ++i) {
            currentComplementary +=
                oracle.gValue(static_cast<ComponentId>(i));
            ++valueCalls;
        }
    }

    if (oracle.selectedSize() < totalBudget) {
        auto fill = detail::continueTotalGreedy(
            oracle, totalBudget - oracle.selectedSize());
        greedy.marginalCalls += fill.marginalCalls;
        greedy.selectionOrder.insert(
            greedy.selectionOrder.end(), fill.selectionOrder.begin(),
            fill.selectionOrder.end());
    }

    Solution solution;
    solution.elements = std::move(greedy.selectionOrder);
    solution.completionCount = completionCount;
    solution.fValue = oracle.fValue();
    ++valueCalls;
    solution.complementaryValue = 0.0;
    for (std::size_t i = 0; i < oracle.componentCount(); ++i) {
        solution.complementaryValue +=
            oracle.gValue(static_cast<ComponentId>(i));
        ++valueCalls;
    }
    solution.totalValue = solution.fValue + solution.complementaryValue;
    solution.marginalCalls = greedy.marginalCalls;
    solution.valueCalls = valueCalls;
    solution.runtimeSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    return solution;
}

Solution costAwareCompletion(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    GreedyStrategy strategy) {
    return costAwareCompletionImpl(
        oracle, k, augmentation, tau, strategy);
}

Solution adaptiveCompletion(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t maximumAugmentation,
    std::size_t tau,
    GreedyStrategy strategy) {
    const auto start = std::chrono::steady_clock::now();
    auto threshold = thresholdTargetedCompletion(
        oracle, k, maximumAugmentation, tau, strategy);
    if (threshold.elements.size() < k) {
        throw std::logic_error("Threshold completion returned fewer than k elements");
    }
    const std::size_t adaptiveAugmentation = threshold.elements.size() - k;
    auto filled = greedyCompletion(
        oracle, k, adaptiveAugmentation, tau, strategy);
    auto costAware = costAwareCompletion(
        oracle, k, adaptiveAugmentation, tau, strategy);

    const std::uint64_t marginalCalls =
        threshold.marginalCalls + filled.marginalCalls +
        costAware.marginalCalls;
    const std::uint64_t valueCalls =
        threshold.valueCalls + filled.valueCalls + costAware.valueCalls;
    threshold.winningBranch = "threshold";
    filled.winningBranch = "original";
    costAware.winningBranch = "cost-aware";
    Solution result = std::move(threshold);
    if (isBetter(filled, result)) {
        result = std::move(filled);
    }
    if (isBetter(costAware, result)) {
        result = std::move(costAware);
    }

    oracle.reset();
    for (const ElementId element : result.elements) {
        oracle.add(element);
    }
    result.marginalCalls = marginalCalls;
    result.valueCalls = valueCalls;
    result.runtimeSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    return result;
}

Solution adaptiveCompletionFromFOrder(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t maximumAugmentation,
    std::size_t tau,
    std::span<const ElementId> baseFOrder,
    std::span<const ElementId> augmentedFOrder,
    GreedyStrategy fillStrategy) {
    validateCompletionInput(oracle, k, maximumAugmentation, tau);
    if (baseFOrder.size() < k) {
        throw std::invalid_argument("Certified base f-order is too short");
    }
    const auto start = std::chrono::steady_clock::now();
    auto threshold = thresholdTargetedCompletionFromOrder(
        oracle, k, maximumAugmentation, tau, fillStrategy, baseFOrder);
    const std::size_t adaptiveAugmentation = threshold.elements.size() - k;
    if (augmentedFOrder.size() < k + adaptiveAugmentation) {
        throw std::invalid_argument(
            "Certified augmented f-order is shorter than the realized budget");
    }
    auto filled = greedyCompletionImpl(
        oracle, k, adaptiveAugmentation, tau, fillStrategy, true,
        augmentedFOrder);
    auto costAware = costAwareCompletionImpl(
        oracle, k, adaptiveAugmentation, tau, fillStrategy, augmentedFOrder);

    const std::uint64_t marginalCalls = threshold.marginalCalls +
        filled.marginalCalls + costAware.marginalCalls;
    const std::uint64_t valueCalls = threshold.valueCalls +
        filled.valueCalls + costAware.valueCalls;
    threshold.winningBranch = "threshold";
    filled.winningBranch = "original";
    costAware.winningBranch = "cost-aware";
    Solution result = std::move(threshold);
    if (isBetter(filled, result)) result = std::move(filled);
    if (isBetter(costAware, result)) result = std::move(costAware);
    oracle.reset();
    for (ElementId element : result.elements) oracle.add(element);
    result.marginalCalls = marginalCalls;
    result.valueCalls = valueCalls;
    result.runtimeSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    return result;
}

Solution thresholdTargetedCompletionFromFOrder(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    std::span<const ElementId> fOrder) {
    return thresholdTargetedCompletionFromOrder(
        oracle, k, augmentation, tau, GreedyStrategy::Lazy, fOrder);
}

Solution costAwareCompletionFromFOrder(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    std::span<const ElementId> fOrder) {
    return costAwareCompletionImpl(
        oracle, k, augmentation, tau, GreedyStrategy::Lazy, fOrder);
}

}  // namespace nsmsc
