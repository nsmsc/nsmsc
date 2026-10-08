#include "nsmsc/greedy.hpp"

#include "greedy_internal.hpp"

#include <chrono>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nsmsc {
namespace {

struct HeapEntry {
    double gain;
    ElementId element;
};

struct HeapEntryLess {
    bool operator()(const HeapEntry& lhs, const HeapEntry& rhs) const {
        if (lhs.gain != rhs.gain) {
            return lhs.gain < rhs.gain;
        }
        return lhs.element > rhs.element;
    }
};

void validateGain(double gain) {
    if (!std::isfinite(gain)) {
        throw std::runtime_error("Oracle returned a non-finite marginal gain");
    }
}

detail::GreedyRun runStandardGreedy(
    ObjectiveOracle& oracle,
    std::size_t budget) {
    detail::GreedyRun run;
    run.selectionOrder.reserve(budget);

    for (std::size_t step = 0; step < budget; ++step) {
        ElementId bestElement = 0;
        double bestGain = -std::numeric_limits<double>::infinity();
        bool found = false;

        for (std::size_t index = 0; index < oracle.groundSize(); ++index) {
            const auto element = static_cast<ElementId>(index);
            if (oracle.selected(element)) {
                continue;
            }
            const double gain = oracle.fMarginal(element);
            ++run.marginalCalls;
            validateGain(gain);
            if (!found || gain > bestGain ||
                (gain == bestGain && element < bestElement)) {
                bestElement = element;
                bestGain = gain;
                found = true;
            }
        }

        if (!found) {
            throw std::logic_error("No unselected element remains");
        }
        oracle.add(bestElement);
        run.selectionOrder.push_back(bestElement);
    }
    return run;
}

detail::GreedyRun runLazyGreedy(
    ObjectiveOracle& oracle,
    std::size_t budget) {
    detail::GreedyRun run;
    run.selectionOrder.reserve(budget);
    if (budget == 0) {
        return run;
    }
    std::priority_queue<HeapEntry, std::vector<HeapEntry>, HeapEntryLess> heap;

    for (std::size_t index = 0; index < oracle.groundSize(); ++index) {
        const auto element = static_cast<ElementId>(index);
        const double gain = oracle.fMarginal(element);
        ++run.marginalCalls;
        validateGain(gain);
        heap.push({gain, element});
    }

    while (run.selectionOrder.size() < budget) {
        if (heap.empty()) {
            throw std::logic_error("No unselected element remains");
        }

        HeapEntry candidate = heap.top();
        heap.pop();
        if (oracle.selected(candidate.element)) {
            continue;
        }

        candidate.gain = oracle.fMarginal(candidate.element);
        ++run.marginalCalls;
        validateGain(candidate.gain);

        const bool wins = heap.empty() ||
            candidate.gain > heap.top().gain ||
            (candidate.gain == heap.top().gain &&
             candidate.element < heap.top().element);
        if (wins) {
            oracle.add(candidate.element);
            run.selectionOrder.push_back(candidate.element);
        } else {
            heap.push(candidate);
        }
    }
    return run;
}

Solution evaluateGreedySolution(
    ObjectiveOracle& oracle,
    detail::GreedyRun run) {
    Solution solution;
    solution.elements = std::move(run.selectionOrder);
    solution.fValue = oracle.fValue();
    ++solution.valueCalls;
    for (std::size_t i = 0; i < oracle.componentCount(); ++i) {
        solution.complementaryValue +=
            oracle.gValue(static_cast<ComponentId>(i));
        ++solution.valueCalls;
    }
    solution.totalValue = solution.fValue + solution.complementaryValue;
    solution.marginalCalls = run.marginalCalls;
    return solution;
}

}  // namespace

namespace detail {

void validateBudget(const ObjectiveOracle& oracle, std::size_t budget) {
    if (oracle.groundSize() >
        static_cast<std::size_t>(std::numeric_limits<ElementId>::max()) + 1ULL) {
        throw std::invalid_argument("Ground set is too large for ElementId");
    }
    if (budget > oracle.groundSize()) {
        throw std::invalid_argument("Budget exceeds the ground-set size");
    }
}

GreedyRun runGreedy(
    ObjectiveOracle& oracle,
    std::size_t budget,
    GreedyStrategy strategy) {
    validateBudget(oracle, budget);
    if (oracle.selectedSize() != 0) {
        throw std::logic_error("Greedy requires an empty oracle state");
    }
    if (strategy == GreedyStrategy::Standard) {
        return runStandardGreedy(oracle, budget);
    }
    return runLazyGreedy(oracle, budget);
}

GreedyRun continueGreedy(
    ObjectiveOracle& oracle,
    std::size_t count,
    GreedyStrategy strategy) {
    validateBudget(oracle, oracle.selectedSize() + count);
    if (strategy == GreedyStrategy::Standard) {
        return runStandardGreedy(oracle, count);
    }
    return runLazyGreedy(oracle, count);
}

GreedyRun continueTotalGreedy(
    ObjectiveOracle& oracle,
    std::size_t count) {
    validateBudget(oracle, oracle.selectedSize() + count);
    GreedyRun run;
    run.selectionOrder.reserve(count);

    for (std::size_t step = 0; step < count; ++step) {
        ElementId bestElement = 0;
        double bestGain = -std::numeric_limits<double>::infinity();
        bool found = false;
        for (std::size_t index = 0; index < oracle.groundSize(); ++index) {
            const auto element = static_cast<ElementId>(index);
            if (oracle.selected(element)) {
                continue;
            }
            const double gain = oracle.totalMarginal(element);
            ++run.marginalCalls;
            validateGain(gain);
            if (!found || gain > bestGain ||
                (gain == bestGain && element < bestElement)) {
                bestElement = element;
                bestGain = gain;
                found = true;
            }
        }
        if (!found) {
            throw std::logic_error("No unselected element remains");
        }
        oracle.add(bestElement);
        run.selectionOrder.push_back(bestElement);
    }
    return run;
}

}  // namespace detail

Solution greedyOnF(
    ObjectiveOracle& oracle,
    std::size_t budget,
    GreedyStrategy strategy) {
    const auto start = std::chrono::steady_clock::now();
    oracle.reset();
    auto run = detail::runGreedy(oracle, budget, strategy);
    auto solution = evaluateGreedySolution(oracle, std::move(run));
    const auto end = std::chrono::steady_clock::now();
    solution.runtimeSeconds = std::chrono::duration<double>(end - start).count();
    return solution;
}

}  // namespace nsmsc
