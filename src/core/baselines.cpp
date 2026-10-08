#include "nsmsc/baselines.hpp"

#include "greedy_internal.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numbers>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

namespace nsmsc {
namespace {

Solution evaluateCurrent(
    ObjectiveOracle& oracle,
    std::vector<ElementId> elements,
    std::uint64_t marginalCalls,
    const std::chrono::steady_clock::time_point start) {
    Solution solution;
    solution.elements = std::move(elements);
    solution.fValue = oracle.fValue();
    ++solution.valueCalls;
    for (std::size_t i = 0; i < oracle.componentCount(); ++i) {
        solution.complementaryValue +=
            oracle.gValue(static_cast<ComponentId>(i));
        ++solution.valueCalls;
    }
    solution.totalValue = solution.fValue + solution.complementaryValue;
    solution.marginalCalls = marginalCalls;
    solution.runtimeSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    return solution;
}

}  // namespace

Solution greedyOnTotal(ObjectiveOracle& oracle, std::size_t budget) {
    detail::validateBudget(oracle, budget);
    const auto start = std::chrono::steady_clock::now();
    oracle.reset();
    std::vector<ElementId> elements;
    elements.reserve(budget);
    std::uint64_t marginalCalls = 0;

    for (std::size_t step = 0; step < budget; ++step) {
        bool found = false;
        ElementId bestElement = 0;
        double bestGain = -std::numeric_limits<double>::infinity();
        for (std::size_t index = 0; index < oracle.groundSize(); ++index) {
            const auto element = static_cast<ElementId>(index);
            if (oracle.selected(element)) {
                continue;
            }
            const double gain = oracle.totalMarginal(element);
            ++marginalCalls;
            if (!std::isfinite(gain)) {
                throw std::runtime_error(
                    "Oracle returned a non-finite total marginal gain");
            }
            if (!found || gain > bestGain ||
                (gain == bestGain && element < bestElement)) {
                found = true;
                bestElement = element;
                bestGain = gain;
            }
        }
        if (!found) {
            throw std::logic_error("No unselected element remains");
        }
        oracle.add(bestElement);
        elements.push_back(bestElement);
    }
    return evaluateCurrent(oracle, std::move(elements), marginalCalls, start);
}

Solution thresholdGreedyOnTotal(
    ObjectiveOracle& oracle, std::size_t budget, double epsilon) {
    detail::validateBudget(oracle, budget);
    if (!std::isfinite(epsilon) || epsilon <= 0 || epsilon >= 1) {
        throw std::invalid_argument("epsilon must be in (0,1)");
    }
    // Algorithm 1 assumes k >= 2; solve the smaller budgets directly.
    if (budget < 2) return greedyOnTotal(oracle, budget);
    const auto start = std::chrono::steady_clock::now();
    oracle.reset();
    std::vector<ElementId> elements;
    std::uint64_t calls = 0;
    auto marginal = [&](ElementId e) {
        const double gain = oracle.totalMarginal(e);
        ++calls;
        if (!std::isfinite(gain) || gain < 0)
            throw std::runtime_error("Expected finite nonnegative marginal");
        return gain;
    };
    double d = 0;
    for (std::size_t i = 0; i < oracle.groundSize(); ++i)
        d = std::max(d, marginal(static_cast<ElementId>(i)));
    // A zero initial threshold would never decrease. No singleton is positive.
    if (d == 0) return evaluateCurrent(oracle, {}, calls, start);
    const double logK = std::log(static_cast<double>(budget));
    const double eta = epsilon / (2 + logK);
    const double decay = 1 - eta;
    const double stop = eta * eta * decay * d / (4 * budget * logK);
    double threshold = d / decay;
    if (decay == 1 || stop <= 0 || !std::isfinite(threshold))
        throw std::invalid_argument("Threshold parameters exceed floating-point range");
    while (threshold >= stop && elements.size() < budget) {
        for (std::size_t i = 0; i < oracle.groundSize() && elements.size() < budget; ++i) {
            const auto e = static_cast<ElementId>(i);
            if (!oracle.selected(e) && marginal(e) >= threshold) {
                oracle.add(e);
                elements.push_back(e);
            }
        }
        threshold *= decay;
    }
    return evaluateCurrent(oracle, std::move(elements), calls, start);
}

Solution fGreedyLocalSearch(ObjectiveOracle& oracle, std::size_t budget, double epsilon) {
    if (!std::isfinite(epsilon) || epsilon <= 0 || epsilon >= 1)
        throw std::invalid_argument("LS epsilon must be in (0,1)");
    const auto start = std::chrono::steady_clock::now();
    auto initial = greedyOnTotal(oracle, budget);
    auto elements = initial.elements;
    const double eta = epsilon / std::max<std::size_t>(1, oracle.groundSize());
    std::uint64_t marginalCalls = initial.marginalCalls;
    std::uint64_t valueCalls = initial.valueCalls;
    std::uint64_t evaluations = 0;
    std::size_t swaps = 0;
    bool improved;
    do {
        improved = false;
        // Scan slots in initialization order; incoming elements in ID order.
        for (std::size_t slot = 0; slot < elements.size(); ++slot) {
            const auto outgoing = elements[slot];
            const double current = oracle.totalValue();
            ++valueCalls;
            oracle.remove(outgoing);
            const double reduced = oracle.totalValue();
            ++valueCalls;
            bool accepted = false;
            for (std::size_t i = 0; i < oracle.groundSize(); ++i) {
                const auto incoming = static_cast<ElementId>(i);
                if (incoming == outgoing || oracle.selected(incoming)) continue;
                const double gain = oracle.totalMarginal(incoming);
                ++marginalCalls;
                ++evaluations;
                if (!std::isfinite(gain)) {
                    oracle.add(outgoing);
                    throw std::runtime_error("Non-finite exchange marginal");
                }
                if (reduced + gain - current > eta * current) {
                    oracle.add(incoming);
                    elements[slot] = incoming;
                    ++swaps;
                    accepted = improved = true;
                    break;
                }
            }
            if (!accepted) oracle.add(outgoing);
        }
    } while (improved);
    auto result = evaluateCurrent(oracle, std::move(elements), marginalCalls, start);
    result.valueCalls += valueCalls;
    result.localSearchSwaps = swaps;
    result.localSearchEvaluations = evaluations;
    return result;
}

Solution localSearchFromOrder(
    ObjectiveOracle& oracle,
    std::span<const ElementId> order,
    std::size_t budget,
    double epsilon,
    std::uint64_t maxExchangeEvaluations,
    double maxSeconds) {
    detail::validateBudget(oracle, budget);
    if (!std::isfinite(epsilon) || epsilon <= 0 || epsilon >= 1 ||
        !std::isfinite(maxSeconds) || maxSeconds <= 0 ||
        order.size() < budget) {
        throw std::invalid_argument("Invalid external-order LS parameters");
    }
    const auto start = std::chrono::steady_clock::now();
    oracle.reset();
    std::vector<ElementId> elements(order.begin(), order.begin() + budget);
    for (ElementId element : elements) {
        if (static_cast<std::size_t>(element) >= oracle.groundSize() ||
            oracle.selected(element)) {
            throw std::invalid_argument("External LS order is invalid");
        }
        oracle.add(element);
    }
    const double eta = epsilon / std::max<std::size_t>(1, oracle.groundSize());
    std::uint64_t marginalCalls = 0;
    std::uint64_t valueCalls = 0;
    std::uint64_t evaluations = 0;
    std::size_t swaps = 0;
    bool improved = false;
    bool limited = false;
    do {
        improved = false;
        for (std::size_t slot = 0; slot < elements.size(); ++slot) {
            const auto outgoing = elements[slot];
            const double current = oracle.totalValue();
            ++valueCalls;
            oracle.remove(outgoing);
            const double reduced = oracle.totalValue();
            ++valueCalls;
            bool accepted = false;
            for (std::size_t i = 0; i < oracle.groundSize(); ++i) {
                if (evaluations >= maxExchangeEvaluations ||
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - start).count() >=
                        maxSeconds) {
                    limited = true;
                    break;
                }
                const auto incoming = static_cast<ElementId>(i);
                if (incoming == outgoing || oracle.selected(incoming)) continue;
                const double gain = oracle.totalMarginal(incoming);
                ++marginalCalls;
                ++evaluations;
                if (!std::isfinite(gain)) {
                    oracle.add(outgoing);
                    throw std::runtime_error("Non-finite exchange marginal");
                }
                if (reduced + gain - current > eta * current) {
                    oracle.add(incoming);
                    elements[slot] = incoming;
                    ++swaps;
                    accepted = improved = true;
                    break;
                }
            }
            if (!accepted) oracle.add(outgoing);
            if (limited) break;
        }
    } while (improved && !limited);
    auto result = evaluateCurrent(
        oracle, std::move(elements), marginalCalls, start);
    result.valueCalls += valueCalls;
    result.localSearchSwaps = swaps;
    result.localSearchEvaluations = evaluations;
    result.winningBranch = limited ? "limit" : "local-optimum";
    return result;
}

Solution randomSelection(
    ObjectiveOracle& oracle,
    std::size_t budget,
    std::uint64_t seed) {
    detail::validateBudget(oracle, budget);
    const auto start = std::chrono::steady_clock::now();
    oracle.reset();
    std::vector<ElementId> candidates(oracle.groundSize());
    std::iota(candidates.begin(), candidates.end(), ElementId{0});
    std::mt19937_64 random(seed);
    std::shuffle(candidates.begin(), candidates.end(), random);
    candidates.resize(budget);
    for (const ElementId element : candidates) {
        oracle.add(element);
    }
    return evaluateCurrent(oracle, std::move(candidates), 0, start);
}

Solution residualRandomGreedy(
    ObjectiveOracle& oracle,
    std::size_t budget,
    std::uint64_t seed) {
    detail::validateBudget(oracle, budget);
    const auto start = std::chrono::steady_clock::now();
    oracle.reset();
    std::vector<ElementId> elements;
    elements.reserve(budget);
    std::mt19937_64 random(seed);
    std::uint64_t marginalCalls = 0;

    struct Candidate {
        double gain;
        ElementId element;
    };
    for (std::size_t step = 0; step < budget; ++step) {
        std::vector<Candidate> candidates;
        candidates.reserve(oracle.groundSize() - step);
        for (std::size_t index = 0; index < oracle.groundSize(); ++index) {
            const auto element = static_cast<ElementId>(index);
            if (oracle.selected(element)) {
                continue;
            }
            const double gain = oracle.totalMarginal(element);
            ++marginalCalls;
            if (!std::isfinite(gain)) {
                throw std::runtime_error(
                    "Oracle returned a non-finite total marginal gain");
            }
            candidates.push_back({gain, element});
        }

        const std::size_t residualRank = budget - step;
        std::partial_sort(
            candidates.begin(), candidates.begin() + residualRank,
            candidates.end(), [](const Candidate& lhs, const Candidate& rhs) {
                return lhs.gain > rhs.gain ||
                    (lhs.gain == rhs.gain && lhs.element < rhs.element);
            });
        std::uniform_int_distribution<std::size_t> choice(0, residualRank - 1);
        const ElementId selected = candidates[choice(random)].element;
        oracle.add(selected);
        elements.push_back(selected);
    }
    return evaluateCurrent(
        oracle, std::move(elements), marginalCalls, start);
}

Solution fastThresholdGreedyOnTotal(
    ObjectiveOracle& oracle,
    std::size_t budget,
    double epsilon,
    double alpha) {
    detail::validateBudget(oracle, budget);
    if (!std::isfinite(epsilon) || epsilon <= 0.0 || epsilon >= 1.0) {
        throw std::invalid_argument("FTG epsilon must be in (0,1)");
    }
    if (!std::isfinite(alpha) || alpha < 1.0) {
        throw std::invalid_argument("FTG alpha must be finite and at least 1");
    }
    const auto start = std::chrono::steady_clock::now();
    oracle.reset();
    if (budget == 0) {
        return evaluateCurrent(oracle, {}, 0, start);
    }

    std::uint64_t marginalCalls = 0;
    double estimateValue = 0.0;
    // Algorithm 1 with uniform cost c(e)=1/budget. This auxiliary set is
    // unconstrained and is used only to estimate OPT's value scale.
    for (std::size_t index = 0; index < oracle.groundSize(); ++index) {
        const auto element = static_cast<ElementId>(index);
        const double gain = oracle.totalMarginal(element);
        ++marginalCalls;
        if (!std::isfinite(gain) || gain < 0.0) {
            throw std::runtime_error(
                "FTG expected a finite nonnegative total marginal");
        }
        if (gain * static_cast<double>(budget) >= estimateValue) {
            oracle.add(element);
            estimateValue += gain;
        }
    }
    const double gamma = estimateValue / 4.0;

    oracle.reset();
    std::vector<ElementId> elements;
    elements.reserve(budget);
    if (gamma == 0.0) {
        return evaluateCurrent(
            oracle, std::move(elements), marginalCalls, start);
    }
    double threshold = 8.0 * alpha * gamma;
    const double stop = (1.0 - epsilon) * gamma / std::numbers::e;
    while (threshold > stop && elements.size() < budget) {
        for (std::size_t index = 0;
             index < oracle.groundSize() && elements.size() < budget;
             ++index) {
            const auto element = static_cast<ElementId>(index);
            if (oracle.selected(element)) {
                continue;
            }
            const double gain = oracle.totalMarginal(element);
            ++marginalCalls;
            if (!std::isfinite(gain) || gain < 0.0) {
                throw std::runtime_error(
                    "FTG expected a finite nonnegative total marginal");
            }
            if (gain * static_cast<double>(budget) >= threshold) {
                oracle.add(element);
                elements.push_back(element);
            }
        }
        threshold *= 1.0 - epsilon;
    }
    return evaluateCurrent(
        oracle, std::move(elements), marginalCalls, start);
}

}  // namespace nsmsc
