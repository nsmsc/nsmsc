#pragma once

#include "nsmsc/greedy.hpp"
#include "nsmsc/objective_oracle.hpp"
#include "nsmsc/solution.hpp"

#include <cstddef>
#include <cstdint>
#include <span>

namespace nsmsc {

// Completes witnesses and fills unused augmented budget with f-greedy elements.
[[nodiscard]] Solution greedyCompletion(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    GreedyStrategy strategy = GreedyStrategy::Lazy);

// Implements Algorithm 1 without filling budget unused by witness completion.
[[nodiscard]] Solution greedyCompletionOriginal(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    GreedyStrategy strategy = GreedyStrategy::Lazy);

// Runs GCC only at q = ceil((1 - 1/e)h), clipped to the feasible q range.
[[nodiscard]] Solution greedyCompletionFixedQ(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    GreedyStrategy strategy = GreedyStrategy::Lazy);

// Enumerates q as GCC does, but replaces residual-value ranking with a
// reproducibly shuffled component order for each q.
[[nodiscard]] Solution greedyCompletionRandomRanking(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    std::uint64_t seed,
    GreedyStrategy strategy = GreedyStrategy::Lazy);

// Selects k elements by f-greedy, then completes witnesses until the
// complementary value reaches (1 - 1/e) of its global maximum or b is spent.
[[nodiscard]] Solution thresholdTargetedCompletion(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    GreedyStrategy strategy = GreedyStrategy::Lazy);

// Completes the affordable witness with the largest exact F-gain per newly
// added element, then fills unused budget with exact F-greedy elements.
[[nodiscard]] Solution costAwareCompletion(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    GreedyStrategy strategy = GreedyStrategy::Lazy);

// Uses threshold completion to determine b-hat, compares threshold, filled,
// and cost-aware candidates at that budget, and charges all three phases.
[[nodiscard]] Solution adaptiveCompletion(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t maximumAugmentation,
    std::size_t tau,
    GreedyStrategy strategy = GreedyStrategy::Lazy);

// Uses prefixes of an externally certified f-greedy order in GCC.
[[nodiscard]] Solution greedyCompletionFromFOrder(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    std::span<const ElementId> fOrder);

// Uses prefixes of an externally certified f-greedy order for every
// from-empty initialization performed by adaptive completion.
[[nodiscard]] Solution adaptiveCompletionFromFOrder(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t maximumAugmentation,
    std::size_t tau,
    std::span<const ElementId> baseFOrder,
    std::span<const ElementId> augmentedFOrder,
    GreedyStrategy fillStrategy = GreedyStrategy::Lazy);

[[nodiscard]] Solution thresholdTargetedCompletionFromFOrder(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    std::span<const ElementId> fOrder);

[[nodiscard]] Solution costAwareCompletionFromFOrder(
    ObjectiveOracle& oracle,
    std::size_t k,
    std::size_t augmentation,
    std::size_t tau,
    std::span<const ElementId> fOrder);

}  // namespace nsmsc
