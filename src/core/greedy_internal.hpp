#pragma once

#include "nsmsc/greedy.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace nsmsc::detail {

struct GreedyRun {
    std::vector<ElementId> selectionOrder;
    std::uint64_t marginalCalls = 0;
};

// Assumes the oracle is empty and leaves it at the greedy solution.
[[nodiscard]] GreedyRun runGreedy(
    ObjectiveOracle& oracle,
    std::size_t budget,
    GreedyStrategy strategy);

// Adds exactly count elements to the oracle's current selection.
[[nodiscard]] GreedyRun continueGreedy(
    ObjectiveOracle& oracle,
    std::size_t count,
    GreedyStrategy strategy);

// Adds exactly count elements using exact marginals of F = f + sum_i g_i.
[[nodiscard]] GreedyRun continueTotalGreedy(
    ObjectiveOracle& oracle,
    std::size_t count);

void validateBudget(const ObjectiveOracle& oracle, std::size_t budget);

}  // namespace nsmsc::detail
