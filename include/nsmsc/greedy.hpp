#pragma once

#include "nsmsc/objective_oracle.hpp"
#include "nsmsc/solution.hpp"

#include <cstddef>

namespace nsmsc {

enum class GreedyStrategy {
    Standard,
    Lazy,
};

// Resets the oracle and selects exactly budget elements by greedily maximizing f.
[[nodiscard]] Solution greedyOnF(
    ObjectiveOracle& oracle,
    std::size_t budget,
    GreedyStrategy strategy = GreedyStrategy::Lazy);

}  // namespace nsmsc
