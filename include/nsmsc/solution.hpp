#pragma once

#include "nsmsc/objective_oracle.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace nsmsc {

struct Solution {
    std::vector<ElementId> elements;
    double fValue = 0.0;
    double complementaryValue = 0.0;
    double totalValue = 0.0;
    std::size_t completionCount = 0;
    std::size_t localSearchSwaps = 0;
    std::uint64_t localSearchEvaluations = 0;
    std::uint64_t marginalCalls = 0;
    std::uint64_t valueCalls = 0;
    double runtimeSeconds = 0.0;
    std::string winningBranch;
};

}  // namespace nsmsc
