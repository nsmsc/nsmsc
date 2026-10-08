#pragma once

#include "nsmsc/coverage_instance.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace nsmsc {

struct MipResult {
    bool feasible = false;
    bool optimal = false;
    std::string status;
    std::vector<ElementId> elements;
    double fValue = 0.0;
    double totalValue = 0.0;
    double upperBound = 0.0;
    double gap = 0.0;
    double buildSeconds = 0.0;
    double solveSeconds = 0.0;
};

MipResult solveMip(const CoverageInstance& coverage,
                   const ComplementarityConfig& config,
                   std::size_t budget,
                   double timeLimitSeconds);

}  // namespace nsmsc
