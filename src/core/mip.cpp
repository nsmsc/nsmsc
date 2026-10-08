#include "nsmsc/mip.hpp"

#include "nsmsc/coverage_oracle.hpp"

#include "Highs.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace nsmsc {

MipResult solveMip(const CoverageInstance& coverage,
                   const ComplementarityConfig& config,
                   std::size_t budget,
                   double timeLimitSeconds) {
    using Clock = std::chrono::steady_clock;

    const auto start = Clock::now();
    coverage.validate();
    config.validate(coverage);
    if (budget > coverage.elementCount() ||
        !std::isfinite(timeLimitSeconds) || timeLimitSeconds <= 0.0) {
        throw std::invalid_argument("Invalid budget or time limit");
    }
    const auto elementCount = coverage.elementCount();
    const auto targetCount = coverage.targetCount();
    const auto componentCount = config.components.size();
    std::vector<std::size_t> alternativeColumnOffsets(componentCount + 1, 0);
    for (std::size_t component = 0; component < componentCount; ++component) {
        alternativeColumnOffsets[component + 1] = alternativeColumnOffsets[component];
        if (config.components[component].kind ==
            ComponentKind::AlternativeWitnesses) {
            alternativeColumnOffsets[component + 1] +=
                config.components[component].alternativeCount();
        }
    }
    const std::size_t alternativeCount = alternativeColumnOffsets.back();
    if (elementCount + targetCount + componentCount + alternativeCount >=
        static_cast<std::size_t>(std::numeric_limits<HighsInt>::max())) {
        throw std::invalid_argument("Model exceeds HiGHS index range");
    }

    HighsModel model;
    auto& lp = model.lp_;
    lp.num_col_ = static_cast<HighsInt>(
        elementCount + targetCount + componentCount + alternativeCount);
    lp.sense_ = ObjSense::kMaximize;
    lp.col_lower_.assign(lp.num_col_, 0.0);
    lp.col_upper_.assign(lp.num_col_, 1.0);
    lp.col_cost_.assign(lp.num_col_, 0.0);
    lp.integrality_.assign(lp.num_col_, HighsVarType::kInteger);
    for (std::size_t target = 0; target < targetCount; ++target) {
        lp.col_cost_[elementCount + target] = coverage.targetWeights[target];
    }

    std::vector<std::vector<HighsInt>> coveringElements(targetCount);
    for (std::size_t element = 0; element < elementCount; ++element) {
        for (const TargetId target :
             coverage.targets(static_cast<ElementId>(element))) {
            coveringElements[target].push_back(static_cast<HighsInt>(element));
        }
    }

    auto& a = lp.a_matrix_;
    a.format_ = MatrixFormat::kRowwise;
    a.start_.assign(1, 0);
    auto entry = [&](HighsInt col, double value) {
        a.index_.push_back(col);
        a.value_.push_back(value);
    };
    auto endRow = [&](double upper) {
        if (a.index_.size() >=
            static_cast<std::size_t>(std::numeric_limits<HighsInt>::max())) {
            throw std::invalid_argument("Model matrix exceeds HiGHS index range");
        }
        a.start_.push_back(static_cast<HighsInt>(a.index_.size()));
        lp.row_lower_.push_back(-kHighsInf);
        lp.row_upper_.push_back(upper);
    };

    for (std::size_t element = 0; element < elementCount; ++element) {
        entry(static_cast<HighsInt>(element), 1.0);
    }
    endRow(static_cast<double>(budget));

    for (std::size_t target = 0; target < targetCount; ++target) {
        entry(static_cast<HighsInt>(elementCount + target), 1.0);
        for (const HighsInt element : coveringElements[target]) {
            entry(element, -1.0);
        }
        endRow(0.0);
    }

    // Nonnegative rewards make one-sided activation constraints sufficient.
    for (std::size_t component = 0; component < componentCount; ++component) {
        const auto& entryConfig = config.components[component];
        const auto componentColumn = static_cast<HighsInt>(
            elementCount + targetCount + component);
        lp.col_cost_[componentColumn] = entryConfig.value;
        if (entryConfig.kind == ComponentKind::ThresholdCoverage) {
            entry(componentColumn, static_cast<double>(entryConfig.threshold));
            for (const HighsInt element : coveringElements[entryConfig.target]) {
                entry(element, -1.0);
            }
            endRow(0.0);
            continue;
        }

        // z_w may be one only when every element of witness w is selected.
        // Positive reward then permits g_i only when at least one z_w is one.
        entry(componentColumn, 1.0);
        for (std::size_t witness = 0;
             witness < entryConfig.alternativeCount(); ++witness) {
            const auto zColumn = static_cast<HighsInt>(
                elementCount + targetCount + componentCount +
                alternativeColumnOffsets[component] + witness);
            entry(zColumn, -1.0);
        }
        endRow(0.0);
        for (std::size_t witness = 0;
             witness < entryConfig.alternativeCount(); ++witness) {
            const auto zColumn = static_cast<HighsInt>(
                elementCount + targetCount + componentCount +
                alternativeColumnOffsets[component] + witness);
            for (ElementId element : entryConfig.alternative(witness)) {
                entry(zColumn, 1.0);
                entry(static_cast<HighsInt>(element), -1.0);
                endRow(0.0);
            }
        }
    }

    lp.num_row_ = static_cast<HighsInt>(lp.row_upper_.size());
    a.num_col_ = lp.num_col_;
    a.num_row_ = lp.num_row_;

    Highs solver;
    const auto check = [](HighsStatus status) {
        if (status == HighsStatus::kError) {
            throw std::runtime_error("HiGHS API error");
        }
    };
    check(solver.setOptionValue("output_flag", false));
    check(solver.setOptionValue("threads", 1));
    check(solver.setOptionValue("random_seed", 1));
    check(solver.setOptionValue("time_limit", timeLimitSeconds));
    check(solver.setOptionValue("mip_rel_gap", 0.0));
    check(solver.setOptionValue("mip_abs_gap", 0.0));
    check(solver.passModel(model));
    const auto built = Clock::now();
    check(solver.run());
    const auto finished = Clock::now();

    MipResult r;
    r.buildSeconds = std::chrono::duration<double>(built - start).count();
    r.solveSeconds = std::chrono::duration<double>(finished - built).count();
    r.status = solver.modelStatusToString(solver.getModelStatus());
    r.optimal = solver.getModelStatus() == HighsModelStatus::kOptimal;
    const auto& info = solver.getInfo();
    r.feasible = info.primal_solution_status == kSolutionStatusFeasible;
    r.upperBound = info.valid
        ? info.mip_dual_bound
        : std::numeric_limits<double>::quiet_NaN();
    r.gap = r.feasible
        ? info.mip_gap
        : std::numeric_limits<double>::quiet_NaN();
    if (r.feasible) {
        CoverageOracle oracle(
            std::make_shared<CoverageInstance>(coverage), config);
        const auto& solution = solver.getSolution();
        for (std::size_t element = 0; element < elementCount; ++element) {
            if (solution.col_value[element] > 0.5) {
                r.elements.push_back(static_cast<ElementId>(element));
                oracle.add(static_cast<ElementId>(element));
            }
        }
        if (r.elements.size() > budget) {
            throw std::runtime_error("MIP violated budget");
        }
        r.fValue = oracle.fValue();
        r.totalValue = oracle.totalValue();
        if (std::isfinite(r.upperBound)) {
            r.gap = r.totalValue != 0.0
                ? std::max(0.0, r.upperBound - r.totalValue) /
                    std::abs(r.totalValue)
                : (r.upperBound == 0.0
                       ? 0.0
                       : std::numeric_limits<double>::infinity());
        }
        if (r.optimal &&
            std::abs(r.totalValue - info.objective_function_value) >
                1e-6 * (1.0 + std::abs(r.totalValue))) {
            throw std::runtime_error(
                "MIP and coverage oracle disagree: budget=" +
                std::to_string(budget) + " oracle=" +
                std::to_string(r.totalValue) + " solver=" +
                std::to_string(info.objective_function_value));
        }
    }
    return r;
}

}  // namespace nsmsc
