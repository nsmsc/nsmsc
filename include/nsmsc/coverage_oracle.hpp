#pragma once

#include "nsmsc/coverage_instance.hpp"
#include "nsmsc/objective_oracle.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <vector>

namespace nsmsc {

class CoverageOracle final : public ObjectiveOracle {
public:
    CoverageOracle(
        std::shared_ptr<const CoverageInstance> instance,
        ComplementarityConfig config);

    [[nodiscard]] std::size_t groundSize() const override;
    [[nodiscard]] std::size_t componentCount() const override;

    void reset() override;
    [[nodiscard]] bool selected(ElementId element) const override;
    [[nodiscard]] std::size_t selectedSize() const override;

    [[nodiscard]] double fValue() const override;
    [[nodiscard]] double fMarginal(ElementId element) const override;
    [[nodiscard]] double totalMarginal(ElementId element) const override;
    void add(ElementId element) override;
    void remove(ElementId element) override;

    [[nodiscard]] double gValue(ComponentId component) const override;
    [[nodiscard]] double gMaximum(ComponentId component) const override;
    [[nodiscard]] std::span<const ElementId> witness(
        ComponentId component) const override;
    [[nodiscard]] std::size_t witnessCount(
        ComponentId component) const override;
    [[nodiscard]] std::span<const ElementId> witness(
        ComponentId component, std::size_t alternative) const override;

    [[nodiscard]] double totalValue() const override;

    [[nodiscard]] const CoverageInstance& instance() const;
    [[nodiscard]] const ComplementarityConfig& config() const;

private:
    void validateElement(ElementId element) const;
    void validateComponent(ComponentId component) const;

    std::shared_ptr<const CoverageInstance> instance_;
    ComplementarityConfig config_;
    std::vector<std::vector<std::size_t>> componentsByTarget_;
    struct AlternativeOccurrence {
        ComponentId component = 0;
        std::uint32_t alternative = 0;
    };
    std::vector<std::vector<AlternativeOccurrence>> alternativesByElement_;
    std::vector<std::vector<std::uint32_t>> alternativeSelectedCounts_;
    std::vector<std::uint32_t> satisfiedAlternativeCounts_;
    std::vector<std::uint8_t> selected_;
    std::vector<std::uint32_t> targetCoverCounts_;
    std::size_t selectedSize_ = 0;
    double fValue_ = 0.0;
};

}  // namespace nsmsc
