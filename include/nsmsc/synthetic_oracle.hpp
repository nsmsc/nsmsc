#pragma once

#include "nsmsc/objective_oracle.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace nsmsc {

struct ComplementaryComponent {
    double value = 0.0;
    std::vector<ElementId> witness;
};

// A weighted coverage function plus all-or-nothing complementary components.
class SyntheticOracle final : public ObjectiveOracle {
public:
    SyntheticOracle(
        std::vector<std::vector<std::size_t>> coveredTargets,
        std::vector<double> targetWeights,
        std::vector<ComplementaryComponent> components);

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

    [[nodiscard]] double totalValue() const override;

private:
    void validateElement(ElementId element) const;
    void validateComponent(ComponentId component) const;

    std::vector<std::vector<std::size_t>> coveredTargets_;
    std::vector<double> targetWeights_;
    std::vector<ComplementaryComponent> components_;

    std::vector<bool> selected_;
    std::vector<std::size_t> targetCoverCounts_;
    std::size_t selectedSize_ = 0;
    double fValue_ = 0.0;
};

}  // namespace nsmsc
