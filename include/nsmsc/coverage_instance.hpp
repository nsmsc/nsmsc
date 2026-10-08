#pragma once

#include "nsmsc/objective_oracle.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace nsmsc {

using TargetId = std::uint32_t;

struct CoverageInstance {
    std::vector<double> targetWeights;
    std::vector<std::uint32_t> elementCosts;
    std::vector<std::uint64_t> offsets;
    std::vector<TargetId> targetIds;

    [[nodiscard]] std::size_t targetCount() const;
    [[nodiscard]] std::size_t elementCount() const;
    [[nodiscard]] std::size_t incidenceCount() const;
    [[nodiscard]] std::span<const TargetId> targets(ElementId element) const;
    void validate() const;
};

enum class ComponentKind : std::uint32_t {
    ThresholdCoverage = 0,
    AlternativeWitnesses = 1,
};

struct ThresholdCoverageComponent {
    ThresholdCoverageComponent() = default;
    ThresholdCoverageComponent(TargetId targetValue, std::uint32_t thresholdValue,
                               double componentValue,
                               std::vector<ElementId> canonicalWitness)
        : target(targetValue), threshold(thresholdValue), value(componentValue),
          witness(std::move(canonicalWitness)) {}

    TargetId target = 0;
    std::uint32_t threshold = 0;
    double value = 0.0;
    std::vector<ElementId> witness;
    ComponentKind kind = ComponentKind::ThresholdCoverage;
    std::vector<std::uint64_t> alternativeOffsets;
    std::vector<ElementId> alternativeElements;

    [[nodiscard]] std::size_t alternativeCount() const;
    [[nodiscard]] std::span<const ElementId> alternative(
        std::size_t index) const;
};

struct ComplementarityConfig {
    std::uint64_t coverageFingerprint = 0;
    std::uint64_t seed = 0;
    std::size_t tau = 0;
    std::vector<ThresholdCoverageComponent> components;

    void validate(const CoverageInstance& instance) const;
};

[[nodiscard]] std::uint64_t instanceFingerprint(
    const CoverageInstance& instance);

}  // namespace nsmsc
