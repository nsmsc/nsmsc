#pragma once

#include "nsmsc/objective_oracle.hpp"
#include "nsmsc/synthetic_oracle.hpp"

#include <cstdint>
#include <filesystem>
#include <random>
#include <span>
#include <vector>

namespace nsmsc {

struct InfluenceGraph {
    std::vector<std::uint64_t> externalIds;
    std::vector<std::vector<ElementId>> incoming;
};

struct RRSetCollection {
    std::size_t nodeCount = 0;
    std::uint64_t seed = 0;
    std::vector<std::vector<ElementId>> sets;
};

struct RRComplementarityConfig {
    std::size_t nodeCount = 0;
    std::size_t tau = 0;
    double rho = 0.0;
    std::vector<ComplementaryComponent> components;
};

InfluenceGraph readHiggsRetweetGraph(const std::filesystem::path& path);

class WeightedCascadeSampler {
public:
    WeightedCascadeSampler(
        const InfluenceGraph& graph,
        std::uint64_t seed,
        bool subsimJumpSampling = true);
    void append(RRSetCollection& collection, std::size_t count);

private:
    const InfluenceGraph& graph_;
    std::mt19937_64 random_;
    std::vector<std::uint32_t> marks_;
    std::uint32_t generation_ = 0;
    bool subsimJumpSampling_ = true;
};

void writeRRSetsBinary(
    const std::filesystem::path& path, const RRSetCollection& collection);
RRSetCollection readRRSetsBinary(const std::filesystem::path& path);
void writeRRComplementarityBinary(
    const std::filesystem::path& path, const RRComplementarityConfig& config);
RRComplementarityConfig readRRComplementarityBinary(
    const std::filesystem::path& path, std::size_t expectedNodeCount);

// Influence spread is estimated as n times the fraction of covered RR sets.
class RRSetOracle final : public ObjectiveOracle {
public:
    RRSetOracle(RRSetCollection collection,
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
    [[nodiscard]] std::span<const ElementId> witness(ComponentId component) const override;
    [[nodiscard]] double totalValue() const override;

private:
    void validateElement(ElementId element) const;
    RRSetCollection collection_;
    std::vector<ComplementaryComponent> components_;
    std::vector<std::vector<std::uint32_t>> setsByElement_;
    std::vector<std::vector<ComponentId>> componentsByElement_;
    std::vector<std::uint32_t> coverCount_;
    std::vector<std::uint8_t> selected_;
    std::size_t selectedSize_ = 0;
    std::size_t coveredSets_ = 0;
};

}  // namespace nsmsc
