#include "nsmsc/coverage_oracle.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace nsmsc {

CoverageOracle::CoverageOracle(
    std::shared_ptr<const CoverageInstance> instance,
    ComplementarityConfig config)
    : instance_(std::move(instance)), config_(std::move(config)) {
    if (!instance_) {
        throw std::invalid_argument("Coverage instance cannot be null");
    }
    config_.validate(*instance_);
    componentsByTarget_.resize(instance_->targetCount());
    alternativesByElement_.resize(instance_->elementCount());
    alternativeSelectedCounts_.resize(config_.components.size());
    satisfiedAlternativeCounts_.assign(config_.components.size(), 0);
    for (std::size_t i = 0; i < config_.components.size(); ++i) {
        const auto& component = config_.components[i];
        if (component.kind == ComponentKind::ThresholdCoverage) {
            componentsByTarget_[component.target].push_back(i);
            continue;
        }
        alternativeSelectedCounts_[i].assign(component.alternativeCount(), 0);
        for (std::size_t j = 0; j < component.alternativeCount(); ++j) {
            for (ElementId element : component.alternative(j)) {
                alternativesByElement_[element].push_back(
                    {static_cast<ComponentId>(i), static_cast<std::uint32_t>(j)});
            }
        }
    }
    selected_.assign(instance_->elementCount(), 0);
    targetCoverCounts_.assign(instance_->targetCount(), 0);
}

std::size_t CoverageOracle::groundSize() const {
    return instance_->elementCount();
}

std::size_t CoverageOracle::componentCount() const {
    return config_.components.size();
}

void CoverageOracle::reset() {
    std::fill(selected_.begin(), selected_.end(), 0);
    std::fill(targetCoverCounts_.begin(), targetCoverCounts_.end(), 0);
    for (auto& counts : alternativeSelectedCounts_) {
        std::fill(counts.begin(), counts.end(), 0);
    }
    std::fill(satisfiedAlternativeCounts_.begin(),
              satisfiedAlternativeCounts_.end(), 0);
    selectedSize_ = 0;
    fValue_ = 0.0;
}

bool CoverageOracle::selected(ElementId element) const {
    validateElement(element);
    return selected_[element] != 0;
}

std::size_t CoverageOracle::selectedSize() const {
    return selectedSize_;
}

double CoverageOracle::fValue() const {
    return fValue_;
}

double CoverageOracle::fMarginal(ElementId element) const {
    validateElement(element);
    if (selected_[element] != 0) {
        return 0.0;
    }
    double gain = 0.0;
    for (const TargetId target : instance_->targets(element)) {
        if (targetCoverCounts_[target] == 0) {
            gain += instance_->targetWeights[target];
        }
    }
    return gain;
}

double CoverageOracle::totalMarginal(ElementId element) const {
    validateElement(element);
    if (selected_[element] != 0) {
        return 0.0;
    }
    double gain = fMarginal(element);
    for (const TargetId target : instance_->targets(element)) {
        for (const std::size_t componentIndex : componentsByTarget_[target]) {
            const auto& component = config_.components[componentIndex];
            if (targetCoverCounts_[target] < component.threshold &&
                targetCoverCounts_[target] + 1 >= component.threshold) {
                gain += component.value;
            }
        }
    }
    std::vector<std::uint8_t> gained(config_.components.size(), 0);
    for (const auto occurrence : alternativesByElement_[element]) {
        const auto componentIndex = static_cast<std::size_t>(occurrence.component);
        if (satisfiedAlternativeCounts_[componentIndex] != 0 ||
            gained[componentIndex] != 0) continue;
        const auto& component = config_.components[componentIndex];
        if (alternativeSelectedCounts_[componentIndex][occurrence.alternative] + 1 ==
            component.alternative(occurrence.alternative).size()) {
            gain += component.value;
            gained[componentIndex] = 1;
        }
    }
    return gain;
}

void CoverageOracle::add(ElementId element) {
    validateElement(element);
    if (selected_[element] != 0) {
        throw std::invalid_argument("Cannot add an already selected element");
    }
    selected_[element] = 1;
    ++selectedSize_;
    for (const TargetId target : instance_->targets(element)) {
        if (targetCoverCounts_[target]++ == 0) {
            fValue_ += instance_->targetWeights[target];
        }
    }
    for (const auto occurrence : alternativesByElement_[element]) {
        auto& count = alternativeSelectedCounts_[occurrence.component]
                                               [occurrence.alternative];
        ++count;
        if (count == config_.components[occurrence.component]
                         .alternative(occurrence.alternative).size()) {
            ++satisfiedAlternativeCounts_[occurrence.component];
        }
    }
}

void CoverageOracle::remove(ElementId element) {
    validateElement(element);
    if (!selected_[element]) {
        throw std::invalid_argument("Cannot remove an unselected element");
    }
    selected_[element] = false;
    --selectedSize_;
    for (const auto target : instance_->targets(element)) {
        if (--targetCoverCounts_[target] == 0) {
            fValue_ -= instance_->targetWeights[target];
        }
    }
    for (const auto occurrence : alternativesByElement_[element]) {
        auto& count = alternativeSelectedCounts_[occurrence.component]
                                               [occurrence.alternative];
        if (count == config_.components[occurrence.component]
                         .alternative(occurrence.alternative).size()) {
            --satisfiedAlternativeCounts_[occurrence.component];
        }
        --count;
    }
}

double CoverageOracle::gValue(ComponentId component) const {
    validateComponent(component);
    const auto& entry = config_.components[component];
    if (entry.kind == ComponentKind::AlternativeWitnesses) {
        return satisfiedAlternativeCounts_[component] != 0 ? entry.value : 0.0;
    }
    return targetCoverCounts_[entry.target] >= entry.threshold ? entry.value : 0.0;
}

double CoverageOracle::gMaximum(ComponentId component) const {
    validateComponent(component);
    return config_.components[component].value;
}

std::span<const ElementId> CoverageOracle::witness(
    ComponentId component) const {
    validateComponent(component);
    return config_.components[component].witness;
}

std::size_t CoverageOracle::witnessCount(ComponentId component) const {
    validateComponent(component);
    const auto& entry = config_.components[component];
    return entry.kind == ComponentKind::AlternativeWitnesses
        ? entry.alternativeCount()
        : 1;
}

std::span<const ElementId> CoverageOracle::witness(
    ComponentId component, std::size_t alternative) const {
    validateComponent(component);
    const auto& entry = config_.components[component];
    if (entry.kind == ComponentKind::AlternativeWitnesses) {
        return entry.alternative(alternative);
    }
    if (alternative != 0) throw std::out_of_range("Witness is out of range");
    return entry.witness;
}

double CoverageOracle::totalValue() const {
    double result = fValue_;
    for (std::size_t i = 0; i < componentCount(); ++i) {
        result += gValue(static_cast<ComponentId>(i));
    }
    return result;
}

const CoverageInstance& CoverageOracle::instance() const {
    return *instance_;
}

const ComplementarityConfig& CoverageOracle::config() const {
    return config_;
}

void CoverageOracle::validateElement(ElementId element) const {
    if (static_cast<std::size_t>(element) >= groundSize()) {
        throw std::out_of_range("ElementId is out of range");
    }
}

void CoverageOracle::validateComponent(ComponentId component) const {
    if (static_cast<std::size_t>(component) >= componentCount()) {
        throw std::out_of_range("ComponentId is out of range");
    }
}

}  // namespace nsmsc
