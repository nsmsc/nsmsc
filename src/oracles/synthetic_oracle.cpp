#include "nsmsc/synthetic_oracle.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace nsmsc {

SyntheticOracle::SyntheticOracle(
    std::vector<std::vector<std::size_t>> coveredTargets,
    std::vector<double> targetWeights,
    std::vector<ComplementaryComponent> components)
    : coveredTargets_(std::move(coveredTargets)),
      targetWeights_(std::move(targetWeights)),
      components_(std::move(components)),
      selected_(coveredTargets_.size(), false),
      targetCoverCounts_(targetWeights_.size(), 0) {
    for (const double weight : targetWeights_) {
        if (!std::isfinite(weight) || weight < 0.0) {
            throw std::invalid_argument(
                "Target weights must be finite and nonnegative");
        }
    }
    for (const auto& targets : coveredTargets_) {
        std::unordered_set<std::size_t> uniqueTargets;
        for (const std::size_t target : targets) {
            if (target >= targetWeights_.size()) {
                throw std::invalid_argument("Coverage target is out of range");
            }
            if (!uniqueTargets.insert(target).second) {
                throw std::invalid_argument(
                    "An element contains a duplicate coverage target");
            }
        }
    }
    for (const auto& component : components_) {
        if (!std::isfinite(component.value) || component.value < 0.0) {
            throw std::invalid_argument(
                "Component values must be finite and nonnegative");
        }
        if (component.witness.empty()) {
            throw std::invalid_argument("Component witnesses must be nonempty");
        }
        std::unordered_set<ElementId> uniqueElements;
        for (const ElementId element : component.witness) {
            validateElement(element);
            if (!uniqueElements.insert(element).second) {
                throw std::invalid_argument(
                    "A witness contains a duplicate element");
            }
        }
    }
}

std::size_t SyntheticOracle::groundSize() const {
    return coveredTargets_.size();
}

std::size_t SyntheticOracle::componentCount() const {
    return components_.size();
}

void SyntheticOracle::reset() {
    std::fill(selected_.begin(), selected_.end(), false);
    std::fill(targetCoverCounts_.begin(), targetCoverCounts_.end(), 0);
    selectedSize_ = 0;
    fValue_ = 0.0;
}

bool SyntheticOracle::selected(ElementId element) const {
    validateElement(element);
    return selected_[element];
}

std::size_t SyntheticOracle::selectedSize() const {
    return selectedSize_;
}

double SyntheticOracle::fValue() const {
    return fValue_;
}

double SyntheticOracle::fMarginal(ElementId element) const {
    validateElement(element);
    if (selected_[element]) {
        return 0.0;
    }
    double gain = 0.0;
    for (const std::size_t target : coveredTargets_[element]) {
        if (targetCoverCounts_[target] == 0) {
            gain += targetWeights_[target];
        }
    }
    return gain;
}

double SyntheticOracle::totalMarginal(ElementId element) const {
    validateElement(element);
    if (selected_[element]) {
        return 0.0;
    }
    double gain = fMarginal(element);
    for (const auto& component : components_) {
        bool containsCandidate = false;
        bool allOtherWitnessElementsSelected = true;
        for (const ElementId witnessElement : component.witness) {
            if (witnessElement == element) {
                containsCandidate = true;
            } else if (!selected_[witnessElement]) {
                allOtherWitnessElementsSelected = false;
            }
        }
        if (containsCandidate && allOtherWitnessElementsSelected) {
            gain += component.value;
        }
    }
    return gain;
}

void SyntheticOracle::add(ElementId element) {
    validateElement(element);
    if (selected_[element]) {
        throw std::invalid_argument("Cannot add an already selected element");
    }
    selected_[element] = true;
    ++selectedSize_;
    for (const std::size_t target : coveredTargets_[element]) {
        if (targetCoverCounts_[target]++ == 0) {
            fValue_ += targetWeights_[target];
        }
    }
}

void SyntheticOracle::remove(ElementId element) {
    validateElement(element);
    if (!selected_[element]) {
        throw std::invalid_argument("Cannot remove an unselected element");
    }
    selected_[element] = false;
    --selectedSize_;
    for (const auto target : coveredTargets_[element]) {
        if (--targetCoverCounts_[target] == 0) {
            fValue_ -= targetWeights_[target];
        }
    }
}

double SyntheticOracle::gValue(ComponentId component) const {
    validateComponent(component);
    const auto& entry = components_[component];
    const bool saturated = std::all_of(
        entry.witness.begin(), entry.witness.end(),
        [this](ElementId element) { return selected_[element]; });
    return saturated ? entry.value : 0.0;
}

double SyntheticOracle::gMaximum(ComponentId component) const {
    validateComponent(component);
    return components_[component].value;
}

std::span<const ElementId> SyntheticOracle::witness(
    ComponentId component) const {
    validateComponent(component);
    return components_[component].witness;
}

double SyntheticOracle::totalValue() const {
    double value = fValue_;
    for (std::size_t i = 0; i < components_.size(); ++i) {
        value += gValue(static_cast<ComponentId>(i));
    }
    return value;
}

void SyntheticOracle::validateElement(ElementId element) const {
    if (static_cast<std::size_t>(element) >= coveredTargets_.size()) {
        throw std::out_of_range("ElementId is out of range");
    }
}

void SyntheticOracle::validateComponent(ComponentId component) const {
    if (static_cast<std::size_t>(component) >= components_.size()) {
        throw std::out_of_range("ComponentId is out of range");
    }
}

}  // namespace nsmsc
