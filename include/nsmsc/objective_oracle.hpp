#pragma once

#include <cstddef>
#include <cstdint>
#include <span>

namespace nsmsc {

using ElementId = std::uint32_t;
using ComponentId = std::uint32_t;

class ObjectiveOracle {
public:
    virtual ~ObjectiveOracle() = default;

    [[nodiscard]] virtual std::size_t groundSize() const = 0;
    [[nodiscard]] virtual std::size_t componentCount() const = 0;

    virtual void reset() = 0;
    [[nodiscard]] virtual bool selected(ElementId element) const = 0;
    [[nodiscard]] virtual std::size_t selectedSize() const = 0;

    [[nodiscard]] virtual double fValue() const = 0;
    [[nodiscard]] virtual double fMarginal(ElementId element) const = 0;
    [[nodiscard]] virtual double totalMarginal(ElementId element) const = 0;
    virtual void add(ElementId element) = 0;
    virtual void remove(ElementId element) = 0;

    [[nodiscard]] virtual double gValue(ComponentId component) const = 0;
    [[nodiscard]] virtual double gMaximum(ComponentId component) const = 0;
    [[nodiscard]] virtual std::span<const ElementId> witness(
        ComponentId component) const = 0;
    [[nodiscard]] virtual std::size_t witnessCount(
        ComponentId component) const {
        static_cast<void>(component);
        return 1;
    }
    [[nodiscard]] virtual std::span<const ElementId> witness(
        ComponentId component, std::size_t alternative) const {
        if (alternative != 0) return {};
        return witness(component);
    }

    [[nodiscard]] virtual double totalValue() const = 0;
};

}  // namespace nsmsc
