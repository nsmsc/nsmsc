#include "nsmsc/experiment_protocol.hpp"

#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace nsmsc {

std::size_t theoreticalAugmentation(std::size_t h, std::size_t tau) {
    if (tau == 0) {
        throw std::invalid_argument("tau must be positive");
    }
    const double fraction = 1.0 - 1.0 / std::numbers::e;
    const auto completedComponents = static_cast<std::size_t>(
        std::ceil(fraction * static_cast<double>(h)));
    if (completedComponents != 0 &&
        tau > std::numeric_limits<std::size_t>::max() / completedComponents) {
        throw std::overflow_error("Theoretical augmentation overflows size_t");
    }
    return tau * completedComponents;
}

}  // namespace nsmsc
