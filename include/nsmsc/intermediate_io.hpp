#pragma once

#include "nsmsc/coverage_instance.hpp"

#include <filesystem>

namespace nsmsc {

void writeCoverageBinary(
    const std::filesystem::path& path,
    const CoverageInstance& instance);

[[nodiscard]] CoverageInstance readCoverageBinary(
    const std::filesystem::path& path);

void writeComplementarityBinary(
    const std::filesystem::path& path,
    const ComplementarityConfig& config,
    const CoverageInstance& instance);

[[nodiscard]] ComplementarityConfig readComplementarityBinary(
    const std::filesystem::path& path,
    const CoverageInstance& instance);

}  // namespace nsmsc
