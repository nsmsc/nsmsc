cmake_minimum_required(VERSION 3.20)

project(gcc_experiments LANGUAGES CXX)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING
        "Flags used by the C++ compiler for Release builds" FORCE)
endif()

add_library(nsmsc_core
    src/core/baselines.cpp
    src/core/complementary_completion.cpp
    src/core/experiment_protocol.cpp
    src/core/greedy.cpp
    src/data/coverage_instance.cpp
    src/data/intermediate_io.cpp
    src/oracles/coverage_oracle.cpp
    src/oracles/rr_set_oracle.cpp
    src/oracles/synthetic_oracle.cpp
)

target_compile_features(nsmsc_core PUBLIC cxx_std_20)
target_include_directories(nsmsc_core PUBLIC include)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(nsmsc_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(run_coverage src/tools/run_coverage.cpp)
target_link_libraries(run_coverage PRIVATE nsmsc_core)

add_executable(run_experiment src/tools/run_experiment.cpp)
target_link_libraries(run_experiment PRIVATE nsmsc_core)

add_executable(run_higgs_experiment src/tools/run_higgs_experiment.cpp)
target_link_libraries(run_higgs_experiment PRIVATE nsmsc_core)

option(NSMSC_WITH_HIGHS "Build the optional MIP-HiGHS baseline" OFF)
set(NSMSC_HIGHS_SOURCE_DIR "" CACHE PATH
    "Path to an external HiGHS source tree")

if(NSMSC_WITH_HIGHS)
    if(NSMSC_HIGHS_SOURCE_DIR)
        if(NOT EXISTS "${NSMSC_HIGHS_SOURCE_DIR}/CMakeLists.txt")
            message(FATAL_ERROR
                "NSMSC_HIGHS_SOURCE_DIR does not contain CMakeLists.txt")
        endif()
        set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
        set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(BUILD_CXX_EXE OFF CACHE BOOL "" FORCE)
        set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
        add_subdirectory(
            "${NSMSC_HIGHS_SOURCE_DIR}"
            "${CMAKE_BINARY_DIR}/_deps/highs-build"
            EXCLUDE_FROM_ALL)
    else()
        find_package(highs CONFIG REQUIRED)
    endif()

    if(TARGET highs)
        set(NSMSC_HIGHS_TARGET highs)
    elseif(TARGET highs::highs)
        set(NSMSC_HIGHS_TARGET highs::highs)
    else()
        message(FATAL_ERROR "A HiGHS CMake target was not found")
    endif()

    add_library(nsmsc_mip src/core/mip.cpp)
    target_compile_features(nsmsc_mip PUBLIC cxx_std_20)
    target_link_libraries(nsmsc_mip PUBLIC nsmsc_core ${NSMSC_HIGHS_TARGET})

    add_executable(run_mip src/tools/run_mip.cpp)
    target_link_libraries(run_mip PRIVATE nsmsc_mip)

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(nsmsc_mip PRIVATE -Wall -Wextra -Wpedantic)
        target_compile_options(run_mip PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()
