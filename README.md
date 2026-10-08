# GCC Experiment Code

C++20 implementation of Greedy with Complementary Completion (GCC), its
ablations, and baseline algorithms for coverage and influence maximization.

## Build

Requires CMake 3.20+ and a C++20 compiler. Release builds use `-O3 -DNDEBUG`
with GCC and Clang.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Executables: `run_experiment`, `run_coverage`, and `run_higgs_experiment`.

## Coverage Experiments

Provide a prepared coverage instance and its matching complementarity file:

```bash
METHOD_FILTER=proposed ./build/run_experiment \
  /path/to/coverage.bin /path/to/config.bin \
  100 /path/to/output.jsonl 1 \
  --augmentation=48 --omit-selected-elements
```

Syntax:

```text
run_experiment <coverage-input> <config-input> <k> <result-output>
               <seed|time> [threshold-epsilon] [ls-epsilon]
               [ftg-epsilon] [ftg-alpha] [--augmentation=<b>]
               [--omit-selected-elements]
```

`k` is the base budget; `--augmentation=b` sets the extra budget. Omitting
the option uses the theoretical augmentation. Use an integer seed for
reproducibility or `time` for a time-derived seed. Optional algorithm parameters
default to `0.1`, `0.1`, `0.1`, and `1.0`, respectively.

`METHOD_FILTER` selects one method; leaving it unset runs all coverage methods:

- GCC: `proposed`.
- Ablations: `proposed-original` (no budget filling),
  `proposed-no-enumeration`, `proposed-no-ranking`.
- Variants: `proposed-threshold`, `proposed-cost-aware`, `proposed-adaptive`.
- Baselines: `f_greedy`, `F_greedy`, `F-Greedy-LS`, `FTG-2022`,
  `threshold_greedy_2024`, `RRG`, `random`.

Minimal GCC-only command:

```bash
./build/run_coverage \
  /path/to/coverage.bin /path/to/config.bin 100 48 lazy
```

## Influence-Maximization Experiments

Use independent optimization and evaluation RR-set pools and a matching
component file. The executable supports datasets beyond Higgs:

```bash
METHOD_FILTER=proposed-HIST ./build/run_higgs_experiment \
  /path/to/optimization.rrbin \
  /path/to/evaluation.rrbin \
  /path/to/components.bin 100 /path/to/output.jsonl 1 \
  --dataset=Higgs --augmentation=128 \
  --hist-k=/path/to/hist_k.seed \
  --hist-aug=/path/to/hist_augmented.seed
```

HIST-based methods require external seed orders via `--hist-k` and
`--hist-aug`. OPIM-C uses `--opimc-k` and `--opimc-aug`. These files must
contain orders for the corresponding base and augmented budgets.

## Optional MIP Baseline

Build `run_mip` with an external HiGHS source tree:

```bash
cmake -S . -B build-highs \
  -DCMAKE_BUILD_TYPE=Release \
  -DNSMSC_WITH_HIGHS=ON \
  -DNSMSC_HIGHS_SOURCE_DIR=/path/to/HiGHS
cmake --build build-highs -j --target run_mip
```

For an installed HiGHS CMake package, omit `NSMSC_HIGHS_SOURCE_DIR`.

```bash
./build-highs/run_mip \
  /path/to/coverage.bin /path/to/config.bin \
  148 1800 /path/to/mip.jsonl
```

Here `148` is the total budget and `1800` is the solver time limit in seconds.

## Results

Experiment runners write JSONL records with objective values, budgets, selected
sizes, seeds, and runtimes. Coverage results use `total_value`; influence
results use `evaluation_total_value` from the evaluation RR pool.
`--omit-selected-elements` suppresses selected IDs in coverage results.

For comparisons, reuse the same input files and evaluation pool, fix random
seeds, and use equal allowed budgets. Coverage/configuration compatibility is
checked using a stored fingerprint.
