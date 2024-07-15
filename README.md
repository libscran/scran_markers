# Marker detection for groups of cells

![Unit tests](https://github.com/libscran/scran_markers/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/libscran/scran_markers/actions/workflows/doxygenate.yaml/badge.svg)
[![Codecov](https://codecov.io/gh/libscran/scran_markers/graph/badge.svg?token=iL6GuHkCjz)](https://codecov.io/gh/libscran/scran_markers)

## Overview

This library contains functions for detecting group-specific markers (e.g., for clusters or cell types) from a single-cell expression matrix.
It performs differential analyses between pairs of groups, computing a variety of effect sizes like Cohen's d and the AUC.
The effect sizes are then summarized for each group to obtain some rankings for prioritizing interesting genes.
The code itself was originally derived from the [**scran** R package](https://bioconductor.org/packages/scran),
factored out into a separate C++ library for easier re-use.

## Quick start

Given a [`tatami::Matrix`](https://github.com/tatami-inc/tatami) and an array of group assignments,
the `score_markers_summary()` function will compute the aggregate statistics across all genes for each group.

```cpp
#include "scran_markers/scran_markers.hpp"

// Expression matrix, usually log-normalized.
const tatami::Matrix<double, int>& matrix = some_data_source();

// Array containing integer assignments to groups 0, 1, 2, etc.
std::vector<int> groupings = some_groupings();

scran_markers::ScoreMarkersSummaryOptions opt;
auto res = scran::score_markers_summary(matrix, groupings.data(), opt);

res.mean[0]; // mean of each gene in the first group.
res.detected[0]; // detected proportion of each gene in the first group.
```

The most interesting part of this result is the effect size summary for each group.
For each group, we compute the Cohen's d, AUC, delta-mean and delta-detected by comparing that group against every other group.
We then combine the effect sizes from all comparisons into some summary statistics like the mean or median.
Ranking by these summaries yields a list of potential group-specific marker genes for further examination.
Picking a different summary statistic and/or effect size will favor different types of markers in the ranking.

```cpp
res.cohens_d[0].mean; // mean Cohen's d across all genes for the first group
res.auc[0].median; // median AUC across all genes for the first group.
```

If the dataset contains some uninteresting variation (e.g., batches, samples),
we can ensure that it does not affect the effect size calculation by blocking on that factor.
This performs the comparisons within each level of the blocking factor so as to ignore the irrelevant variation.

```cpp
// Array containing integer assignments to blocks 0, 1, 2, etc.
std::vector<int> blocks = some_blocks();

auto block_res = scran::score_markers_summary_blocked(
    matrix,
    groupings.data(), 
    blocks.data(),
    opt
);
```

If more detail is necessary, we can obtain effect sizes from all pairwise comparisons using the `score_markers_pairwise()` function.

```cpp
scran_markers::ScoreMarkersPairwiseOptions popt;

auto pair_res = scran::score_markers_pairwise(
    matrix,
    groupings.data(),
    popt
);

// 3D arrays of effect sizes for each pairwise comparison and gene.
pair_res.cohens_d;
pair_res.auc; 
```

Alternatively, if we already have an array effect sizes, we can use the `summarize_effects()` function to obtain summaries for each group.
In fact, `score_markers_summary()` is just a more memory-efficient version of `score_markers_pairwise()` followed by `summarize_effects()`.

```cpp
auto cohen_summary = scran_markers::summarize_effects(
    matrix.nrows(),
    pair_res.mean.size(), // i.e., number of groups
    pair_res.cohens_d.data()
);
```

Check out the [reference documentation](https://libscran.github.io/scran_markers) for more details.

## Building projects

### CMake with `FetchContent`

If you're using CMake, you just need to add something like this to your `CMakeLists.txt`:

```cmake
include(FetchContent)

FetchContent_Declare(
  scran_markers
  GIT_REPOSITORY https://github.com/libscran/scran_markers
  GIT_TAG master # or any version of interest
)

FetchContent_MakeAvailable(scran_markers)
```

Then you can link to **scran_markers** to make the headers available during compilation:

```cmake
# For executables:
target_link_libraries(myexe libscran::scran_markers)

# For libaries
target_link_libraries(mylib INTERFACE libscran::scran_markers)
```

### CMake with `find_package()`

```cmake
find_package(libscran_scran_markers CONFIG REQUIRED)
target_link_libraries(mylib INTERFACE libscran::scran_markers)
```

To install the library, use:

```sh
mkdir build && cd build
cmake .. -DSCRAN_MARKERS_TESTS=OFF
cmake --build . --target install
```

By default, this will use `FetchContent` to fetch all external dependencies.
If you want to install them manually, use `-DSCRAN_MARKERS_FETCH_EXTERN=OFF`.
See the commit hashes in [`extern/CMakeLists.txt`](extern/CMakeLists.txt) to find compatible versions of each dependency.

### Manual

If you're not using CMake, the simple approach is to just copy the files in `include/` - either directly or with Git submodules - and include their path during compilation with, e.g., GCC's `-I`.
This requires the external dependencies listed in [`extern/CMakeLists.txt`](extern/CMakeLists.txt), which also need to be made available during compilation.
