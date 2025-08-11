#ifndef SCRAN_MARKERS_CREATE_COMBINATIONS_HPP
#define SCRAN_MARKERS_CREATE_COMBINATIONS_HPP

#include <vector>
#include <cstddef>

#include "sanisizer/sanisizer.hpp"

namespace scran_markers {

namespace internal {

// When we combine 'group' and 'block' into a single 'combinations' factor, the
// resulting combinations can be considered to index into a 2-dimensional array
// of dimension 'ngroups * nblocks' where the group is the faster-changing
// dimension. This 2D array layout is used for all 'combo_*'-prefixed arrays
// like 'combo_weights', 'combo_means', etc.
template<typename Group_, typename Block_, typename Index_>
std::vector<std::size_t> create_combinations(std::size_t ngroups, const Group_* group, const Block_* block, Index_ NC) {
    auto combinations = sanisizer::create<std::vector<std::size_t> >(NC);
    for (Index_ c = 0; c < NC; ++c) {
        combinations[c] = sanisizer::nd_offset<std::size_t>(group[c], ngroups, block[c]); // group is the faster changing dimension.
    }
    return combinations;
}

// We can't just use tatami_stats::tabulate_groups as downstream is expecting a 'ngroups * nblocks' array;
// tabulate_groups() will not report the full length if not all combinations are observed.
template<typename Index_>
std::vector<Index_> tabulate_combinations(std::size_t ngroups, std::size_t nblocks, const std::vector<std::size_t>& combinations) {
    std::vector<Index_> output(sanisizer::product<typename std::vector<Index_>::size_type>(ngroups, nblocks));
    for (auto c : combinations) {
        ++output[c];
    }
    return output;
}

}

}

#endif
