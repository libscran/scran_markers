#ifndef SCRAN_MARKERS_CREATE_COMBINATIONS_HPP
#define SCRAN_MARKERS_CREATE_COMBINATIONS_HPP

#include <vector>

namespace scran_markers {

namespace internal {

// When we combine 'group' and 'block' into a single 'combinations' factor, the
// resulting combinations can be considered to index into a 2-dimensional array
// of dimension 'ngroups * nblocks' where the group is the faster-changing
// dimension. This 2D array layout is used for all 'combo_*'-prefixed arrays
// like 'combo_weights', 'combo_means', etc.
template<typename Group_, typename Block_, typename Index_>
std::vector<size_t> create_combinations(size_t ngroups, const Group_* group, const Block_* block, Index_ NC) {
    std::vector<size_t> combinations(NC);
    for (Index_ c = 0; c < NC; ++c) {
        combinations[c] = static_cast<size_t>(block[c]) * ngroups + static_cast<size_t>(group[c]); // group is the faster changing dimension.
    }
    return combinations;
}

}

}

#endif
