#ifndef SCRAN_MARKERS_PRECOMPUTED_PAIRWISE_WEIGHTS_HPP
#define SCRAN_MARKERS_PRECOMPUTED_PAIRWISE_WEIGHTS_HPP

#include <vector>

namespace scran_markers {

namespace internal {

template<typename Weight_>
class PrecomputedPairwiseWeights {
public:
    // 'combo_weights' are expected to be 'ngroups * nblocks' arrays where
    // groups are the faster-changing dimension and the blocks are slower.
    PrecomputedPairwiseWeights(size_t ngroups, size_t nblocks, const Weight_* combo_weights) :
        my_total(ngroups * ngroups),
        my_by_block(my_total.size() * nblocks),
        my_ngroups(ngroups),
        my_nblocks(nblocks)
    {
        for (size_t b = 0; b < nblocks; ++b) {
            for (size_t g1 = 1; g1 < ngroups; ++g1) {
                auto w1 = combo_weights[b * ngroups + g1 /* already size_t's */];
                for (size_t g2 = 0; g2 < g1; ++g2) {
                    Weight_ combined = w1 * combo_weights[b * ngroups + g2 /* already size_t's */];

                    // Storing it as a 3D array where the blocks are the fastest changing, 
                    // and then the two groups are the next fastest changing.
                    size_t out_offset1 = g1 * ngroups + g2;
                    my_by_block[out_offset1 * nblocks + b] = combined;
                    my_by_block[(g2 * ngroups + g1) * nblocks + b] = combined;

                    my_total[out_offset1] += combined;
                }
            }
        }

        // Filling the other side, for completeness.
        for (size_t g1 = 1; g1 < ngroups; ++g1) {
            for (size_t g2 = 0; g2 < g1; ++g2) {
                my_total[g2 * ngroups + g1] = my_total[g1 * ngroups + g2];
            }
        }
    }

public:
    std::pair<const Weight_*, Weight_> get(size_t g1, size_t g2) const {
        size_t offset = g1 * my_ngroups + g2; // no need to cast.
        return std::make_pair(
            my_by_block.data() + offset * my_nblocks,
            my_total[offset]
        );
    }

private:
    std::vector<Weight_> my_total;
    std::vector<Weight_> my_by_block;
    size_t my_ngroups, my_nblocks;
};

}

}

#endif
