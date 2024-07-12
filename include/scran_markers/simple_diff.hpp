#ifndef SCRAN_MARKERS_SIMPLE_DIFF_HPP
#define SCRAN_MARKERS_SIMPLE_DIFF_HPP

#include <limits>

namespace scran_markers {

namespace internal {

// 'values' and 'combo_weights' are expected to be 'ngroups * nblocks' arrays
// where groups are the faster-changing dimension and the blocks are slower.
template<typename Stat_, typename Weight_>
Stat_ compute_pairwise_simple_diff(size_t g1, size_t g2, const Stat_* values, size_t ngroups, size_t nblocks, const Weight_* combo_weights) {
    Stat_ total_weight = 0;
    Stat_ output = 0;

    size_t offset1 = g1, offset2 = g2; 
    for (size_t b = 0; b < nblocks; ++b, offset1 += ngroups, offset2 += ngroups) {
        auto lweight = combo_weights[offset1];
        if (!lweight) {
            continue;
        }

        auto rweight = combo_weights[offset2];
        if (!rweight) {
            continue;
        }

        Stat_ weight = static_cast<Stat_>(lweight) * static_cast<Stat_>(rweight);
        total_weight += weight;

        auto left = values[offset1];
        auto right = values[offset2]; 
        output += (left - right) * weight;
    }

    if (total_weight) {
        output /= total_weight;
    } else {
        output = std::numeric_limits<Stat_>::quiet_NaN();
    }

    return output;
}

template<typename Stat_, typename Weight_>
void compute_pairwise_simple_diff(const Stat_* values, size_t ngroups, size_t nblocks, const Weight_* combo_weights, Stat_* output) {
    for (size_t g1 = 0; g1 < ngroups; ++g1) {
        for (size_t g2 = 0; g2 < g1; ++g2) {
            auto d = compute_pairwise_simple_diff(g1, g2, values, ngroups, nblocks, combo_weights);
            output[g1 * ngroups + g2] = d;
            output[g2 * ngroups + g1] = -d;
        }
    }
}

}

}

#endif
