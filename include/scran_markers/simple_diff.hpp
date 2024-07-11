#ifndef SCRAN_MARKERS_SIMPLE_DIFF_HPP
#define SCRAN_MARKERS_SIMPLE_DIFF_HPP

#include <limits>

namespace scran_markers {

namespace internal {

template<typename Stat_, typename Weight_>
Stat_ compute_pairwise_simple_diff(size_t g1, size_t g2, const Stat_* values, const Weight_* block_weights, size_t nblocks) {
    Stat_ total_weight = 0;
    Stat_ output = 0;

    for (size_t b = 0; b < nblocks; ++b) {
        size_t offset1 = g1 * nblocks + b; // no need to cast, everything's already a size_t.
        auto left = values[offset1];
        auto lweight = block_weights[offset1];
        if (!lweight) {
            continue;
        }

        size_t offset2 = g2 * nblocks + b;
        auto right = values[offset2]; 
        auto rweight = block_weights[offset2];
        if (!rweight) {
            continue;
        }

        Stat_ weight = static_cast<Stat_>(lweight) * static_cast<Stat_>(rweight);
        total_weight += weight;
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
void compute_pairwise_simple_diff(const Stat_* values, const Weight_* block_weights, size_t ngroups, size_t nblocks, Stat_* output) {
    for (size_t g1 = 0; g1 < ngroups; ++g1) {
        for (size_t g2 = 0; g2 < g1; ++g2) {
            auto d = compute_pairwise_simple_diff(g1, g2, values, block_weights, nblocks);
            output[g1 * ngroups + g2] = d;
            output[g2 * ngroups + g1] = -d;
        }
    }
}

}

}

#endif
