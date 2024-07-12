#ifndef SCRAN_MARKERS_AVERAGE_GROUP_STATS_HPP
#define SCRAN_MARKERS_AVERAGE_GROUP_STATS_HPP

#include <vector>
#include <limits>

namespace scran_markers {

namespace internal {

template<typename Weight_>
std::vector<Weight_> compute_total_weight_per_group(size_t ngroups, size_t nblocks, const Weight_* combo_weights) {
    std::vector<Weight_> output(ngroups);
    for (size_t b = 0; b < nblocks; ++b) {
        for (size_t g = 0; g < ngroups; ++g) {
            output[g] += *combo_weights;
            ++combo_weights;
        }
    }
    return output;
}

template<typename Stat_, typename Weight_>
void average_group_stats(
    size_t gene, 
    size_t ngroups,
    size_t nblocks,
    const Stat_* tmp_means,
    const Stat_* tmp_detected,
    const Weight_* combo_weights,
    const Weight_* total_weights,
    const std::vector<Stat_*>& means,
    const std::vector<Stat_*>& detected)
{
    for (size_t g = 0; g < ngroups; ++g) {
        auto& gmean = means[g][gene];
        auto& gdet = detected[g][gene];

        auto total_weight = total_weights[g];
        if (total_weight == 0) {
            gdet = std::numeric_limits<Stat_>::quiet_NaN();
            gmean = std::numeric_limits<Stat_>::quiet_NaN();
            continue;
        }

        gmean = 0;
        gdet = 0;

        size_t offset = g;
        for (size_t b = 0; b < nblocks; ++b, offset += ngroups) {
            const auto& curweight = combo_weights[offset];
            if (curweight) {
                gmean += curweight * tmp_means[offset];
                gdet += curweight * tmp_detected[offset];
            } 
        }

        gmean /= total_weight;
        gdet /= total_weight;
    }
}

}

}

#endif
