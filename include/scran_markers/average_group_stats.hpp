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

        for (size_t b = 0; b < nblocks; ++b) {
            // Remember, blocks are the slower changing dimension, so we need to jump by 'ngroups'.
            size_t offset = b * ngroups + g; // already size_t's.
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

template<typename Stat_>
void fill_average_results(
    size_t ngenes,
    size_t ngroups,
    std::vector<std::vector<Stat_> >& mean_res, 
    std::vector<std::vector<Stat_> >& detected_res, 
    std::vector<Stat_*>& mean_ptrs,
    std::vector<Stat_*>& detected_ptrs)
{
    mean_res.reserve(ngroups);
    detected_res.reserve(ngroups);
    mean_ptrs.reserve(ngroups);
    detected_ptrs.reserve(ngroups);
    for (size_t g = 0; g < ngroups; ++g) {
        mean_res.emplace_back(ngenes);
        detected_res.emplace_back(ngenes);
        mean_ptrs.emplace_back(mean_res.back().data());
        detected_ptrs.emplace_back(detected_res.back().data());
    }
}

}

}

#endif
