#ifndef SCRAN_MARKERS_AVERAGE_GROUP_STATS_HPP
#define SCRAN_MARKERS_AVERAGE_GROUP_STATS_HPP

#include <vector>
#include <limits>
#include <cstddef>

#include "sanisizer/sanisizer.hpp"

#include "utils.hpp"

namespace scran_markers {

namespace internal {

template<typename Weight_>
std::vector<Weight_> compute_total_weight_per_group(const std::size_t ngroups, const std::size_t nblocks, const Weight_* const combo_weights) {
    auto output = sanisizer::create<std::vector<Weight_> >(ngroups);
    for (decltype(I(nblocks)) b = 0; b < nblocks; ++b) {
        for (decltype(I(ngroups)) g = 0; g < ngroups; ++g) {
            output[g] += combo_weights[sanisizer::nd_offset<std::size_t>(g, ngroups, b)];
        }
    }
    return output;
}

template<typename Gene_, typename Stat_, typename Weight_>
void average_group_stats_blockmean(
    const Gene_ gene,
    const std::size_t ngroups,
    const std::size_t nblocks,
    const Stat_* const stats,
    const Weight_* const combo_weights,
    const Weight_* const total_weights,
    const std::vector<Stat_*>& out_stats 
) {
    for (decltype(I(ngroups)) g = 0; g < ngroups; ++g) {
        auto& output = out_stats[g][gene];

        const auto total_weight = total_weights[g];
        if (total_weight == 0) {
            output = std::numeric_limits<Stat_>::quiet_NaN();
            continue;
        }

        Stat_ sum = 0;
        for (decltype(I(nblocks)) b = 0; b < nblocks; ++b) {
            // Remember, blocks are the slower changing dimension, so we need to jump by 'ngroups'.
            const auto offset = sanisizer::nd_offset<std::size_t>(g, ngroups, b);
            const auto& curweight = combo_weights[offset];
            if (curweight) { // check if this is zero and skip it explicitly, as the value would probably be NaN. 
                sum += curweight * tmp_stats[offset];
            }
        }

        output = sum / total_weight;
    }
}

template<typename Gene_, typename Stat_, typename Weight_>
void average_group_stats_blockquantile(
    const Gene_ gene,
    const std::size_t ngroups,
    const std::size_t nblocks,
    const Stat_* const tmp_stats,
    std::vector<Stat_>& buffer,
    QuantileCalculator<Stat_>& qcalc,
    const std::vector<Stat_*>& out_stats 
) {
    for (decltype(I(ngroups)) g = 0; g < ngroups; ++g) {
        buffer.clear();

        for (decltype(I(nblocks)) b = 0; b < nblocks; ++b) {
            // Remember, blocks are the slower changing dimension, so we need to jump by 'ngroups'.
            const auto offset = sanisizer::nd_offset<std::size_t>(g, ngroups, b);
            const auto val = stats[offset];
            if (!std::isnan(val)) {
                buffer.push_back(val);
            }
        }

        out_stats[g][gene] = qcalc.compute(buffer);
    }
}

template<typename Gene_, typename Stat_>
void preallocate_average_results(
    const Gene_ ngenes,
    const std::size_t ngroups,
    std::vector<std::vector<Stat_> >& res, 
    std::vector<Stat_*>& ptrs
) {
    res.reserve(ngroups);
    ptrs.reserve(ngroups);
    for (decltype(I(ngroups)) g = 0; g < ngroups; ++g) {
        res.emplace_back(
            sanisizer::cast<decltype(I(res.front().size()))>(ngenes)
#ifdef SCRAN_MARKERS_TEST_INIT
            , SCRAN_MARKERS_TEST_INIT
#endif
        );
        ptrs.emplace_back(res.back().data());
    }
}

}

}

#endif
