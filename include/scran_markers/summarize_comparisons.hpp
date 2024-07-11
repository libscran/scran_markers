#ifndef SCRAN_MARKERS_SUMMARIZE_COMPARISONS_HPP
#define SCRAN_MARKERS_SUMMARIZE_COMPARISONS_HPP

#include "tatami_stats/tatami_stats.hpp"

#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>

/**
 * @file summarize_comparisons.hpp
 * @brief Utilities for effect summarization.
 */

namespace scran_markers {

/**
 * Choice of summary statistic for the pairwise effects of a given group, see `SummarizeEffects` for details.
 * `n_summaries` is used to denote the number of possible summary statistics. 
 */
enum class Summary : unsigned char {
    MIN,
    MEAN,
    MEDIAN,
    MAX,
    MIN_RANK,
    n_summaries 
};

/**
 * @cond
 */
namespace internal {

template<typename Stat_>
void summarize_comparisons_internal(size_t ngroups, const Stat_* effects, size_t group, size_t gene, std::vector<std::vector<Stat_*> >& output, std::vector<Stat_>& buffer) {
    auto ebegin = buffer.data();
    auto elast = ebegin;	

    // Ignoring the self comparison and pruning out NaNs.
    {
        auto eptr = effects;
        for (size_t r = 0; r < ngroups; ++r, ++eptr) {
            if (r == group || std::isnan(*eptr)) {
                continue;
            }
            *elast = *eptr;
            ++elast;
        }
    }

    size_t ncomps = elast - ebegin;
    constexpr size_t before_minrank = static_cast<size_t>(Summary::MIN_RANK);
    if (ncomps == 0) {
        for (size_t i = 0; i < before_minrank; ++i) {
            if (output[i].size()) {
                output[i][group][gene] = std::numeric_limits<Stat_>::quiet_NaN();
            }
        }
    } else if (ncomps == 1) {
        for (size_t i = 0; i < before_minrank; ++i) { 
            if (output[i].size()) {
                output[i][group][gene] = *ebegin;
            }
        }
    } else {
        if (output[static_cast<size_t>(Summary::MIN)].size()) {
            output[static_cast<size_t>(Summary::MIN)][group][gene] = *std::min_element(ebegin, elast);
        }
        if (output[static_cast<size_t>(Summary::MEAN)].size()) {
            output[static_cast<size_t>(Summary::MEAN)][group][gene] = std::accumulate(ebegin, elast, 0.0) / ncomps; 
        }
        if (output[static_cast<size_t>(Summary::MAX)].size()) {
            output[static_cast<size_t>(Summary::MAX)][group][gene] = *std::max_element(ebegin, elast);
        }
        if (output[static_cast<size_t>(Summary::MEDIAN)].size()) { // this mutates the buffer, so we put this last to avoid surprises.
            output[static_cast<size_t>(Summary::MEDIAN)][group][gene] = tatami_stats::medians::direct(ebegin, ncomps, /* skip_nan = */ false); 
        }
    }
}

template<typename Stat_>
void summarize_comparisons(size_t ngenes, size_t ngroups, const Stat_* effects, std::vector<std::vector<Stat_*> >& output, int threads) {
    size_t shift = ngroups * ngroups;

    tatami::parallelize([&](size_t, size_t start, size_t length) -> void {
        std::vector<Stat_> buffer(ngroups);
        auto effect_ptr = effects + start * shift; // everything is already a size_t, so it's fine.

        for (size_t gene = start, end = start + length; gene < end; ++gene, effect_ptr += shift) {
            auto current_effects = effect_ptr;
            for (size_t l = 0; l < ngroups; ++l, current_effects += ngroups) {
                summarize_comparisons_internal(ngroups, current_effects, l, gene, output, buffer);
            }
        }
    }, ngenes, threads);
}

template<typename Stat_, typename Index_>
Index_ fill_and_sort_rank_buffer(const Stat_* effects, size_t stride, std::vector<std::pair<Stat_, Index_> >& buffer) {
    auto bIt = buffer.begin();
    for (Index_ i = 0, end = buffer.size(); i < end; ++i, effects += stride) {
        if (!std::isnan(*effects)) {
            bIt->first = -*effects; // negative to sort by decreasing value.
            bIt->second = i;
            ++bIt;
        }
    }
    std::sort(buffer.begin(), bIt);
    return bIt - buffer.begin();
}

template<typename Stat_, typename Index_, typename Rank_>
void compute_min_rank_internal(Index_ use, const std::vector<std::pair<Stat_, Index_> >& buffer, Rank_* output) {
    Rank_ counter = 1;
    for (Index_ i = 0; i < use; ++i) {
        auto& current = output[buffer[i].second];
        if (counter < current) {
            current = counter;
        }
        ++counter;
    }
}

template<typename Stat_, typename Index_, typename Rank_>
void compute_min_rank_for_group(Index_ ngenes, size_t ngroups, size_t group, const Stat_* effects, Rank_* output, int threads) {
    std::vector<std::vector<Rank_> > stores(threads - 1);

    tatami::parallelize([&](size_t t, size_t start, size_t length) {
        Rank_* curoutput;
        if (t == 0) {
            curoutput = output;
            std::fill_n(curoutput, ngenes, ngenes + 1);
        } else {
            stores[t - 1].resize(ngenes, ngenes + 1);
            curoutput = stores[t - 1].data();
        }
        std::vector<std::pair<Stat_, Index_> > buffer(ngenes);

        for (size_t g = start, end = start + length; g < end; ++g) {
            if (g == group) {
                continue;
            }
            auto used = fill_and_sort_rank_buffer(effects + g, ngroups, buffer);
            compute_min_rank_internal(used, buffer, curoutput);
        }
    }, ngroups, threads);

    for (const auto& curstore : stores) {
        auto copy = output;
        for (auto x : curstore) {
            if (x < *copy) {
                *copy = x;
            }
            ++copy;
        }
    }
}

template<typename Stat_, typename Index_, typename Rank_>
void compute_min_rank_pairwise(Index_ ngenes, size_t ngroups, const Stat_* effects, std::vector<Rank_*>& output, int threads) {
    size_t shift = ngroups * ngroups;

    tatami::parallelize([&](size_t, size_t start, size_t length) {
        std::vector<std::pair<Stat_, Index_> > buffer(ngenes);
        for (size_t g = start, end = start + length; g < end; ++g) { 
            auto target = output[g];
            std::fill_n(target, ngenes, ngenes + 1); 
            auto base = effects + g * ngroups;

            for (size_t g2 = 0; g2 < ngroups; ++g2) {
                if (g == g2) {
                    continue;
                }
                auto used = fill_and_sort_rank_buffer(base + g2, shift, buffer);
                compute_min_rank_internal(used, buffer, target);
            }
        }
    }, ngroups, threads);
}

}
/**
 * @endcond
 */

}

#endif
