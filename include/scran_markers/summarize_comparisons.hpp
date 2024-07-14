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
 * @brief Pointers to arrays to hold the summary statistics.
 *
 * @tparam Stat_ Floating-point type for the statistics.
 * @tparam Rank_ Numeric type for the rank.
 */
template<typename Stat_ = double, typename Rank_ = int>
struct SummaryBuffers {
    /**
     * Pointer to an array of length equal to the number of genes,
     * to be filled with the minimum effect size for each gene.
     * If NULL, the minimum is not computed.
     */ 
    Stat_* min = NULL;

    /**
     * Pointer to an array of length equal to the number of genes,
     * to be filled with the mean effect size for each gene.
     * If NULL, the mean is not computed.
     */ 
    Stat_* mean = NULL;

    /**
     * Pointer to an array of length equal to the number of genes,
     * to be filled with the median effect size for each gene.
     * If NULL, the median is not computed.
     */ 
    Stat_* median = NULL;

    /**
     * Pointer to an array of length equal to the number of genes,
     * to be filled with the maximum effect size for each gene.
     * If NULL, the maximum is not computed.
     */ 
    Stat_* max = NULL;

    /**
     * Pointer to an array of length equal to the number of genes,
     * to be filled with the minimum rank of the effect sizes for each gene.
     * If NULL, the minimum rank is not computed.
     */ 
    Rank_* min_rank = NULL;
};

/**
 * @brief Container for the summary statistics.
 *
 * @tparam Stat_ Floating-point type for the statistics.
 * @tparam Rank_ Numeric type for the rank.
 */
template<typename Stat_ = double, typename Rank_ = int>
struct SummaryResults {
    /**
     * Vector of length equal to the number of genes,
     * to be filled with the minimum effect size for each gene.
     */ 
    std::vector<Stat_> min;

    /**
     * Vector of length equal to the number of genes,
     * to be filled with the mean effect size for each gene.
     */ 
    std::vector<Stat_> mean;

    /**
     * Vector of length equal to the number of genes,
     * to be filled with the median effect size for each gene.
     */
    std::vector<Stat_> median;

    /**
     * Vector of length equal to the number of genes,
     * to be filled with the maximum effect size for each gene.
     */
    std::vector<Stat_> max;

    /**
     * Vector of length equal to the number of genes,
     * to be filled with the minimum rank of the effect sizes for each gene.
     */ 
    std::vector<Rank_> min_rank;
};

/**
 * @cond
 */
namespace internal {

template<typename Stat_, typename Rank_>
void summarize_comparisons(size_t ngroups, const Stat_* effects, size_t group, size_t gene, const SummaryBuffers<Stat_, Rank_>& output, std::vector<Stat_>& buffer) {
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
    if (ncomps <= 1) {
        Stat_ val = (ncomps == 0 ? std::numeric_limits<Stat_>::quiet_NaN() : *ebegin);
        if (output.min) {
            output.min[gene] = val;
        }
        if (output.mean) {
            output.mean[gene] = val;
        }
        if (output.max) {
            output.max[gene] = val;
        }
        if (output.median) {
            output.median[gene] = val;
        }

    } else {
        if (output.min) {
            output.min[gene] = *std::min_element(ebegin, elast);
        }
        if (output.mean) {
            output.mean[gene] = std::accumulate(ebegin, elast, 0.0) / ncomps; 
        }
        if (output.max) {
            output.max[gene] = *std::max_element(ebegin, elast);
        }
        if (output.median) { // this mutates the buffer, so we put this last to avoid surprises.
            output.median[gene] = tatami_stats::medians::direct(ebegin, ncomps, /* skip_nan = */ false); 
        }
    }
}

template<typename Stat_, typename Rank_>
void summarize_comparisons(size_t ngenes, size_t ngroups, const Stat_* effects, const std::vector<SummaryBuffers<Stat_, Rank_> >& output, int threads) {
    size_t shift = ngroups * ngroups;

    tatami::parallelize([&](size_t, size_t start, size_t length) -> void {
        std::vector<Stat_> buffer(ngroups);
        auto effect_ptr = effects + start * shift; // everything is already a size_t, so it's fine.

        for (size_t gene = start, end = start + length; gene < end; ++gene, effect_ptr += shift) {
            auto current_effects = effect_ptr;
            for (size_t l = 0; l < ngroups; ++l, current_effects += ngroups) {
                summarize_comparisons(ngroups, current_effects, l, gene, output[l], buffer);
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
void compute_min_rank_pairwise(Index_ ngenes, size_t ngroups, const Stat_* effects, const std::vector<SummaryBuffers<Stat_, Rank_> >& output, int threads) {
    size_t shift = ngroups * ngroups;

    tatami::parallelize([&](size_t, size_t start, size_t length) {
        std::vector<std::pair<Stat_, Index_> > buffer(ngenes);
        for (size_t g = start, end = start + length; g < end; ++g) { 
            auto target = output[g].min_rank;
            if (target == NULL) {
                continue;
            }

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

template<typename Stat_, typename Rank_>
SummaryBuffers<Stat_, Rank_> fill_summary_results(
    size_t ngenes,
    SummaryResults<Stat_, Rank_>& out, 
    bool compute_min,
    bool compute_mean,
    bool compute_median,
    bool compute_max,
    bool compute_min_rank) 
{
    SummaryBuffers<Stat_, Rank_> ptr;

    if (compute_min) {
        out.min.resize(ngenes);
        ptr.min = out.min.data();
    }
    if (compute_mean) {
        out.mean.resize(ngenes);
        ptr.mean = out.mean.data();
    }
    if (compute_median) {
        out.median.resize(ngenes);
        ptr.median = out.median.data();
    }
    if (compute_max) {
        out.max.resize(ngenes);
        ptr.max = out.max.data();
    }
    if (compute_min_rank) {
        out.min_rank.resize(ngenes);
        ptr.min_rank = out.min_rank.data();
    }

    return ptr;
}

template<typename Stat_, typename Rank_>
std::vector<SummaryBuffers<Stat_, Rank_> > fill_summary_results(
    size_t ngenes,
    size_t ngroups,
    std::vector<SummaryResults<Stat_, Rank_> >& outputs, 
    bool compute_min,
    bool compute_mean,
    bool compute_median,
    bool compute_max,
    bool compute_min_rank) 
{
    outputs.resize(ngroups);
    std::vector<SummaryBuffers<Stat_, Rank_> > ptrs;
    ptrs.reserve(ngroups);
    for (size_t g = 0; g < ngroups; ++g) {
        ptrs.emplace_back(fill_summary_results(ngenes, outputs[g], compute_min, compute_mean, compute_median, compute_max, compute_min_rank));
    }
    return ptrs;
}

}
/**
 * @endcond
 */

}

#endif
