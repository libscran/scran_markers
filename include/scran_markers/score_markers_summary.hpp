#ifndef SCRAN_SCORE_MARKERS_HPP
#define SCRAN_SCORE_MARKERS_HPP

#include "scran_blocks/scran_blocks.hpp"
#include "tatami/tatami.hpp"
#include "sanisizer/sanisizer.hpp"

#include <array>
#include <map>
#include <vector>

#include "scan_matrix.hpp"
#include "cohens_d.hpp"
#include "simple_diff.hpp"
#include "summarize_comparisons.hpp"
#include "average_group_stats.hpp"
#include "create_combinations.hpp"
#include "utils.hpp"

/**
 * @file score_markers_summary.hpp
 * @brief Score potential markers by summaries of effect sizes between pairs of groups of cells.
 */

namespace scran_markers {

/**
 * @brief Options for `score_markers_summary()` and friends.
 */
struct ScoreMarkersSummaryOptions {
    /**
     * Threshold on the differences in expression values, used to adjust the Cohen's d and AUC calculations.
     * This should be non-negative.
     * Higher thresholds will favor genes with large differences at the expense of those with low variance. 
     */
    double threshold = 0;

    /**
     * Number of threads to use. 
     * The parallelization scheme is determined by `tatami::parallelize()`.
     */
    int num_threads = 1;

    /** 
     * Size of the cache, in terms of the number of pairwise comparisons.
     * Larger values speed up the comparisons at the cost of higher memory consumption.
     */
    int cache_size = 100;

    /**
     * Whether to compute Cohen's d. 
     * This only affects the `score_markers_summary()` overload that returns a `ScoreMarkersSummaryResults`.
     */
    bool compute_cohens_d = true;

    /**
     * Whether to compute the AUC.
     * This only affects the `score_markers_summary()` overload that returns a `ScoreMarkersSummaryResults`.
     */
    bool compute_auc = true;

    /**
     * Whether to compute the difference in means.
     * This only affects the `score_markers_summary()` overload that returns a `ScoreMarkersSummaryResults`.
     */
    bool compute_delta_mean = true;

    /**
     * Whether to compute the difference in the detected proportion.
     * This only affects the `score_markers_summary()` overload that returns a `ScoreMarkersSummaryResults`.
     */
    bool compute_delta_detected = true;

    /**
     * Whether to report the minimum of the effect sizes for each group.
     * Only affects the `score_markers_summary()` overload that returns a `ScoreMarkersSummaryResults`.
     */
    bool compute_min = true;

    /**
     * Whether to report the mean of the effect sizes for each group.
     * Only affects the `score_markers_summary()` overload that returns a `ScoreMarkersSummaryResults`.
     */
    bool compute_mean = true;

    /**
     * Whether to report the median of the effect sizes for each group.
     * Only affects the `score_markers_summary()` overload that returns a `ScoreMarkersSummaryResults`.
     */
    bool compute_median = true;

    /**
     * Whether to report the maximum of the effect sizes for each group.
     * Only affects the `score_markers_summary()` overload that returns a `ScoreMarkersSummaryResults`.
     */
    bool compute_max = true;

    /**
     * Whether to report the minimum rank of the effect sizes for each group.
     * Only affects the `score_markers_summary()` overload that returns a `ScoreMarkersSummaryResults`.
     */
    bool compute_min_rank = true;

    /**
     * Policy to use for weighting blocks when computing average statistics/effect sizes across blocks.
     *
     * The default of `scran_blocks::WeightPolicy::VARIABLE` is to define equal weights for blocks once they reach a certain size
     * (see `ScoreMarkersPairwiseOptions::variable_block_weight_parameters`).
     * For smaller blocks, the weight is linearly proportional to its size to avoid outsized contributions from very small blocks.
     *
     * Other options include `scran_blocks::WeightPolicy::EQUAL`, where all blocks are equally weighted regardless of size;
     * and `scran_blocks::WeightPolicy::NONE`, where the contribution of each block is proportional to its size.
     */
    scran_blocks::WeightPolicy block_weight_policy = scran_blocks::WeightPolicy::VARIABLE;

    /**
     * Parameters for the variable block weights, including the threshold at which blocks are considered to be large enough to have equal weight.
     * Only used when `ScoreMarkersSummaryOptions::block_weight_policy = scran_blocks::WeightPolicy::VARIABLE`.
     */
    scran_blocks::VariableWeightParameters variable_block_weight_parameters;
};

/**
 * @brief Buffers for `score_markers_summary()` and friends.
 * @tparam Stat_ Floating-point type of the output statistics.
 * @tparam Rank_ Numeric type of the rank.
 */
template<typename Stat_, typename Rank_>
struct ScoreMarkersSummaryBuffers {
    /**
     * Vector of length equal to the number of groups.
     * Each pointer corresponds to a group and points to an array of length equal to the number of genes,
     * to be filled with the mean expression of each gene in that group. 
     */
    std::vector<Stat_*> mean;

    /**
     * Vector of length equal to the number of groups.
     * Each pointer corresponds to a group and points to an array of length equal to the number of genes,
     * to be filled with the proportion of cells with detected expression in that group. 
     */
    std::vector<Stat_*> detected;

    /**
     * Vector of length equal to the number of groups.
     * Each entry contains the buffers in which to store the corresponding group's summary statistics for Cohen's d.
     *
     * Any of the pointers in any of the `SummaryBuffers` may be NULL, in which case the corresponding summary statistic is not computed.
     * This vector may also be empty, in which case no summary statistics are computed for this effect size.
     */
    std::vector<SummaryBuffers<Stat_, Rank_> > cohens_d;

    /**
     * Vector of length equal to the number of groups.
     * Each entry contains the buffers in which to store the corresponding group's summary statistics for the AUC.
     *
     * Any of the pointers in any of the `SummaryBuffers` may be NULL, in which case the corresponding summary statistic is not computed.
     * This vector may also be empty, in which case no summary statistics are computed for this effect size.
     */
    std::vector<SummaryBuffers<Stat_, Rank_> > auc;

    /**
     * Vector of length equal to the number of groups.
     * Each entry contains the buffers in which to store the corresponding group's summary statistics for the difference in means.
     *
     * Any of the pointers in any of the `SummaryBuffers` may be NULL, in which case the corresponding summary statistic is not computed.
     * This vector may also be empty, in which case no summary statistics are computed for this effect size.
     */
    std::vector<SummaryBuffers<Stat_, Rank_> > delta_mean;

    /**
     * Vector of length equal to the number of groups.
     * Each entry contains the buffers in which to store the corresponding group's summary statistics for the difference in the detected proportions.
     *
     * Any of the pointers in any of the `SummaryBuffers` may be NULL, in which case the corresponding summary statistic is not computed.
     * This vector may also be empty, in which case no summary statistics are computed for this effect size.
     */
    std::vector<SummaryBuffers<Stat_, Rank_> > delta_detected;
};

/**
 * @brief Results for `score_markers_summary()` and friends.
 * @tparam Stat_ Floating-point type of the output statistics.
 * @tparam Rank_ Numeric type of the rank.
 */
template<typename Stat_, typename Rank_>
struct ScoreMarkersSummaryResults {
    /**
     * Vector of length equal to the number of groups.
     * Each inner vector corresponds to a group and contains the mean expression of each gene in that group. 
     */
    std::vector<std::vector<Stat_> > mean;

    /**
     * Vector of length equal to the number of groups.
     * Each inner vector corresponds to a group and contains the mean expression of each gene in that group. 
     */
    std::vector<std::vector<Stat_> > detected;

    /**
     * Vector of length equal to the number of groups, containing the summaries of the Cohen's d for each group.
     * This may be an empty vector if `ScoreMarkersSummaryOptions::compute_cohens_d = false`.
     *
     * Individual vectors inside the `SummaryResults` may also be empty if specified by the relevant option,
     * e.g., `ScoreMarkersSummaryOptions::compute_min = false` will cause `SummaryResults::min` to be empty.
     */
    std::vector<SummaryResults<Stat_, Rank_> > cohens_d;

    /**
     * Vector of length equal to the number of groups, containing the summaries of the AUC for each group.
     * This may be an empty vector if `ScoreMarkersSummaryOptions::compute_auc = false`.
     *
     * Individual vectors inside the `SummaryResults` may also be empty if specified by the relevant option,
     * e.g., `ScoreMarkersSummaryOptions::compute_min = false` will cause `SummaryResults::min` to be empty.
     */
    std::vector<SummaryResults<Stat_, Rank_> > auc;

    /**
     * Vector of length equal to the number of groups, containing the summaries of the differences in means for each group.
     * This may be an empty vector if `ScoreMarkersSummaryOptions::compute_delta_mean = false`.
     *
     * Individual vectors inside the `SummaryResults` may also be empty if specified by the relevant option,
     * e.g., `ScoreMarkersSummaryOptions::compute_min = false` will cause `SummaryResults::min` to be empty.
     */
    std::vector<SummaryResults<Stat_, Rank_> > delta_mean;

    /**
     * Vector of length equal to the number of groups, containing the summaries of the differences in detected proportions for each group.
     * This may be an empty vector if `ScoreMarkersSummaryOptions::compute_delta_detected = false`.
     *
     * Individual vectors inside the `SummaryResults` may also be empty if specified by the relevant option,
     * e.g., `ScoreMarkersSummaryOptions::compute_min = false` will cause `SummaryResults::min` to be empty.
     */
    std::vector<SummaryResults<Stat_, Rank_> > delta_detected;
};

/**
 * @cond
 */
namespace internal {

enum class CacheAction : unsigned char { SKIP, COMPUTE, CACHE };

// Safely cap the cache size at the maximum possible number of comparisons (i.e., ngroups * (ngroups - 1) / 2).
inline std::size_t cap_cache_size(const std::size_t cache_size, const std::size_t ngroups) {
    if (ngroups < 2) {
        return 0;
    } else if (ngroups % 2 == 0) {
        const auto denom = ngroups / 2;
        const auto ratio = cache_size / denom;
        if (ratio <= ngroups - 1) {
            return cache_size;
        } else {
            return denom * (ngroups - 1);
        }
    } else {
        const auto denom = (ngroups - 1) / 2;
        const auto ratio = cache_size / denom;
        if (ratio <= ngroups) {
            return cache_size;
        } else {
            return denom * ngroups;
        }
    }
}

/* 
 * We compute effect sizes in a pairwise fashion with some nested loops, i.e.,
 * iterating from g1 in [1, G) and then for g2 in [0, g1). When we compute the
 * effect size for g1 vs g2, we sometimes get a free effect size for g2 vs g1,
 * which we call the "reverse effect". However, we can't use the reverse effect
 * size until we get around to summarizing effects for g2, hence the caching.
 * 
 * This cache tries to store as many of the reverse effects as possible before
 * it starts evicting. Evictions are based on the principle that it is better
 * to store effects that will be re-used quickly, thus freeing up the cache for
 * future stores. The 'speed' of reusability of each cache entry depends on the
 * first group in the comparison corresponding to each cached effect size; the
 * smaller the first group, the sooner it will be reached when iterating across
 * groups in the ScoreMarkers function.
 *
 * So, the policy is to evict cache entries when the identity of the first
 * group in the cached entry is larger than the identity of the first group for
 * the incoming entry. Given that, if the cache is full, we have to throw away
 * one of these effects anyway, I'd prefer to hold onto the one we're using
 * soon, because at least it'll be freed up rapidly.
 */
template<typename Index_, typename Stat_>
class EffectsCacher {
public:
    EffectsCacher(const Index_ ngenes, const std::size_t ngroups, const std::size_t cache_size) :
        my_ngenes(ngenes),
        my_ngroups(ngroups),
        my_cache_size(cap_cache_size(cache_size, ngroups)),
        my_actions(sanisizer::cast<decltype(I(my_actions.size()))>(ngroups)),
        my_common_cache(sanisizer::product<decltype(I(my_common_cache.size()))>(ngenes, my_cache_size)),
        my_staging_cache(sanisizer::cast<decltype(I(my_staging_cache.size()))>(ngroups))
    {
        my_unused_pool.reserve(my_cache_size);
        for (decltype(I(my_cache_size)) c = 0; c < my_cache_size; ++c) {
            my_unused_pool.push_back(my_common_cache.data() + sanisizer::product_unsafe<std::size_t>(c, ngenes));
        }
    }

private:
    Index_ my_ngenes;
    std::size_t my_ngroups;
    std::size_t my_cache_size;

    std::vector<CacheAction> my_actions;

    // 'common_cache' contains allocation so that we don't have to do
    // a lot of fiddling with move constructors, swaps, etc.
    std::vector<Stat_> my_common_cache;

    // 'staging_cache' contains the set of cached effects in the other
    // direction, i.e., all other groups compared to the current group. This is
    // only used to avoid repeated look-ups in 'cached' while filling the
    // effect size vectors; they will ultimately be transferred to cached after
    // the processing for the current group is complete.
    std::vector<Stat_*> my_staging_cache;

    // 'unused_pool' contains the currently-unused set of pointers to free 
    // subarrays in 'my_common_cache'
    std::vector<Stat_*> my_unused_pool;

    // 'cached' contains the cached effect size vectors from previous groups. Note
    // that the use of a map is deliberate as we need the sorting.
    std::map<std::pair<std::size_t, std::size_t>, Stat_*> my_cached;

public:
    void clear() {
        my_cached.clear();
    }

public:
    void fill_effects_from_cache(const std::size_t group, std::vector<double>& full_effects) {
        // During calculation of effects, the current group (i.e., 'group') is
        // the first group in the comparison and 'other' is the second group.
        // However, remember that we want to cache the reverse effects, so in
        // the cached entry, 'group' is second and 'other' is first.
        for (decltype(I(my_ngroups)) other = 0; other < my_ngroups; ++other) {
            if (other == group) {
                my_actions[other] = CacheAction::SKIP;
                continue;
            }

            if (my_cache_size == 0) {
                my_actions[other] = CacheAction::COMPUTE;
                continue;
            }

            // If 'other' is later than 'group', it's a candidate to be cached,
            // as it will be used when this method is called with 'group = other'.
            // Note that the ACTUAL caching decision is made in the refinement step below.
            if (other > group) {
                my_actions[other] = CacheAction::CACHE;
                continue;
            }

            // Need to recompute cache entries that were previously evicted. We
            // do so if the cache is empty or the first group of the first cached
            // entry has a higher index than the current group (and thus the 
            // desired comparison's effects cannot possibly exist in the cache).
            if (my_cached.empty()) { 
                my_actions[other] = CacheAction::COMPUTE;
                continue;
            }

            const auto& front = my_cached.begin()->first;
            if (front.first > group || front.second > other) { 
                // Technically, the second clause should be (front.first == group && front.second > other).
                // However, less-thans should be impossible as they should have been used up during processing
                // of previous 'group' values. Thus, equality is already implied if the first condition fails.
                my_actions[other] = CacheAction::COMPUTE;
                continue;
            }

            // If we got past the previous clause, this implies that the first cache entry
            // contains the effect sizes for the desired comparison (i.e., 'group - other').
            // We thus transfer the cached vector to the full_set.
            my_actions[other] = CacheAction::SKIP;
            const auto curcache = my_cached.begin()->second;
            for (decltype(I(my_ngenes)) i = 0; i < my_ngenes; ++i) {
                full_effects[sanisizer::nd_offset<std::size_t>(other, my_ngroups, i)] = curcache[i];
            }

            my_unused_pool.push_back(curcache);
            my_cached.erase(my_cached.begin());
        }

        // Refining our choice of cacheable entries by doing a dummy run and
        // seeing whether eviction actually happens. If it doesn't, we won't
        // bother caching, because that would be a waste of memory accesses.
        for (decltype(I(my_ngroups)) other = 0; other < my_ngroups; ++other) {
            if (my_actions[other] != CacheAction::CACHE) {
                continue;
            }

            const std::pair<std::size_t, std::size_t> key(other, group);
            if (my_cached.size() < my_cache_size) {
                const auto ptr = my_unused_pool.back();
                my_cached[key] = ptr;
                my_staging_cache[other] = ptr;
                my_unused_pool.pop_back();
                continue;
            }

            // Looking at the last cache entry. If the first group of this
            // entry is larger than the first group of the incoming entry, we
            // evict it, as the incoming entry has faster reusability.
            auto it = my_cached.end();
            --it;
            if ((it->first).first > other) {
                const auto ptr = it->second;
                my_cached[key] = ptr;
                my_staging_cache[other] = ptr;
                my_cached.erase(it);
            } else {
                // Otherwise, if we're not going to do any evictions, we
                // indicate that we shouldn't even bother computing the
                // reverse, because we're won't cache the incoming entry.
                my_actions[other] = CacheAction::COMPUTE;
            }
        }
    }

public:
    CacheAction get_action(std::size_t other) const {
        return my_actions[other];
    }

    Stat_* get_cache_location(std::size_t other) const {
        return my_staging_cache[other];
    }
};

template<typename Index_, typename Stat_, typename Rank_>
void process_simple_summary_effects(
    const Index_ ngenes,
    const std::size_t ngroups,
    const std::size_t nblocks,
    const std::size_t ncombos,
    const std::vector<Stat_>& combo_means,
    const std::vector<Stat_>& combo_vars,
    const std::vector<Stat_>& combo_detected,
    const ScoreMarkersSummaryBuffers<Stat_, Rank_>& output,
    const std::vector<Stat_>& combo_weights,
    const double threshold,
    const std::size_t cache_size,
    const int num_threads)
{
    // First, computing the pooled averages to get that out of the way.
    {
        std::vector<Stat_> total_weights_per_group;
        auto total_weights_ptr = combo_weights.data();
        if (nblocks > 1) {
            total_weights_per_group = compute_total_weight_per_group(ngroups, nblocks, combo_weights.data());
            total_weights_ptr = total_weights_per_group.data();
        }

        tatami::parallelize([&](const int, const Index_ start, const Index_ length) -> void {
            for (Index_ gene = start, end = start + length; gene < end; ++gene) {
                const auto in_offset = sanisizer::product_unsafe<std::size_t>(gene, ncombos);
                const auto tmp_means = combo_means.data() + in_offset;
                const auto tmp_detected = combo_detected.data() + in_offset;
                average_group_stats(gene, ngroups, nblocks, tmp_means, tmp_detected, combo_weights.data(), total_weights_ptr, output.mean, output.detected);
            }
        }, ngenes, num_threads);
    }

    PrecomputedPairwiseWeights<Stat_> preweights(ngroups, nblocks, combo_weights.data());
    EffectsCacher<Index_, Stat_> cache(ngenes, ngroups, cache_size);
    std::vector<Stat_> full_effects(sanisizer::product<typename std::vector<Stat_>::size_type>(ngroups, ngenes));
    auto effect_buffers = sanisizer::create<std::vector<std::vector<Stat_> > >(num_threads);
    for (auto& ef : effect_buffers) {
        sanisizer::resize(ef, ngroups);
    }

    if (output.cohens_d.size()) {
        cache.clear();
        for (decltype(I(ngroups)) group = 0; group < ngroups; ++group) {
            cache.fill_effects_from_cache(group, full_effects);

            tatami::parallelize([&](const int t, const Index_ start, const Index_ length) -> void {
                auto& effect_buffer = effect_buffers[t];

                for (Index_ gene = start, end = start + length; gene < end; ++gene) {
                    const auto in_offset = sanisizer::product_unsafe<std::size_t>(gene, ncombos);
                    const auto my_means = combo_means.data() + in_offset;
                    const auto my_variances = combo_vars.data() + in_offset;
                    const auto store_ptr = full_effects.data() + sanisizer::product_unsafe<std::size_t>(gene, ngroups);

                    for (decltype(I(ngroups)) other = 0; other < ngroups; ++other) {
                        auto cache_action = cache.get_action(other);
                        if (cache_action == internal::CacheAction::COMPUTE) {
                            store_ptr[other] = compute_pairwise_cohens_d_one_sided(group, other, my_means, my_variances, ngroups, nblocks, preweights, threshold);
                        } else if (cache_action == internal::CacheAction::CACHE) {
                            auto tmp = compute_pairwise_cohens_d_two_sided(group, other, my_means, my_variances, ngroups, nblocks, preweights, threshold);
                            store_ptr[other] = tmp.first;
                            cache.get_cache_location(other)[gene] = tmp.second;
                        }
                    }
                    summarize_comparisons(ngroups, store_ptr, group, gene, output.cohens_d[group], effect_buffer);
                }
            }, ngenes, num_threads);

            const auto mr = output.cohens_d[group].min_rank;
            if (mr) {
                compute_min_rank_for_group(ngenes, ngroups, group, full_effects.data(), mr, num_threads);
            }
        }
    }

    if (output.delta_mean.size()) {
        cache.clear();
        for (decltype(I(ngroups)) group = 0; group < ngroups; ++group) {
            cache.fill_effects_from_cache(group, full_effects);

            tatami::parallelize([&](const int t, const Index_ start, const Index_ length) -> void {
                auto& effect_buffer = effect_buffers[t];

                for (Index_ gene = start, end = start + length; gene < end; ++gene) {
                    const auto my_means = combo_means.data() + sanisizer::product_unsafe<std::size_t>(gene, ncombos);
                    const auto store_ptr = full_effects.data() + sanisizer::product_unsafe<std::size_t>(gene, ngroups);

                    for (decltype(I(ngroups)) other = 0; other < ngroups; ++other) {
                        const auto cache_action = cache.get_action(other);
                        if (cache_action != internal::CacheAction::SKIP) {
                            const auto val = compute_pairwise_simple_diff(group, other, my_means, ngroups, nblocks, preweights);
                            store_ptr[other] = val;
                            if (cache_action == CacheAction::CACHE) {
                                cache.get_cache_location(other)[gene] = -val;
                            } 
                        }
                    }
                    summarize_comparisons(ngroups, store_ptr, group, gene, output.delta_mean[group], effect_buffer);
                }
            }, ngenes, num_threads);

            const auto mr = output.delta_mean[group].min_rank;
            if (mr) {
                compute_min_rank_for_group(ngenes, ngroups, group, full_effects.data(), mr, num_threads);
            }
        }
    }

    if (output.delta_detected.size()) {
        cache.clear();
        for (decltype(I(ngroups)) group = 0; group < ngroups; ++group) {
            cache.fill_effects_from_cache(group, full_effects);

            tatami::parallelize([&](const int t, const Index_ start, const Index_ length) -> void {
                auto& effect_buffer = effect_buffers[t];

                for (Index_ gene = start, end = start + length; gene < end; ++gene) {
                    const auto my_detected = combo_detected.data() + sanisizer::product_unsafe<std::size_t>(gene, ncombos);
                    const auto store_ptr = full_effects.data() + sanisizer::product_unsafe<std::size_t>(gene, ngroups);

                    for (decltype(I(ngroups)) other = 0; other < ngroups; ++other) {
                        const auto cache_action = cache.get_action(other);
                        if (cache_action != CacheAction::SKIP) {
                            const auto val = compute_pairwise_simple_diff(group, other, my_detected, ngroups, nblocks, preweights);
                            store_ptr[other] = val;
                            if (cache_action == CacheAction::CACHE) {
                                cache.get_cache_location(other)[gene] = -val;
                            } 
                        }
                    }
                    summarize_comparisons(ngroups, store_ptr, group, gene, output.delta_detected[group], effect_buffer);
                }
            }, ngenes, num_threads);

            const auto mr = output.delta_detected[group].min_rank;
            if (mr) {
                compute_min_rank_for_group(ngenes, ngroups, group, full_effects.data(), mr, num_threads);
            }
        }
    }
}

template<typename Index_, typename Stat_, typename Rank_>
ScoreMarkersSummaryBuffers<Stat_, Rank_> fill_summary_results(
    const Index_ ngenes,
    const std::size_t ngroups,
    ScoreMarkersSummaryResults<Stat_, Rank_>& store,
    const ScoreMarkersSummaryOptions& options)
{
    ScoreMarkersSummaryBuffers<Stat_, Rank_> output;

    internal::fill_average_results(ngenes, ngroups, store.mean, store.detected, output.mean, output.detected);

    if (options.compute_cohens_d) {
        output.cohens_d = internal::fill_summary_results(
            ngenes,
            ngroups,
            store.cohens_d,
            options.compute_min,
            options.compute_mean,
            options.compute_median,
            options.compute_max,
            options.compute_min_rank
        );
    }

    if (options.compute_auc) {
        output.auc = internal::fill_summary_results(
            ngenes,
            ngroups,
            store.auc,
            options.compute_min,
            options.compute_mean,
            options.compute_median,
            options.compute_max,
            options.compute_min_rank
        );
    }

    if (options.compute_delta_mean) {
        output.delta_mean = internal::fill_summary_results(
            ngenes,
            ngroups,
            store.delta_mean,
            options.compute_min,
            options.compute_mean,
            options.compute_median,
            options.compute_max,
            options.compute_min_rank
        );
    }

    if (options.compute_delta_detected) {
        output.delta_detected = internal::fill_summary_results(
            ngenes,
            ngroups,
            store.delta_detected,
            options.compute_min,
            options.compute_mean,
            options.compute_median,
            options.compute_max,
            options.compute_min_rank
        );
    }

    return output;
}

template<typename Stat_, typename Rank_>
bool auc_needs_minrank(const std::vector<SummaryBuffers<Stat_, Rank_> >& auc_summaries) {
    for (const auto& outs : auc_summaries) {
        if (outs.min_rank) {
            return true;
        }
    }
    return false;
}

template<
    bool single_block_,
    typename Value_,
    typename Index_,
    typename Group_,
    typename Block_,
    typename Stat_,
    typename Rank_
>
void score_markers_summary(
    const tatami::Matrix<Value_, Index_>& matrix, 
    const std::size_t ngroups,
    const Group_* const group, 
    const std::size_t nblocks,
    const Block_* const block,
    const std::size_t ncombos,
    const std::size_t* const combo,
    const std::vector<Index_>& combo_sizes,
    const ScoreMarkersSummaryOptions& options,
    const ScoreMarkersSummaryBuffers<Stat_, Rank_>& output
) {
    const auto ngenes = matrix.nrow();
    const auto payload_size = sanisizer::product<typename std::vector<Stat_>::size_type>(ngenes, ncombos);
    std::vector<Stat_> combo_means(payload_size), combo_vars(payload_size), combo_detected(payload_size);

    // For a single block, this usually doesn't really matter, but we do it for consistency with the multi-block case,
    // and to account for variable weighting where non-zero block sizes get zero weight.
    const auto combo_weights = scran_blocks::compute_weights<Stat_>(
        combo_sizes,
        options.block_weight_policy,
        options.variable_block_weight_parameters
    );

    const bool do_auc = !output.auc.empty();

    if (do_auc || matrix.prefer_rows()) {
        if (do_auc && !auc_needs_minrank(output.auc)) {
            // If we don't need the min-rank, we can compute summaries for the AUCs directly.
            // This means that we don't have to store the full 3D array of AUCs across all genes.
            struct AucResultWorkspace {
                AucResultWorkspace() = default;
                AucResultWorkspace(const std::size_t ngroups) :
                    pairwise_buffer(sanisizer::product<typename std::vector<Stat_>::size_type>(ngroups, ngroups)),
                    summary_buffer(sanisizer::cast<typename std::vector<Stat_>::size_type>(ngroups))
                {};
                std::vector<Stat_> pairwise_buffer;
                std::vector<Stat_> summary_buffer;
            };

            scan_matrix_by_row_custom_auc<single_block_>(
                matrix, 
                ngroups,
                group,
                nblocks,
                block,
                ncombos,
                combo,
                combo_means,
                combo_vars,
                combo_detected,
                /* do_auc = */ true,
                /* auc_result_initialize = */ [&]() -> AucResultWorkspace {
                    return AucResultWorkspace(ngroups);
                },
                /* auc_result_process = */ [&](
                    const Index_ gene,
                    AucScanWorkspace<Value_, Group_, Index_, Stat_>& auc_work,
                    AucResultWorkspace& res_work
                ) -> void {
                    process_auc_for_rows(auc_work, ngroups, nblocks, options.threshold, res_work.pairwise_buffer.data());
                    for (decltype(I(ngroups)) l = 0; l < ngroups; ++l) {
                        const auto current_effects = res_work.pairwise_buffer.data() + sanisizer::product_unsafe<std::size_t>(l, ngroups);
                        summarize_comparisons(ngroups, current_effects, l, gene, output.auc[l], res_work.summary_buffer);
                    }
                },
                combo_sizes,
                combo_weights,
                options.num_threads
            );

        } else {
            std::vector<Stat_> tmp_auc;
            Stat_* auc_ptr = NULL;
            if (do_auc) {
                tmp_auc.resize(sanisizer::product<decltype(I(tmp_auc.size()))>(ngroups, ngroups, ngenes));
                auc_ptr = tmp_auc.data();
            } 

            scan_matrix_by_row_full_auc<single_block_>(
                matrix, 
                ngroups,
                group,
                nblocks,
                block,
                ncombos,
                combo,
                combo_means,
                combo_vars,
                combo_detected,
                auc_ptr,
                combo_sizes,
                combo_weights, 
                options.threshold,
                options.num_threads
            );

            if (do_auc) {
                summarize_comparisons(ngenes, ngroups, auc_ptr, output.auc, options.num_threads);
                compute_min_rank_pairwise(ngenes, ngroups, auc_ptr, output.auc, options.num_threads);
            }
        }

    } else {
        internal::scan_matrix_by_column(
            matrix,
            [&]{
                if constexpr(single_block_) {
                    return ngroups;
                } else {
                    return ncombos;
                }
            }(),
            [&]{
                if constexpr(single_block_) {
                    return group;
                } else {
                    return combo;
                }
            }(),
            combo_means,
            combo_vars,
            combo_detected,
            combo_sizes,
            options.num_threads
        );
    }
 
    internal::process_simple_summary_effects(
        matrix.nrow(),
        ngroups,
        nblocks,
        ncombos,
        combo_means,
        combo_vars,
        combo_detected,
        output,
        combo_weights,
        options.threshold,
        options.cache_size,
        options.num_threads
    );
}

}
/**
 * @endcond
 */

/**
 * Score each gene as a candidate marker for each group of cells, based on summaries of effect sizes from pairwise comparisons between groups.
 *
 * Markers are identified by differential expression analyses between pairs of groups of cells (e.g., clusters, cell types).
 * Given \f$N\f$ groups, each group is involved in \f$N - 1\f$ pairwise comparisons and thus has \f$N - 1\f$ effect sizes for each gene.
 * We summarize each group's effect sizes into a small set of desriptive statistics like the minimum, median or mean.
 * Users can then sort genes by any of these summaries to obtain a ranking of potential markers for the group.
 *
 * The choice of effect size and summary statistic determines the characteristics of the marker ranking.
 * The effect sizes include Cohen's d, the area under the curve (AUC), the delta-mean and the delta-detected (see `score_markers_pairwise()`).
 * The summary statistics include the minimum, mean, median, maximum and min-rank of the effect sizes across each group's pairwise comparisons (see `summarize_effects()`).
 * For example, ranking by the delta-detected with the minimum summary will promote markers that are silent in every other group.
 *
 * This behavior of this function is equivalent to - but more efficient than - calling `score_markers_pairwise()` followed by `summarize_effects()` on each array of effect sizes.
 *
 * @tparam Value_ Matrix data type.
 * @tparam Index_ Matrix index type.
 * @tparam Group_ Integer type of the group assignments.
 * @tparam Stat_ Floating-point type to store the statistics.
 * @tparam Rank_ Numeric type to store the minimum rank.
 *
 * @param matrix A matrix of expression values, typically normalized and log-transformed.
 * Rows should contain genes while columns should contain cells.
 * @param[in] group Pointer to an array of length equal to the number of columns in `matrix`, containing the group assignments.
 * Group identifiers should be 0-based and should contain all integers in \f$[0, N)\f$ where \f$N\f$ is the number of unique groups.
 * @param options Further options.
 * @param[out] output Collection of buffers in which to store the computed statistics.
 * Each buffer is filled with the corresponding statistic for each group or pairwise comparison.
 * Any of `ScoreMarkersSummaryBuffers::cohens_d`, 
 * `ScoreMarkersSummaryBuffers::auc`, 
 * `ScoreMarkersSummaryBuffers::delta_mean` or
 * `ScoreMarkersSummaryBuffers::delta_detected`
 * may be empty, in which case the corresponding statistic is not computed or summarized.
 */
template<typename Value_, typename Index_, typename Group_, typename Stat_, typename Rank_>
void score_markers_summary(
    const tatami::Matrix<Value_, Index_>& matrix, 
    const Group_* const group, 
    const ScoreMarkersSummaryOptions& options,
    const ScoreMarkersSummaryBuffers<Stat_, Rank_>& output
) {
    const auto NC = matrix.ncol();
    const auto group_sizes = tatami_stats::tabulate_groups(group, NC); 
    const auto ngroups = sanisizer::cast<std::size_t>(group_sizes.size());

    internal::score_markers_summary<true>(
        matrix,
        ngroups,
        group,
        1,
        static_cast<int*>(NULL),
        ngroups,
        static_cast<std::size_t*>(NULL),
        group_sizes,
        options,
        output
    );
}

/**
 * Score potential marker genes by computing summary statistics across pairwise comparisons between groups, accounting for any blocking factor in the dataset.
 * Comparisons are only performed between the groups of cells in the same level of the blocking factor, as described in `score_markers_pairwise_blocked()`.
 * This strategy avoids most problems related to batch effects as we never directly compare across different blocking levels.
 * The block-specific effect sizes are combined into a single aggregate value per comparison, which are in turn summarized as described in `summarize_effects()`.
 * This behavior of this function is equivalent to - but more efficient than - calling `score_markers_pairwise_blocked()` followed by `summarize_effects()` on each array of effect sizes.
 *
 * @tparam Value_ Matrix data type.
 * @tparam Index_ Matrix index type.
 * @tparam Group_ Integer type of the group assignments.
 * @tparam Stat_ Floating-point type to store the statistics.
 * @tparam Rank_ Numeric type to store the minimum rank.
 *
 * @param matrix A matrix of expression values, typically normalized and log-transformed.
 * Rows should contain genes while columns should contain cells.
 * @param[in] group Pointer to an array of length equal to the number of columns in `matrix`, containing the group assignments.
 * Group identifiers should be 0-based and should contain all integers in \f$[0, N)\f$ where \f$N\f$ is the number of unique groups.
 * @param[in] block Pointer to an array of length equal to the number of columns in `matrix`, containing the blocking factor.
 * Block identifiers should be 0-based and should contain all integers in \f$[0, B)\f$ where \f$B\f$ is the number of unique blocking levels.
 * @param options Further options.
 * @param[out] output Collection of buffers in which to store the computed statistics.
 * Each buffer is filled with the corresponding statistic for each group or pairwise comparison.
 * Any of `ScoreMarkersSummaryBuffers::cohens_d`, 
 * `ScoreMarkersSummaryBuffers::auc`, 
 * `ScoreMarkersSummaryBuffers::delta_mean` or
 * `ScoreMarkersSummaryBuffers::delta_detected`
 * may be empty, in which case the corresponding statistic is not computed or summarized.
 */
template<typename Value_, typename Index_, typename Group_, typename Block_, typename Stat_, typename Rank_>
void score_markers_summary_blocked(
    const tatami::Matrix<Value_, Index_>& matrix, 
    const Group_* const group, 
    const Block_* const block,
    const ScoreMarkersSummaryOptions& options,
    const ScoreMarkersSummaryBuffers<Stat_, Rank_>& output) 
{
    const auto NC = matrix.ncol();
    const auto ngroups = output.mean.size();
    const auto nblocks = tatami_stats::total_groups(block, NC); 

    const auto combinations = internal::create_combinations(ngroups, group, block, NC);
    const auto combo_sizes = internal::tabulate_combinations<Index_>(ngroups, nblocks, combinations);
    const auto ncombos = combo_sizes.size();

    internal::score_markers_summary<false>(
        matrix,
        sanisizer::cast<std::size_t>(ngroups),
        group,
        sanisizer::cast<std::size_t>(nblocks),
        block,
        sanisizer::cast<std::size_t>(ncombos),
        combinations.data(),
        combo_sizes,
        options,
        output
    );
}


/**
 * Overload of `score_markers_pairwise()` that allocates memory for the output statistics.
 *
 * @tparam Stat_ Floating-point type to store the statistics.
 * @tparam Rank_ Numeric type to store the minimum rank.
 * @tparam Value_ Matrix data type.
 * @tparam Index_ Matrix index type.
 * @tparam Group_ Integer type of the group assignments.
 *
 * @param matrix A matrix of expression values, typically normalized and log-transformed.
 * Rows should contain genes while columns should contain cells.
 * @param[in] group Pointer to an array of length equal to the number of columns in `matrix`, containing the group assignments.
 * Group identifiers should be 0-based and should contain all integers in \f$[0, N)\f$ where \f$N\f$ is the number of unique groups.
 * @param options Further options.
 *
 * @return Object containing the summary statistics and the other per-group statistics.
 */
template<typename Stat_ = double, typename Rank_ = int, typename Value_, typename Index_, typename Group_>
ScoreMarkersSummaryResults<Stat_, Rank_> score_markers_summary(
    const tatami::Matrix<Value_, Index_>& matrix,
    const Group_* const group,
    const ScoreMarkersSummaryOptions& options)
{
    const auto ngroups = tatami_stats::total_groups(group, matrix.ncol());
    ScoreMarkersSummaryResults<Stat_, Rank_> output;
    const auto buffers = internal::fill_summary_results(matrix.nrow(), ngroups, output, options);
    score_markers_summary(matrix, group, options, buffers);
    return output;
}

/**
 * Overload of `score_markers_pairwise_blocked()` that allocates memory for the output statistics.
 *
 * @tparam Stat_ Floating-point type to store the statistics.
 * @tparam Rank_ Numeric type to store the minimum rank.
 * @tparam Value_ Matrix data type.
 * @tparam Index_ Matrix index type.
 * @tparam Group_ Integer type of the group assignments.
 * @tparam Block_ Integer type of the block assignments. 
 *
 * @param matrix A matrix of expression values, typically normalized and log-transformed.
 * Rows should contain genes while columns should contain cells.
 * @param[in] group Pointer to an array of length equal to the number of columns in `matrix`, containing the group assignments.
 * Group identifiers should be 0-based and should contain all integers in \f$[0, N)\f$ where \f$N\f$ is the number of unique groups.
 * @param[in] block Pointer to an array of length equal to the number of columns in `matrix`, containing the blocking factor.
 * Block identifiers should be 0-based and should contain all integers in \f$[0, B)\f$ where \f$B\f$ is the number of unique blocking levels.
 * @param options Further options.
 *
 * @return Object containing the pairwise effects, plus the mean expression and detected proportion in each group.
 */
template<typename Stat_ = double, typename Rank_ = int, typename Value_, typename Index_, typename Group_, typename Block_>
ScoreMarkersSummaryResults<Stat_, Rank_> score_markers_summary_blocked(
    const tatami::Matrix<Value_, Index_>& matrix,
    const Group_* const group,
    const Block_* const block,
    const ScoreMarkersSummaryOptions& options)
{
    const auto ngroups = tatami_stats::total_groups(group, matrix.ncol());
    ScoreMarkersSummaryResults<Stat_, Rank_> output;
    const auto buffers = internal::fill_summary_results(matrix.nrow(), ngroups, output, options);
    score_markers_summary_blocked(matrix, group, block, options, buffers);
    return output;
}

}

#endif
