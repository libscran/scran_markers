#ifndef SCRAN_SCORE_MARKERS_HPP
#define SCRAN_SCORE_MARKERS_HPP

#include "scan_matrix.hpp"
#include "cohens_d.hpp"
#include "simple_diff.hpp"
#include "summarize_comparisons.hpp"

#include "scran_blocks/scran_blocks.hpp"
#include "tatami/tatami.hpp"

#include <array>
#include <map>
#include <vector>

/**
 * @file ScoreMarkers.hpp
 *
 * @brief Compute marker scores for each gene in each group of cells.
 */

namespace scran_markers {

/**
 * @brief Options for `score_markers_summary()` and friends.
 */
struct ScoreMarkersSummaryOptions {
    /**
     * Threshold on the differences in expression values, used to adjust the Cohen's D and AUC calculations.
     * This should be non-negative.
     */
    double threshold = 0;

    /**
     * Number of threads to use. 
     */
    int num_threads = 1;

    /** 
     * Size of the cache, in terms of the number of pairwise comparisons.
     * Larger values speed up the comparisons at the cost of higher memory consumption.
     */
    int cache_size = 100;

    /**
     * Whether to compute Cohen's d. 
     * This only affects the `score_markers_summary()` overload that return `Results`.
     */
    bool compute_cohens_d = true;

    /**
     * Whether to compute the AUC.
     * This only affects the `score_markers_summary()` overload that return `Results`.
     */
    bool compute_auc = true;

    /**
     * Whether to compute the difference in means.
     * This only affects the `score_markers_summary()` overload that return `Results`.
     */
    bool compute_delta_mean = true;

    /**
     * Whether to compute the difference in the detected proportion.
     * This only affects the `score_markers_summary()` overload that return `Results`.
     */
    bool compute_delta_detected = true;

    /**
     * Policy to use for weighting blocks when computing average statistics/effect sizes across blocks.
     */
    scran_blocks::WeightPolicy block_weight_policy = scran_blocks::WeightPolicy::VARIABLE;

    /**
     * Parameters for the variable block weights. 
     * Only used when the block weight policy is set to `scran_blocks::WeightPolicy::VARIABLE`.
     */
    scran_blocks::VariableWeightParameters variable_block_weight_parameters;
};

/**
 * @brief Buffers for `score_markers_summary()` and friends.
 * @tparam Stat_ Floating-point type for the output statistics.
 * @tparam Rank_ Numeric type for the rank.
 */
template<typename Stat_>
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
     */
    std::vector<SummaryBuffers<Stat_, Rank_> > cohens_d;

    /**
     * Vector of length equal to the number of groups.
     * Each entry contains the buffers in which to store the corresponding group's summary statistics for the AUC.
     */
    std::vector<SummaryBuffers<Stat_, Rank_> > auc;

    /**
     * Vector of length equal to the number of groups.
     * Each entry contains the buffers in which to store the corresponding group's summary statistics for the delta-mean.AUC.
     */
    std::vector<SummaryBuffers<Stat_, Rank_> > delta_mean;

    /**
     * Vector of length equal to the number of groups.
     * Each entry contains the buffers in which to store the corresponding group's summary statistics for the delta-detected.
     */
    std::vector<SummaryBuffers<Stat_, Rank_> > delta_detected;
};

/**
 * @cond
 */
namespace internal {

enum class CacheAction : unsigned char { SKIP, COMPUTE, CACHE };

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
template<typename Stat_>
class EffectsCacher {
public:
    EffectsCacher(size_t ngenes, size_t ngroups, size_t cache_size) :
        my_ngenes(ngenes),
        my_ngroups(ngroups),
        my_cache_size(std::min(cache_size, ngroups * (ngroups - 1) / 2)), // cap it at the maximum possible number of comparisons.
        my_common_cache(ngenes * my_cache_size),
        my_actions(ngroups),
        my_staging_cache(ngroups)
    {
        my_unused_pool.reserve(my_cache_size);
        my_cached.reserve(my_cache_size);
        auto ptr = my_common_cache.data();
        for (size_t c = 0; c < my_cache_size; ++c, ptr += ngenes) {
            my_unused_pool.push_back(ptr);
        }
    }

private:
    size_t my_ngenes;
    size_t my_ngroups;
    size_t my_cache_size;

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
    std::map<std::pair<size_t, size_t>, Stat_*> my_cached;

public:
    void clear() {
        my_cached.clear();
    }

public:
    void fill_effects_from_cache(size_t group, std::vector<double>& full_effects) {
        // During calculation of effects, the current group (i.e., 'group') is
        // the first group in the comparison and 'other' is the second group.
        // However, remember that we want to cache the reverse effects, so in
        // the cached entry, 'group' is second and 'other' is first.
        for (size_t other = 0; other < my_ngroups; ++other) {
            if (other == group) {
                my_actions[other] = CacheAction::SKIP;
                continue;
            }

            if (cache_size == 0) {
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
            auto curcache = my_cached.begin()->second;
            size_t offset = other;
            for (size_t i = 0; i < my_ngenes; ++i, offset += my_ngroups) {
                full_effects[offset /* = other + ngroups * i */] = curcache[i];
            }

            my_unused_pool.push_back(x);
            my_cached.erase(cached.begin());
        }

        // Refining our choice of cacheable entries by doing a dummy run and
        // seeing whether eviction actually happens. If it doesn't, we won't
        // bother caching, because that would be a waste of memory accesses.
        for (size_t other = 0; other < ngroups; ++other) {
            if (my_actions[other] != CacheAction::CACHE) {
                continue;
            }

            std::pair<size_t, size_t> key(other, group);
            if (cached.size() < cache_size) {
                auto ptr = my_vector_pool.back();
                my_cached[key] = ptr;
                my_staging_cache[other] = ptr;
                my_vector_pool.pop_back();
                continue;
            }

            // Looking at the last cache entry. If the first group of this
            // entry is larger than the first group of the incoming entry, we
            // evict it, as the incoming entry has faster reusability.
            auto it = my_cached.end();
            --it;
            if ((it->first).first > other) {
                auto ptr = it->second;
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
    CacheAction get_action(size_t other) const {
        return my_actions[other];
    }

    Stat_* get_cache_location(size_t other) const {
        return my_staging_cache[other];
    }
};

template<typename Index_, typename Stat_, typename Rank_>
void process_simple_summary_effects(
    Index_ ngenes, 
    size_t ngroups,
    size_t nblocks,
    size_t ncombos,
    std::vector<Stat_>& combo_means,
    std::vector<Stat_>& combo_vars,
    std::vector<Stat_>& combo_detected,
    const ScoreMarkersSummaryBuffers<Stat_, Rank_>& output,
    const std::vector<Stat_>& combo_weights,
    double threshold,
    int num_threads)
{
    std::vector<Stat_> total_weights_per_group;
    const Stat_* total_weights_ptr = combo_weights.data();
    if (nblocks > 1) {
        total_weights_per_group = compute_total_weight_per_group(ngroups, nblocks, combo_weights.data());
        total_weights_ptr = total_weights_per_group.data();
    }
    PrecomputedPairwiseWeights<Stat_> preweights(ngroups, nblocks, combo_weights.data());

    EffectsCacher<Stat_> cache(ngenes, ngroups, cache_size);
    std::vector<Stat_> full_effects(ngroups * ngenes);
    std::vector<std::vector<Stat_> > effect_buffers(num_threads);
    for (auto& ef : effect_buffer) {
        ef.resize(ngroups);
    }

    if (output.cohens_d.size()) {
        cache.clear();
        for (size_t group = 0; group < ngroups; ++group) {
            cache.fill_effects_from_cache(group, full_effects);

            tatami::parallelize([&](size_t t, Index_ start, Index_ length) -> void {
                size_t in_offset = ncombos * static_cast<size_t>(start); // cast to avoid oerflow.
                auto my_means = tmp_means + in_offset;
                auto my_variances = tmp_variances + in_offset;

                auto store_ptr = full_effects.data() + static_cast<size_t>(start) * ngroups; // cast to avoid overflow.
                Index_ end = start + length;
                auto copy_store = store_ptr;

                auto cache_action = cache.get_action(other);
                auto cache_location = cache.get_cache_location(other);
                for (Index_ gene = start; gene < end; ++gene, my_means += ncombos, my_variances += ncombos, copy_store += ngroups) {
                    for (size_t other = 0; other < ngroups; ++other) {
                        if (cache_action == internal::CacheAction::COMPUTE) {
                            copy_store[other] = internal::compute_pairwise_cohens_d_one_sided(group, other, my_means, my_variances, ngroups, nblocks, preweights, threshold);
                        } else if (cache_action == internal::CacheAction::CACHE) {
                            auto tmp = differential_analysis::compute_pairwise_cohens_d_two_sided(group, other, my_means, my_variances, ngroups, nblocks, preweights, threshold);
                            copy_store[other] = tmp.first;
                            cache_location[gene] = tmp.second;
                        }
                    }
                }

                copy_store = store_ptr;
                for (Index_ gene = start; gene < end; ++gene, copy_store += ngroups) {
                    internal::summarize_comparisons(ngroups, copy_store, group, gene, output.cohens_d[group], effect_buffers[t]);
                }
            }, ngenes, num_threads);

            auto mr = output.cohens_d[group].min_rank;
            if (mr) {
                differential_analysis::compute_min_rank(ngenes, ngroups, group, full_effects.data(), mr, num_threads);
            }
        }
    }

    if (lfc.size()) {
        cache.clear();
        for (int group = 0; group < ngroups; ++group) {
            cache.configure(group, full_set);

            tatami::parallelize([&](size_t, size_t start, size_t length) -> void {
                auto my_means = tmp_means + nlevels * start;

                const auto& actions = cache.actions;
                auto& staging_cache = cache.staging_cache;

                auto lfc_ptr = full_set.data() + start * ngroups;
                std::vector<double> effect_buffer(ngroups);

                for (size_t gene = start, end = start + length; gene < end; ++gene, my_means += nlevels, lfc_ptr += ngroups) {
                    for (int other = 0; other < ngroups; ++other) {
                        if (actions[other] == differential_analysis::CacheAction::SKIP) {
                            continue;
                        }

                        auto val = differential_analysis::compute_pairwise_simple_diff(group, other, my_means, level_weight, ngroups, nblocks);
                        lfc_ptr[other] = val;
                        if (actions[other] == differential_analysis::CacheAction::CACHE) {
                            staging_cache[other][gene] = -val;
                        } 
                    }

                    differential_analysis::summarize_comparisons(ngroups, lfc_ptr, group, gene, lfc, effect_buffer);
                }
            }, ngenes, nthreads);

            if (lfc[differential_analysis::summary::MIN_RANK].size()) {
                differential_analysis::compute_min_rank(ngenes, ngroups, group, full_set.data(), lfc[differential_analysis::summary::MIN_RANK][group], nthreads);
            }

            cache.transfer(group);
        }
    }

    if (delta_detected.size()) {
        cache.clear();
        for (int group = 0; group < ngroups; ++group) {
            cache.configure(group, full_set);

            tatami::parallelize([&](size_t, size_t start, size_t length) -> void {
                auto my_detected = tmp_detected + nlevels * start;

                const auto& actions = cache.actions;
                auto& staging_cache = cache.staging_cache;

                auto delta_detected_ptr = full_set.data() + start * ngroups;
                std::vector<double> effect_buffer(ngroups);

                for (size_t gene = start, end = start + length; gene < end; ++gene, my_detected += nlevels, delta_detected_ptr += ngroups) {
                    for (int other = 0; other < ngroups; ++other) {
                        if (actions[other] == differential_analysis::CacheAction::SKIP) {
                            continue;
                        }

                        auto val = differential_analysis::compute_pairwise_simple_diff(group, other, my_detected, level_weight, ngroups, nblocks);
                        delta_detected_ptr[other] = val;
                        if (actions[other] == differential_analysis::CacheAction::CACHE) {
                            staging_cache[other][gene] = -val;
                        } 
                    }

                    differential_analysis::summarize_comparisons(ngroups, delta_detected_ptr, group, gene, delta_detected, effect_buffer);
                }
            }, ngenes, nthreads);

            if (delta_detected[differential_analysis::summary::MIN_RANK].size()) {
                differential_analysis::compute_min_rank(ngenes, ngroups, group, full_set.data(), delta_detected[differential_analysis::summary::MIN_RANK][group], nthreads);
            }

            cache.transfer(group);
        }
    }
}

template<typename Stat_>
void summarize_auc(
    size_t ngenes, 
    size_t ngroups,
    const differential_analysis::MatrixCalculator::State& state, 
    std::vector<std::vector<Stat_*> >& auc,
    std::vector<Stat_>& auc_buffer) 
const {
    // If we need the min-rank AUCs, there's no choice but to hold everything in memory.
    if (auc.size()) {
        differential_analysis::summarize_comparisons(ngenes, ngroups, auc_buffer.data(), auc, nthreads);
        if (auc[differential_analysis::summary::MIN_RANK].size()) {
            differential_analysis::compute_min_rank(ngenes, ngroups, auc_buffer.data(), auc[differential_analysis::summary::MIN_RANK], nthreads);
        }
    }
}

}
/**
 * @endcond
 */

/**
 * @brief Score each gene as a candidate marker for each group of cells.
 *
 * Markers are identified by differential expression analyses between pairs of groups of cells (e.g., clusters, cell types).
 * Given `n` groups, each group is involved in `n - 1` pairwise comparisons and thus has `n - 1` effect sizes.
 * For each group, we compute summary statistics - e.g., median, mean - of the effect sizes across all of that group's comparisons.
 * Users can then sort by any of these summaries to obtain a ranking of potential marker genes for each group.
 *
 * The choice of effect size and summary statistic determines the characteristics of the resulting marker set.
 * For the effect sizes: we compute Cohen's d, the area under the curve (AUC), the log-fold change and the delta-detected,
 * which are described in more detail in the documentation for `score_markers_pairwise()`.
 * For the summary statistics: we compute the minimum, mean, median, maximum and min-rank of the effect sizes across each group's pairwise comparisons,
 * which are described in `summarize_effects()`.
 *
 * If the dataset contains blocking factors such as batch or sample, we compute the effect size within each level of the blocking factor.
 * This avoids interference from batch effects or sample-to-sample variation.
 * Users can also adjust the effect size to account for a minimum log-fold change threshold,
 * in order to focus on markers with larger changes in expression.
 * See `PairwiseEffects` for more details. 
 */
template<typename Value_, typename Index_, typename Group_, typename Stat_>
void score_markers_summary(
    const tatami::Matrix<Value_, Index_>& matrix, 
    const Group_* group, 
    const ScoreMarkersPairwiseOptions& options,
    const ScoreMarkersPairwiseBuffers<Stat_>& output) 
{
    Index_ NC = matrix.ncol();
    auto group_sizes = tatami_stats::tabulate_groups(group, NC); 

}

public:
    /**
     * Score potential marker genes by computing summary statistics across pairwise comparisons between groups.
     * On completion, `means`, `detected`, `cohen`, `auc`, `lfc` and `delta_detected` are filled with their corresponding statistics.
     *
     * If `cohen` is of length 0, Cohen's d is not computed.
     * If any of the inner vectors are of length 0, the corresponding summary statistic is not computed.
     * The same applies to `auc`, `lfc` and `delta_detected`.
     * (`set_compute_cohen()` and related functions have no effect here.)
     *
     * @tparam Data_ Matrix data type.
     * @tparam Index_ Matrix index type.
     * @tparam Group_ Integer type for the group assignments.
     * @tparam Stat_ Floating-point type to store the statistics.
     *
     * @param p Pointer to a **tatami** matrix instance.
     * @param[in] group Pointer to an array of length equal to the number of columns in `p`, containing the group assignments.
     * These should be 0-based and consecutive.
     * @param[out] means Pointers to arrays of length equal to the number of rows in `p`,
     * used to store the mean expression of each group.
     * @param[out] detected Pointers to arrays of length equal to the number of rows in `p`,
     * used to store the proportion of detected expression in each group.
     * @param[out] cohen Vector of vector of pointers to arrays of length equal to the number of rows in `p`.
     * Each inner vector corresponds to a summary statistic for Cohen's d, ordered as in `differential_analysis::summary`.
     * Each pointer corresponds to a group and is filled with the relevant summary statistic for that group.
     * @param[out] auc Vector of vector of pointers as described for `cohen`, but instead storing summary statistics for the AUC.
     * @param[out] lfc Vector of vector of pointers as described for `cohen`, but instead storing summary statistics for the log-fold change instead of Cohen's d.
     * @param[out] delta_detected Vector of vector of pointers as described for `cohen`, but instead storing summary statistics for the delta in the detected proportions.
     */
    template<typename Value_, typename Index_, typename Group_, typename Stat_>
    void run(const tatami::Matrix<Value_, Index_>* p, const Group_* group, 
        std::vector<Stat_*> means, 
        std::vector<Stat_*> detected, 
        std::vector<std::vector<Stat_*> > cohen, 
        std::vector<std::vector<Stat_*> > auc,
        std::vector<std::vector<Stat_*> > lfc,
        std::vector<std::vector<Stat_*> > delta_detected) 
    const {
        differential_analysis::MatrixCalculator runner(nthreads, threshold, block_weight_policy, variable_block_weight_parameters);

        size_t ngenes = p->nrow();
        size_t ngroups = means.size();
        Overlord<Stat_> overlord(ngenes, ngroups, auc.empty());
        auto state = runner.run(p, group, ngroups, overlord);

        process_simple_effects(ngenes, ngroups, 1, state, means, detected, cohen, lfc, delta_detected);
        summarize_auc(ngenes, ngroups, state, auc, overlord.auc_buffer);
    }

    /**
     * Score potential marker genes by computing summary statistics across pairwise comparisons between groups in multiple blocks.
     * On completion, `means`, `detected`, `cohen`, `auc`, `lfc` and `delta_detected` are filled with their corresponding statistics.
     *
     * If `cohen` is of length 0, Cohen's d is not computed.
     * If any of the inner vectors are of length 0, the corresponding summary statistic is not computed.
     * The same applies to `auc`, `lfc` and `delta_detected`.
     * (`set_compute_cohen()` and related functions have no effect here.)
     *
     * @tparam Data_ Matrix data type.
     * @tparam Index_ Matrix index type.
     * @tparam Group_ Integer type for the group assignments.
     * @tparam Block_ Integer type for the block assignments.
     * @tparam Stat_ Floating-point type to store the statistics.
     *
     * @param p Pointer to a **tatami** matrix instance.
     * @param[in] group Pointer to an array of length equal to the number of columns in `p`, containing the group assignments.
     * These should be 0-based and consecutive.
     * @param[in] block Pointer to an array of length equal to the number of columns in `p`, containing the blocking factor.
     * Levels should be 0-based and consecutive.
     * If `NULL`, this is ignored and the result is the same as calling `run()`.
     * @param[out] means Vector of vectors of pointers to arrays of length equal to the number of rows in `p`.
     * Each inner vector corresponds to a group and each pointer therein contains the mean expression in a blocking level.
     * @param[out] detected Pointers to arrays of length equal to the number of rows in `p`.
     * Each inner vector corresponds to a group and each pointer therein contains the proportion of detected expression in a blocking level.
     * @param[out] cohen Vector of vector of pointers to arrays of length equal to the number of rows in `p`.
     * Each inner vector corresponds to a summary statistic for Cohen's d, ordered as in `differential_analysis::summary`.
     * Each pointer corresponds to a group and is filled with the relevant summary statistic for that group.
     * @param[out] auc Vector of vector of pointers as described for `cohen`, but instead storing summary statistics for the AUC.
     * @param[out] lfc Vector of vector of pointers as described for `cohen`, but instead storing summary statistics for the log-fold change instead of Cohen's d.
     * @param[out] delta_detected Vector of vector of pointers as described for `cohen`, but instead storing summary statistics for the delta in the detected proportions.
     */
    template<typename Value_, typename Index_, typename Group_, typename Block_, typename Stat_>
    void run_blocked(const tatami::Matrix<Value_, Index_>* p, const Group_* group, const Block_* block, 
        std::vector<Stat_*> means, 
        std::vector<Stat_*> detected, 
        std::vector<std::vector<Stat_*> > cohen,
        std::vector<std::vector<Stat_*> > auc,
        std::vector<std::vector<Stat_*> > lfc,
        std::vector<std::vector<Stat_*> > delta_detected) 
    const {
        if (block == NULL) {
            run(
                p, 
                group, 
                std::move(means),
                std::move(detected),
                std::move(cohen),
                std::move(auc),
                std::move(lfc),
                std::move(delta_detected)
            );
            return;
        }

        differential_analysis::MatrixCalculator runner(nthreads, threshold, block_weight_policy, variable_block_weight_parameters);

        size_t ngenes = p->nrow();
        size_t ngroups = means.size();
        size_t nblocks = count_ids(p->ncol(), block);
        Overlord<Stat_> overlord(ngenes, ngroups, auc.empty());
        auto state = runner.run_blocked(p, group, ngroups, block, nblocks, overlord);

        int ncombos = ngroups * nblocks;
        std::vector<std::vector<Stat_> > means_store(ncombos), detected_store(ncombos);
        std::vector<Stat_*> means2(ncombos), detected2(ncombos);
        for (int c = 0; c < ncombos; ++c) {
            means_store[c].resize(ngenes);
            detected_store[c].resize(ngenes);
            means2[c] = means_store[c].data();
            detected2[c] = detected_store[c].data();
        }

        process_simple_effects(ngenes, ngroups, nblocks, state, means2, detected2, cohen, lfc, delta_detected);
        summarize_auc(ngenes, ngroups, state, auc, overlord.auc_buffer);

        // Averaging the remaining statistics.
        std::vector<double> weights(nblocks);
        std::vector<Stat_*> mstats(nblocks), dstats(nblocks);

        for (int gr = 0; gr < ngroups; ++gr) {
            for (int b = 0; b < nblocks; ++b) {
                size_t offset = gr * static_cast<size_t>(nblocks) + b;
                weights[b] = state.level_weight[offset];
                mstats[b] = means2[offset];
                dstats[b] = detected2[offset];
            }

            average_vectors_weighted(ngenes, mstats, weights.data(), means[gr]);
            average_vectors_weighted(ngenes, dstats, weights.data(), detected[gr]);
        }
    }

private:
    template<typename Stat_>
    class Overlord {
    public:
        Overlord(size_t nr, size_t ng, bool skip_auc) : skipped(skip_auc), auc_buffer(skip_auc ? 0 : nr * ng * ng) {}

        bool needs_auc() const {
            return !skipped;
        }

        bool skipped;
        std::vector<Stat_> auc_buffer;

        Stat_* prepare_auc_buffer(size_t gene, size_t ngroups) { 
            return auc_buffer.data() + gene * ngroups * ngroups;
        }
    };

public:
    /** 
     * @brief Results of the marker scoring.
     * 
     * @tparam Stat_ Floating-point type to store the statistics.
     *
     * Meaningful instances of this object should generally be constructed by calling the `ScoreMarkers::run()` methods.
     * Empty instances can be default-constructed as placeholders.
     */
    template<typename Stat_>
    struct Results {
        /**
         * @cond
         */
        Results() {}

        Results(
            size_t ngenes, 
            int ngroups, 
            const ComputeSummaries& do_cohen, 
            const ComputeSummaries& do_auc, 
            const ComputeSummaries& do_lfc, 
            const ComputeSummaries& do_delta_detected)
        { 
            auto fill_inner = [&](int N, auto& type) {
                type.reserve(N);
                for (int n = 0; n < N; ++n) {
                    type.emplace_back(ngenes);
                }
            };

            fill_inner(ngroups, means);
            fill_inner(ngroups, detected);

            auto fill_effect = [&](const ComputeSummaries& do_this, auto& effect) {
                bool has_any = false;
                for (size_t i = 0; i < do_this.size(); ++i) {
                    if (do_this[i]) {
                        has_any = true;
                        break;
                    }
                }

                if (has_any) {
                    effect.resize(differential_analysis::n_summaries);
                    if (do_this[differential_analysis::MIN]) {
                        fill_inner(ngroups, effect[differential_analysis::MIN]);
                    }
                    if (do_this[differential_analysis::MEAN]) {
                        fill_inner(ngroups, effect[differential_analysis::MEAN]);
                    }
                    if (do_this[differential_analysis::MEDIAN]) {
                        fill_inner(ngroups, effect[differential_analysis::MEDIAN]);
                    }
                    if (do_this[differential_analysis::MAX]) {
                        fill_inner(ngroups, effect[differential_analysis::MAX]);
                    }
                    if (do_this[differential_analysis::MIN_RANK]) {
                        fill_inner(ngroups, effect[differential_analysis::MIN_RANK]);
                    }
                }
                return;
            };

            fill_effect(do_cohen, cohen);
            fill_effect(do_auc, auc);
            fill_effect(do_lfc, lfc);
            fill_effect(do_delta_detected, delta_detected);
            return;
        }
        /**
         * @endcond
         */

        /**
         * Summary statistics for Cohen's d.
         * Elements of the outer vector correspond to the different summary statistics (see `differential_analysis::summary`);
         * elements of the middle vector correspond to the different groups;
         * and elements of the inner vector correspond to individual genes.
         */
        std::vector<std::vector<std::vector<Stat_> > > cohen;

        /**
         * Summary statistics for the AUC.
         * Elements of the outer vector correspond to the different summary statistics (see `differential_analysis::summary`);
         * elements of the middle vector correspond to the different groups;
         * and elements of the inner vector correspond to individual genes.
         */
        std::vector<std::vector<std::vector<Stat_> > > auc;

        /**
         * Summary statistics for the log-fold change.
         * Elements of the outer vector correspond to the different summary statistics (see `differential_analysis::summary`);
         * elements of the middle vector correspond to the different groups;
         * and elements of the inner vector correspond to individual genes.
         */
        std::vector<std::vector<std::vector<Stat_> > > lfc;

        /**
         * Summary statistics for the delta in the detected proportions.
         * Elements of the outer vector correspond to the different summary statistics (see `differential_analysis::summary`);
         * elements of the middle vector correspond to the different groups;
         * and elements of the inner vector correspond to individual genes.
         */
        std::vector<std::vector<std::vector<Stat_> > > delta_detected;

        /**
         * Mean expression in each group.
         * Elements of the outer vector corresponds to the different groups, and elements of the inner vector correspond to individual genes.
         */
        std::vector<std::vector<Stat_> > means;

        /**
         * Proportion of detected expression in each group.
         * Elements of the outer vector corresponds to the different groups, and elements of the inner vector correspond to individual genes.
         */
        std::vector<std::vector<Stat_> > detected;
    };

    /**
     * Score potential marker genes by computing summary statistics across pairwise comparisons between groups. 
     *
     * @tparam Stat_ Floating-point type to store the statistics.
     * @tparam Data_ Matrix data type.
     * @tparam Index_ Matrix index type.
     * @tparam Group_ Integer type for the group assignments.
     *
     * @param p Pointer to a **tatami** matrix instance.
     * @param[in] group Pointer to an array of length equal to the number of columns in `p`, containing the group assignments.
     * These should be 0-based and consecutive.
     *
     * @return A `Results` object containing the summary statistics and the other per-group statistics.
     * Whether particular statistics are computed depends on the configuration from `set_compute_cohen()` and related setters.
     */
    template<typename Stat_ = double, typename Value_, typename Index_, typename Group_>
    Results<Stat_> run(const tatami::Matrix<Value_, Index_>* p, const Group_* group) const {
        auto ngroups = count_ids(p->ncol(), group);
        Results<Stat_> res(p->nrow(), ngroups, do_cohen, do_auc, do_lfc, do_delta_detected); 
        run(
            p, 
            group,
            vector_to_pointers(res.means),
            vector_to_pointers(res.detected),
            vector_to_pointers(res.cohen),
            vector_to_pointers(res.auc),
            vector_to_pointers(res.lfc),
            vector_to_pointers(res.delta_detected)
        );
        return res;
    }

    /**
     * Score potential marker genes by computing summary statistics across pairwise comparisons between groups in multiple blocks.
     *
     * @tparam Stat_ Floating-point type to store the statistics.
     * @tparam Data_ Matrix data type.
     * @tparam Index_ Matrix index type.
     * @tparam Group_ Integer type for the group assignments.
     * @tparam Block_ Integer type for the block assignments.
     *
     * @param p Pointer to a **tatami** matrix instance.
     * @param[in] group Pointer to an array of length equal to the number of columns in `p`, containing the group assignments.
     * These should be 0-based and consecutive.
     * @param[in] block Pointer to an array of length equal to the number of columns in `p`, containing the blocking factor.
     * Levels should be 0-based and consecutive.
     * If `NULL`, this is ignored and the result is the same as calling `run()`.
     *
     * @return A `Results` object containing the summary statistics and the other per-group statistics.
     * Whether particular statistics are computed depends on the configuration from `set_compute_cohen()` and related setters.
     */
    template<typename Stat_ = double, typename Value_, typename Index_, typename Group_, typename Block_>
    Results<Stat_> run_blocked(const tatami::Matrix<Value_, Index_>* p, const Group_* group, const Block_* block) const {
        auto ngroups = count_ids(p->ncol(), group);
        Results<Stat_> res(p->nrow(), ngroups, do_cohen, do_auc, do_lfc, do_delta_detected); 
        run_blocked(
            p, 
            group,
            block,
            vector_to_pointers(res.means),
            vector_to_pointers(res.detected),
            vector_to_pointers(res.cohen),
            vector_to_pointers(res.auc),
            vector_to_pointers(res.lfc),
            vector_to_pointers(res.delta_detected)
        );
        return res;
    }
};

}

#endif
