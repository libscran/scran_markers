#ifndef SCRAN_SCORE_MARKERS_HPP
#define SCRAN_SCORE_MARKERS_HPP

#include "scran_blocks/scran_blocks.hpp"
#include "tatami/tatami.hpp"
#include "sanisizer/sanisizer.hpp"
#include "topicks/topicks.hpp"

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
     * Whether to compute the mean expression in each group.
     * This only affects the `score_markers_summary()` overload that returns a `ScoreMarkersSummaryResults`.
     */
    bool compute_group_mean = true;

    /**
     * Whether to compute the proportion of cells with detected expression in each group.
     * This only affects the `score_markers_summary()` overload that returns a `ScoreMarkersSummaryResults`.
     */
    bool compute_group_detected = true;

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
     * Limit on the reported minimum rank.
     * If a gene has a minimum rank greater than `min_rank_limit`, its reported minimum rank will be set to the number of genes.
     * Lower values improve memory efficiency but discard all ranking information beyond `min_rank_limit`. 
     */
    std::size_t min_rank_limit = 500;

    /**
     * Whether to preserve ties when computing the minimum rank.
     * Tied genes with equal effect sizes receive the same rank within each pairwise comparison. 
     */
    bool min_rank_preserve_ties = false; // TODO: set this to TRUE.

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
     *
     * Alternatively, this vector may be empty, in which case the means are not computed.
     */
    std::vector<Stat_*> mean;

    /**
     * Vector of length equal to the number of groups.
     * Each pointer corresponds to a group and points to an array of length equal to the number of genes,
     * to be filled with the proportion of cells with detected expression in that group. 
     *
     * Alternatively, this vector may be empty, in which case the detected proportions are not computed.
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

template<typename Stat_, typename Index_>
using MinrankTopQueues = std::vector<std::vector<topicks::TopQueue<Stat_, Index_> > >;

template<typename Stat_, typename Index_, typename Rank_>
void preallocate_minrank_queues( 
    const std::size_t ngroups,
    std::vector<MinrankTopQueues<Stat_, Index_> >& all_queues,
    const std::vector<SummaryBuffers<Stat_, Rank_> >& summaries,
    const Index_ limit,
    const bool keep_ties,
    const int num_threads
) { 
    topicks::TopQueueOptions<Stat_> qopt;
    qopt.keep_ties = keep_ties;
    qopt.check_nan = true;

    sanisizer::resize(all_queues, num_threads);
    for (int t = 0; t < num_threads; ++t) {
        sanisizer::resize(all_queues[t], ngroups);

        for (decltype(I(ngroups)) g1 = 0; g1 < ngroups; ++g1) {
            if (summaries[g1].min_rank == NULL) {
                continue;
            }

            for (decltype(I(ngroups)) g2 = 0; g2 < ngroups; ++g2) {
                all_queues[t][g1].emplace_back(limit, true, qopt);
            }
        }
    }
}

template<typename Stat_, typename Index_, typename Rank_>
void compute_summary_stats_per_gene(
    const Index_ gene,
    const std::size_t ngroups,
    const Stat_* const pairwise_buffer_ptr,
    std::vector<Stat_>& summary_buffer,
    MinrankTopQueues<Stat_, Index_>& minrank_queues,
    const std::vector<SummaryBuffers<Stat_, Rank_> >& summaries
) {
    for (decltype(I(ngroups)) gr = 0; gr < ngroups; ++gr) {
        auto& cursummary = summaries[gr];
        const auto in_offset = sanisizer::product_unsafe<std::size_t>(ngroups, gr);
        summarize_comparisons(ngroups, pairwise_buffer_ptr + in_offset, gr, gene, cursummary, summary_buffer);

        if (cursummary.min_rank) {
            auto& cur_queues = minrank_queues[gr];
            for (decltype(I(ngroups)) gr2 = 0; gr2 < ngroups; ++gr2) {
                if (gr != gr2) {
                    cur_queues[gr2].emplace(pairwise_buffer_ptr[in_offset + gr2], gene);
                }
            }
        }
    }
}

template<typename Stat_, typename Index_, typename Rank_>
void report_minrank_from_queues(
    const Index_ ngenes,
    const std::size_t ngroups,
    std::vector<MinrankTopQueues<Stat_, Index_> >& all_queues,
    const std::vector<SummaryBuffers<Stat_, Rank_> >& summaries,
    const int num_threads,
    const bool keep_ties
) {
    tatami::parallelize([&](const int, const std::size_t start, const std::size_t length) -> void {
        std::vector<Index_> tie_buffer;

        for (decltype(I(ngroups)) gr = start, grend = start + length; gr < grend; ++gr) {
            const auto mr_out = summaries[gr].min_rank;
            if (mr_out == NULL) {
                continue;
            }

            // Using the maximum possible rank (i.e., 'ngenes') as the default.
            const auto maxrank_placeholder = sanisizer::cast<Rank_>(ngenes);
            std::fill_n(mr_out, ngenes, maxrank_placeholder);

            for (decltype(I(ngroups)) gr2 = 0; gr2 < ngroups; ++gr2) {
                if (gr == gr2) {
                    continue;
                }
                auto& current_out = all_queues.front()[gr][gr2];

                const auto num_queues = all_queues.size();
                for (decltype(I(num_queues)) q = 1; q < num_queues; ++q) {
                    auto& current_in = all_queues[q][gr][gr2];
                    while (!current_in.empty()) {
                        current_out.push(current_in.top());
                        current_in.pop();
                    }
                }

                // Cast to Rank_ is safe as current_out.size() <= ngenes,
                // and we already checked that ngenes can fit into Rank_ in report_minrank_from_current_outs().
                if (!keep_ties) {
                    while (!current_out.empty()) {
                        auto& mr = mr_out[current_out.top().second];
                        mr = std::min(mr, static_cast<Rank_>(current_out.size()));
                        current_out.pop();
                    }
                } else {
                    while (!current_out.empty()) {
                        tie_buffer.clear();
                        const auto curtop = current_out.top();
                        current_out.pop();

                        while (!current_out.empty() && current_out.top().first == curtop.first) {
                            tie_buffer.push_back(current_out.top().second);
                            current_out.pop();
                        }

                        // Increment is safe as we already reduced the size at least once.
                        const Rank_ tied_rank = current_out.size() + 1;

                        mr_out[curtop.second] = std::min(mr_out[curtop.second], tied_rank);
                        for (const auto t : tie_buffer) {
                            mr_out[t] = std::min(mr_out[t], tied_rank);
                        }
                    }
                }
            }
        }
    }, ngroups, num_threads);
}

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
    const Index_ minrank_limit,
    const bool minrank_keep_ties,
    const int num_threads)
{
    std::vector<MinrankTopQueues<Stat_, Index_> > cohens_d_minrank_all_queues, delta_mean_minrank_all_queues, delta_detected_minrank_all_queues;
    if (output.cohens_d.size()) {
        preallocate_minrank_queues(ngroups, cohens_d_minrank_all_queues, output.cohens_d, minrank_limit, minrank_keep_ties, num_threads);
    }
    if (output.delta_mean.size()) {
        preallocate_minrank_queues(ngroups, delta_mean_minrank_all_queues, output.delta_mean, minrank_limit, minrank_keep_ties, num_threads);
    }
    if (output.delta_detected.size()) {
        preallocate_minrank_queues(ngroups, delta_detected_minrank_all_queues, output.delta_detected, minrank_limit, minrank_keep_ties, num_threads);
    }

    std::vector<Stat_> total_weights_per_group;
    auto total_weights_ptr = combo_weights.data();
    if (!output.mean.empty() || !output.detected.empty()) {
        if (nblocks > 1) {
            total_weights_per_group = compute_total_weight_per_group(ngroups, nblocks, combo_weights.data());
            total_weights_ptr = total_weights_per_group.data();
        }
    }

    PrecomputedPairwiseWeights<Stat_> preweights;
    if (!output.cohens_d.empty() || !output.delta_mean.empty() || !output.delta_detected.empty()) {
        preweights =  PrecomputedPairwiseWeights<Stat_>(ngroups, nblocks, combo_weights.data());
    }

    const auto ngroups2 = sanisizer::product<typename std::vector<Stat_>::size_type>(ngroups, ngroups);
    tatami::parallelize([&](const int t, const Index_ start, const Index_ length) -> void {
        std::vector<Stat_> pairwise_buffer(ngroups2);
        std::vector<Stat_> summary_buffer(ngroups);

        for (Index_ gene = start, end = start + length; gene < end; ++gene) {
            const auto in_offset = sanisizer::product_unsafe<std::size_t>(gene, ncombos);

            if (!output.mean.empty()) {
                const auto tmp_means = combo_means.data() + in_offset;
                average_group_stats(gene, ngroups, nblocks, tmp_means, combo_weights.data(), total_weights_ptr, output.mean);
            }

            if (!output.detected.empty()) {
                const auto tmp_detected = combo_detected.data() + in_offset;
                average_group_stats(gene, ngroups, nblocks, tmp_detected, combo_weights.data(), total_weights_ptr, output.detected);
            }

            if (output.cohens_d.size()) {
                const auto tmp_means = combo_means.data() + in_offset;
                const auto tmp_variances = combo_vars.data() + in_offset;
                compute_pairwise_cohens_d(tmp_means, tmp_variances, ngroups, nblocks, preweights, threshold, pairwise_buffer.data());
                compute_summary_stats_per_gene(gene, ngroups, pairwise_buffer.data(), summary_buffer, cohens_d_minrank_all_queues[t], output.cohens_d);
            }

            if (output.delta_mean.size()) {
                const auto tmp_means = combo_means.data() + in_offset;
                compute_pairwise_simple_diff(tmp_means, ngroups, nblocks, preweights, pairwise_buffer.data());
                compute_summary_stats_per_gene(gene, ngroups, pairwise_buffer.data(), summary_buffer, delta_mean_minrank_all_queues[t], output.delta_mean);
            }

            if (output.delta_detected.size()) {
                const auto tmp_det = combo_detected.data() + in_offset;
                compute_pairwise_simple_diff(tmp_det, ngroups, nblocks, preweights, pairwise_buffer.data());
                compute_summary_stats_per_gene(gene, ngroups, pairwise_buffer.data(), summary_buffer, delta_detected_minrank_all_queues[t], output.delta_detected);
            }
        }
    }, ngenes, num_threads);

    if (output.cohens_d.size()) {
        report_minrank_from_queues(ngenes, ngroups, cohens_d_minrank_all_queues, output.cohens_d, num_threads, minrank_keep_ties);
    }

    if (output.delta_mean.size()) {
        report_minrank_from_queues(ngenes, ngroups, delta_mean_minrank_all_queues, output.delta_mean, num_threads, minrank_keep_ties);
    }

    if (output.delta_detected.size()) {
        report_minrank_from_queues(ngenes, ngroups, delta_detected_minrank_all_queues, output.delta_detected, num_threads, minrank_keep_ties);
    }
}

template<typename Index_, typename Stat_, typename Rank_>
ScoreMarkersSummaryBuffers<Stat_, Rank_> preallocate_summary_results(
    const Index_ ngenes,
    const std::size_t ngroups,
    ScoreMarkersSummaryResults<Stat_, Rank_>& store,
    const ScoreMarkersSummaryOptions& options)
{
    ScoreMarkersSummaryBuffers<Stat_, Rank_> output;

    if (options.compute_group_mean) { 
        preallocate_average_results(ngenes, ngroups, store.mean, output.mean);
    }

    if (options.compute_group_detected) { 
        preallocate_average_results(ngenes, ngroups, store.detected, output.detected);
    }

    if (options.compute_cohens_d) {
        output.cohens_d = fill_summary_results(
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
        output.auc = fill_summary_results(
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
        output.delta_mean = fill_summary_results(
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
        output.delta_detected = fill_summary_results(
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
    std::vector<Stat_> combo_means, combo_vars, combo_detected;
    if (!output.mean.empty() || !output.cohens_d.empty() || !output.delta_mean.empty()) {
        combo_means.resize(payload_size);
    }
    if (!output.cohens_d.empty()) {
        combo_vars.resize(payload_size);
    }
    if (!output.detected.empty() || !output.delta_detected.empty()) {
        combo_detected.resize(payload_size);
    }

    // For a single block, this usually doesn't really matter, but we do it for consistency with the multi-block case,
    // and to account for variable weighting where non-zero block sizes get zero weight.
    const auto combo_weights = scran_blocks::compute_weights<Stat_>(
        combo_sizes,
        options.block_weight_policy,
        options.variable_block_weight_parameters
    );

    const Index_ minrank_limit = sanisizer::cap<Index_>(options.min_rank_limit);

    if (!output.auc.empty()) {
        std::vector<MinrankTopQueues<Stat_, Index_> > auc_minrank_all_queues;
        preallocate_minrank_queues(ngroups, auc_minrank_all_queues, output.auc, minrank_limit, options.min_rank_preserve_ties, options.num_threads);

        struct AucResultWorkspace {
            AucResultWorkspace() = default;
            AucResultWorkspace(const std::size_t ngroups, MinrankTopQueues<Stat_, Index_>& queues) :
                pairwise_buffer(sanisizer::product<typename std::vector<Stat_>::size_type>(ngroups, ngroups)),
                summary_buffer(sanisizer::cast<typename std::vector<Stat_>::size_type>(ngroups)),
                queue_ptr(&queues)
            {};
            std::vector<Stat_> pairwise_buffer;
            std::vector<Stat_> summary_buffer;
            MinrankTopQueues<Stat_, Index_>* queue_ptr;
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
            /* auc_result_initialize = */ [&](const int t) -> AucResultWorkspace {
                return AucResultWorkspace(ngroups, auc_minrank_all_queues[t]);
            },
            /* auc_result_process = */ [&](
                const Index_ gene,
                AucScanWorkspace<Value_, Group_, Index_, Stat_>& auc_work,
                AucResultWorkspace& res_work
            ) -> void {
                process_auc_for_rows(auc_work, ngroups, nblocks, options.threshold, res_work.pairwise_buffer.data());
                compute_summary_stats_per_gene(gene, ngroups, res_work.pairwise_buffer.data(), res_work.summary_buffer, *(res_work.queue_ptr), output.auc);
            },
            combo_sizes,
            combo_weights,
            options.num_threads
        );

        report_minrank_from_queues(ngenes, ngroups, auc_minrank_all_queues, output.auc, options.num_threads, options.min_rank_preserve_ties);

    } else if (matrix.prefer_rows()) {
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
            static_cast<Stat_*>(NULL),
            combo_sizes,
            combo_weights, 
            options.threshold,
            options.num_threads
        );

    } else {
        scan_matrix_by_column(
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

    process_simple_summary_effects(
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
        minrank_limit,
        options.min_rank_preserve_ties,
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
    const auto buffers = internal::preallocate_summary_results(matrix.nrow(), ngroups, output, options);
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
    const auto buffers = internal::preallocate_summary_results(matrix.nrow(), ngroups, output, options);
    score_markers_summary_blocked(matrix, group, block, options, buffers);
    return output;
}

}

#endif
