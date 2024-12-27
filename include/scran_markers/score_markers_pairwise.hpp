#ifndef SCRAN_MARKERS_SCORE_MARKERS_PAIRWISE_HPP
#define SCRAN_MARKERS_SCORE_MARKERS_PAIRWISE_HPP

#include "cohens_d.hpp"
#include "simple_diff.hpp"

#include "scan_matrix.hpp"
#include "average_group_stats.hpp"
#include "PrecomputedPairwiseWeights.hpp"
#include "create_combinations.hpp"

#include <vector>

#include "scran_blocks/scran_blocks.hpp"
#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"

/**
 * @file score_markers_pairwise.hpp
 * @brief Score potential markers by pairwise effect sizes between groups of cells.
 */

namespace scran_markers {

/**
 * @brief Options for `score_markers_pairwise()` and friends.
 */
struct ScoreMarkersPairwiseOptions {
    /**
     * Threshold on the differences in expression values, used to adjust the Cohen's D and AUC calculations.
     * This should be non-negative.
     */
    double threshold = 0;

    /**
     * Number of threads to use. 
     * The parallelization scheme is determined by `tatami::parallelize()`.
     */
    int num_threads = 1;

    /**
     * Whether to compute Cohen's d. 
     * This only affects the `score_markers_pairwise()` overload that returns a `ScoreMarkersPairwiseResults`.
     */
    bool compute_cohens_d = true;

    /**
     * Whether to compute the AUC.
     * This only affects the `score_markers_pairwise()` overload that returns a `ScoreMarkersPairwiseResults`.
     */
    bool compute_auc = true;

    /**
     * Whether to compute the difference in means.
     * This only affects the `score_markers_pairwise()` overload that returns a `ScoreMarkersPairwiseResults`.
     */
    bool compute_delta_mean = true;

    /**
     * Whether to compute the difference in the detected proportion.
     * This only affects the `score_markers_pairwise()` overload that returns a `ScoreMarkersPairwiseResults`.
     */
    bool compute_delta_detected = true;

    /**
     * Policy to use for weighting blocks when computing average statistics/effect sizes across blocks.
     */
    scran_blocks::WeightPolicy block_weight_policy = scran_blocks::WeightPolicy::VARIABLE;

    /**
     * Parameters for the variable block weights.
     * Only used when `ScoreMarkersPairwiseOptions::block_weight_policy = scran_blocks::WeightPolicy::VARIABLE`.
     */
    scran_blocks::VariableWeightParameters variable_block_weight_parameters;
};

/**
 * @brief Buffers for `score_markers_pairwise()` and friends.
 * @tparam Stat_ Floating-point type for the output statistics.
 */
template<typename Stat_>
struct ScoreMarkersPairwiseBuffers {
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
     * Pointer to an array of length equal to \f$GN^2\f$, where \f$G\f$ is the number of genes and \f$N\f$ is the number of groups.
     * This is a 3-dimensional \f$G \times N \times N\f$ array to be filled with the Cohen's D for the comparison between each pair of groups for each gene.
     *
     * The first dimension is the slowest changing, is of length equal to the number of genes, and represents the gene.
     * The second dimension is the second-fastest changing, is of length equal to the number of groups, and represents the first group.
     * The third dimension is the fastest changing, is also of length equal to the number of groups, and represents the second group.
     *
     * Thus, the entry \f$(i, j, k)\f$ (i.e., `effects[i * N * N + j * N + k]`) represents the effect size of gene \f$i\f$ upon comparing group \f$j\f$ against group \f$k\f$.
     * Positive values represent upregulation in group \f$j\f$ compared to \f$k\f$.
     * Note that the comparison of each group to itself is always assigned an effect size of zero, regardless of the `threshold` used in `score_markers_pairwise()`;
     * this is only done to avoid exposing uninitialized memory, and the value should be ignored in downstream steps.
     *
     * Alternatively NULL, in which case the Cohen's D is not stored.
     */
    Stat_* cohens_d = NULL;

    /**
     * Pointer to an array of length equal to \f$GN^2\f$, where \f$G\f$ is the number of genes and \f$N\f$ is the number of groups.
     * This is a 3-dimensional array to be filled with the AUC for the comparison between each pair of groups for each gene;
     * see `ScoreMarkersPairwiseBuffers::cohens_d` for more details.
     *
     * Unlike Cohen's D, all AUC values will lie in \f$[0, 1]\f$.
     * Values above 0.5 represent upregulation in group \f$j\f$ compared to \f$k\f$.
     * The exception to this logic is the comparison of each group to itself, which is always assigned an effect size of zero instead of 0.5;
     * this is only done to avoid exposing uninitialized memory, and the value should be ignored in downstream steps.
     *
     * Alternatively NULL, in which case the AUC is not stored.
     */
    Stat_* auc = NULL;

    /**
     * Pointer to an array of length equal to \f$GN^2\f$, where \f$G\f$ is the number of genes and \f$N\f$ is the number of groups.
     * This is a 3-dimensional array to be filled with the difference in means for the comparison between each pair of groups for each gene;
     * see `ScoreMarkersPairwiseBuffers::cohens_d` for more details.
     * Alternatively NULL, in which case the difference in means is not stored.
     */
    Stat_* delta_mean = NULL;

    /**
     * Pointer to an array of length equal to \f$GN^2\f$, where \f$G\f$ is the number of genes and \f$N\f$ is the number of groups.
     * This is a 3-dimensional array to be filled with the difference in the detected proportions for the comparison between each pair of groups for each gene;
     * see `ScoreMarkersPairwiseBuffers::cohens_d` for more details.
     * Alternatively NULL, in which case the difference in detected proportions is not stored.
     */
    Stat_* delta_detected = NULL;
};

/**
 * @brief Results for `score_markers_pairwise()` and friends.
 * @tparam Stat_ Floating-point type for the output statistics.
 */
template<typename Stat_>
struct ScoreMarkersPairwiseResults {
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
     * Vector of length equal to \f$GN^2\f$, where \f$N\f$ is the number of genes and \f$G\f$ is the number of groups.
     * This is a 3-dimensional array to be filled with the Cohen's D for the comparison between each pair of groups for each gene;
     * see `ScoreMarkersPairwiseBuffers::cohens_d` for details on the layout.
     * Alternatively this may be an empty vector if `ScoreMarkersPairwiseOptions::compute_cohens_d = false`.
     */
    std::vector<Stat_> cohens_d;

    /**
     * Vector of length equal to \f$GN^2\f$, where \f$N\f$ is the number of genes and \f$G\f$ is the number of groups.
     * This is a 3-dimensional array to be filled with the AUC for the comparison between each pair of groups for each gene;
     * see `ScoreMarkersPairwiseBuffers::auc` for details on the layout.
     * Alternatively this may be an empty vector if `ScoreMarkersPairwiseOptions::compute_auc = false`.
     */
    std::vector<Stat_> auc;

    /**
     * Vector of length equal to \f$GN^2\f$, where \f$N\f$ is the number of genes and \f$G\f$ is the number of groups.
     * This is a 3-dimensional array to be filled with the delta-mean for the comparison between each pair of groups for each gene;
     * see `ScoreMarkersPairwiseBuffers::cohens_d` for details on the layout.
     * Alternatively this may be an empty vector if `ScoreMarkersPairwiseOptions::compute_delta_mean = false`.
     */
    std::vector<Stat_> delta_mean;

    /**
     * Vector of length equal to \f$GN^2\f$, where \f$N\f$ is the number of genes and \f$G\f$ is the number of groups.
     * This is a 3-dimensional array to be filled with the delta-detected for the comparison between each pair of groups for each gene;
     * see `ScoreMarkersPairwiseBuffers::cohens_d` for details on the layout.
     * Alternatively this may be an empty vector if `ScoreMarkersPairwiseOptions::compute_delta_detected = false`.
     */
    std::vector<Stat_> delta_detected;
};

/**
 * @cond
 */
namespace internal {

template<typename Index_, typename Stat_>
void process_simple_pairwise_effects(
    Index_ ngenes,
    size_t ngroups,
    size_t nblocks,
    size_t ncombos,
    std::vector<Stat_>& combo_means,
    std::vector<Stat_>& combo_vars,
    std::vector<Stat_>& combo_detected,
    const ScoreMarkersPairwiseBuffers<Stat_>& output,
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

    tatami::parallelize([&](size_t, Index_ start, Index_ length) -> void {
        size_t in_offset = ncombos * static_cast<size_t>(start);
        const auto* tmp_means = combo_means.data() + in_offset;
        const auto* tmp_variances = combo_vars.data() + in_offset;
        const auto* tmp_detected = combo_detected.data() + in_offset;

        size_t squared = ngroups * ngroups;
        for (size_t gene = start, end = start + length; gene < end; ++gene) {
            average_group_stats(gene, ngroups, nblocks, tmp_means, tmp_detected, combo_weights.data(), total_weights_ptr, output.mean, output.detected);

            // Computing the effect sizes.
            size_t out_offset = gene * squared;
            if (output.cohens_d != NULL) {
                internal::compute_pairwise_cohens_d(tmp_means, tmp_variances, ngroups, nblocks, preweights, threshold, output.cohens_d + out_offset);
            }

            if (output.delta_detected != NULL) {
                internal::compute_pairwise_simple_diff(tmp_detected, ngroups, nblocks, preweights, output.delta_detected + out_offset);
            }

            if (output.delta_mean != NULL) {
                internal::compute_pairwise_simple_diff(tmp_means, ngroups, nblocks, preweights, output.delta_mean + out_offset);
            }

            tmp_means += ncombos;
            tmp_variances += ncombos;
            tmp_detected += ncombos;
        }
    }, ngenes, num_threads);
}

template<typename Stat_>
ScoreMarkersPairwiseBuffers<Stat_> fill_pairwise_results(size_t ngenes, size_t ngroups, ScoreMarkersPairwiseResults<Stat_>& store, const ScoreMarkersPairwiseOptions& opt) {
    ScoreMarkersPairwiseBuffers<Stat_> output;

    internal::fill_average_results(ngenes, ngroups, store.mean, store.detected, output.mean, output.detected);

    size_t num_effect_sizes = ngenes * ngroups * ngroups; // already size_t's, no need to cast.

    if (opt.compute_cohens_d) {
        store.cohens_d.resize(num_effect_sizes);
        output.cohens_d = store.cohens_d.data();
    }
    if (opt.compute_auc) {
        store.auc.resize(num_effect_sizes);
        output.auc = store.auc.data();
    }
    if (opt.compute_delta_mean) {
        store.delta_mean.resize(num_effect_sizes);
        output.delta_mean = store.delta_mean.data();
    }
    if (opt.compute_delta_detected) {
        store.delta_detected.resize(num_effect_sizes);
        output.delta_detected = store.delta_detected.data();
    }

    return output;
}

}
/**
 * @endcond
 */

/**
 * Compute the effect sizes for the pairwise comparisons between groups.
 * This can be used to identify marker genes based on a specific comparison between two groups of interest.
 * Alternatively, the pairwise effects can be passed to `summarize_effects()` to obtain summaries for each group
 * (though it would be more efficient to use `score_markers_summary() to do so).
 *
 * @section effect-sizes Choice of effect sizes
 * The delta-mean is the difference in the mean expression between groups.
 * This is fairly straightforward to interpret, where a positive delta-mean corresponds to increased expression in the first group compared to the second. 
 * The delta-mean can also be treated as the log-fold change if the input matrix contains log-transformed normalized expression values.
 *
 * The delta-detected is the difference in the proportion of cells with detected expression between groups.
 * This lies between 1 and -1, with the extremes occurring when a gene is silent in one group and detected in all cells of the other group.
 * For this interpretation, we assume that the input matrix contains non-negative expression values, where a value of zero corresponds to lack of detectable expression.
 *
 * Cohen's d is the standardized difference between two groups.
 * This is defined as the difference in the mean for each group scaled by the average standard deviation across the two groups.
 * (Technically, we should use the pooled variance; however, this introduces some unintuitive asymmetry depending on the variance of the larger group, so we take a simple average instead.)
 * A positive value indicates that the gene has increased expression in the first group compared to the second.
 * Cohen's d is analogous to the t-statistic in a two-sample t-test and avoids spuriously large effect sizes from comparisons between highly variable groups.
 * We can also interpret Cohen's d as the number of standard deviations between the two group means.
 *
 * The area under the curve (AUC) can be interpreted as the probability that a randomly chosen observation in one group is greater than a randomly chosen observation in the other group. 
 * Values greater than 0.5 indicate that a gene is upregulated in the first group.
 * The AUC is closely related to the U-statistic used in the Wilcoxon rank sum test. 
 * The key difference between the AUC and Cohen's d is that the former is less sensitive to the variance within each group, e.g.,
 * if two distributions exhibit no overlap, the AUC is the same regardless of the variance of each distribution. 
 * This may or may not be desirable as it improves robustness to outliers but reduces the information available to obtain a highly resolved ranking. 
 *
 * @section threshold With a minimum change threshold
 * Setting a minimum change threshold (see `ScoreMarkersPairwiseOptions::threshold`) can be helpful as it prioritizes genes with large shifts in expression instead of those with low variances.
 * Currently, only positive thresholds are supported - this focuses on genes that are upregulated in the first group compared to the second.
 * The effect size definitions are generalized when testing against a non-zero threshold.
 *
 * - Cohen's d is redefined as the standardized difference between the difference in means and the specified threshold, analogous to the TREAT method from **limma**.
 *   Large positive values are only obtained when the observed difference in means is significantly greater than the threshold.
 *   For example, if we had a threshold of 2 and we obtained a Cohen's d of 3, this means that the observed difference in means was 3 standard deviations greater than 2.
 *   Importantly, a negative Cohen's d cannot be intepreted as downregulation, as the difference may still be positive but less than the threshold.
 * - The AUC is generalized to the probability of obtaining a random observation in one group that is greater than a random observation plus the threshold in the other group.
 *   For example, if we had a threshold of 2 and we obtained an AUC of 0.8, this means that - 80% of the time - 
 *   the random observation from the first group would be greater than a random observation from the second group by 2 or more.
 *   Again, AUCs below 0.5 cannot be interpreted as downregulation, as it may be caused by a positive shift that is less than the threshold.
 * 
 * @section other Other statistics
 * We report the mean expression of all cells in each group, as well as the proportion of cells with detectable expression in each group.
 * These statistics are useful for quickly interpreting the differences in expression driving the effect size summaries.
 *
 * @tparam Value_ Matrix data type.
 * @tparam Index_ Matrix index type.
 * @tparam Group_ Integer type for the group assignments.
 * @tparam Stat_ Floating-point type to store the statistics.
 *
 * @param matrix A **tatami** matrix instance.
 * @param[in] group Pointer to an array of length equal to the number of columns in `matrix`, containing the group assignments.
 * Group identifiers should be 0-based and should contain all integers in \f$[0, N)\f$ where \f$N\f$ is the number of unique groups.
 * @param options Further options.
 * @param[out] output Collection of buffers in which to store the computed statistics.
 * Each buffer is filled with the corresponding statistic for each group or pairwise comparison.
 * Any of `ScoreMarkersPairwiseBuffers::cohens_d`, 
 * `ScoreMarkersPairwiseBuffers::auc`, 
 * `ScoreMarkersPairwiseBuffers::delta_mean` or
 * `ScoreMarkersPairwiseBuffers::delta_detected`
 * may be NULL, in which case the corresponding statistic is not computed.
 */
template<typename Value_, typename Index_, typename Group_, typename Stat_>
void score_markers_pairwise(
    const tatami::Matrix<Value_, Index_>& matrix, 
    const Group_* group, 
    const ScoreMarkersPairwiseOptions& options,
    const ScoreMarkersPairwiseBuffers<Stat_>& output) 
{
    Index_ NC = matrix.ncol();
    auto group_sizes = tatami_stats::tabulate_groups(group, NC); 

    // In most cases this doesn't really matter, but we do it for consistency with the 1-block case,
    // and to account for variable weighting where non-zero block sizes get zero weight.
    auto group_weights = scran_blocks::compute_weights<Stat_>(group_sizes, options.block_weight_policy, options.variable_block_weight_parameters);

    size_t ngroups = group_sizes.size();
    size_t payload_size = static_cast<size_t>(matrix.nrow()) * ngroups; // cast to size_t to avoid overflow.
    std::vector<Stat_> group_means(payload_size), group_vars(payload_size), group_detected(payload_size);

    if (output.auc != NULL || matrix.prefer_rows()) {
        internal::scan_matrix_by_row<true>(
            matrix, 
            ngroups,
            group,
            1,
            static_cast<int*>(NULL),
            ngroups,
            NULL,
            group_means,
            group_vars,
            group_detected,
            output.auc,
            group_sizes,
            group_weights,
            options.threshold,
            options.num_threads
        );

    } else {
        internal::scan_matrix_by_column(
            matrix,
            ngroups,
            group,
            group_means,
            group_vars,
            group_detected,
            group_sizes,
            options.num_threads
        );
    }

    internal::process_simple_pairwise_effects(
        matrix.nrow(),
        ngroups,
        1,
        ngroups,
        group_means,
        group_vars,
        group_detected,
        output,
        group_weights,
        options.threshold,
        options.num_threads);
}

/**
 * Compute effect sizes for pairwise comparisons between groups, accounting for any blocking factor in the dataset.
 * Comparisons are only performed between the groups of cells in the same level of the blocking factor.
 * The batch-specific effect sizes are then combined into a single aggregate value for output.
 * This strategy avoids most problems related to batch effects as we never directly compare across different blocking levels.
 *
 * Specifically, for each gene and each pair of groups, we obtain one effect size per blocking level.
 * We consolidate these into a single statistic by computing the weighted mean across levels.
 * The weight for each level is defined as the product of the weights of the two groups involved in the comparison,
 * where each weight is computed from the size of the group using the logic described in `scran_blocks::compute_weights()`.
 *
 * Obviously, blocking levels with no cells in either group will not contribute anything to the weighted mean.
 * If two groups never co-occur in the same blocking level, no effect size will be computed and a `NaN` is reported in the output.
 * We do not attempt to reconcile batch effects in a partially confounded scenario.
 *
 * For the mean and detected proportion in each group, we compute a weighted average of each statistic across blocks for each gene.
 * Again, the weight for each block is defined from `scran_blocks::compute_weights()` on the size of the group in that block.
 *
 * @tparam Value_ Matrix data type.
 * @tparam Index_ Matrix index type.
 * @tparam Group_ Integer type for the group assignments.
 * @tparam Block_ Integer type for the block assignments.
 * @tparam Stat_ Floating-point type to store the statistics.
 *
 * @param matrix A **tatami** matrix instance.
 * @param[in] group Pointer to an array of length equal to the number of columns in `matrix`, containing the group assignments.
 * Group identifiers should be 0-based and should contain all integers in \f$[0, N)\f$ where \f$N\f$ is the number of unique groups.
 * @param[in] block Pointer to an array of length equal to the number of columns in `matrix`, containing the blocking factor.
 * Block identifiers should be 0-based and should contain all integers in \f$[0, B)\f$ where \f$B\f$ is the number of unique blocking levels.
 * @param options Further options.
 * @param[out] output Collection of buffers in which to store the computed statistics.
 * Each buffer is filled with the corresponding statistic for each group or pairwise comparison.
 * Any of `ScoreMarkersPairwiseBuffers::cohens_d`, 
 * `ScoreMarkersPairwiseBuffers::auc`, 
 * `ScoreMarkersPairwiseBuffers::delta_mean` or
 * `ScoreMarkersPairwiseBuffers::delta_detected`
 * may be NULL, in which case the corresponding statistic is not computed.
 */
template<typename Value_, typename Index_, typename Group_, typename Block_, typename Stat_>
void score_markers_pairwise_blocked(
    const tatami::Matrix<Value_, Index_>& matrix, 
    const Group_* group, 
    const Block_* block,
    const ScoreMarkersPairwiseOptions& options,
    const ScoreMarkersPairwiseBuffers<Stat_>& output) 
{
    Index_ NC = matrix.ncol();
    size_t ngroups = output.mean.size();
    size_t nblocks = tatami_stats::total_groups(block, NC); 

    auto combinations = internal::create_combinations(ngroups, group, block, NC);
    auto combo_sizes = internal::tabulate_combinations<Index_>(ngroups, nblocks, combinations);
    size_t ncombos = combo_sizes.size();
    auto combo_weights = scran_blocks::compute_weights<Stat_>(combo_sizes, options.block_weight_policy, options.variable_block_weight_parameters);

    size_t payload_size = static_cast<size_t>(matrix.nrow()) * ncombos; // cast to size_t to avoid overflow.
    std::vector<Stat_> combo_means(payload_size), combo_vars(payload_size), combo_detected(payload_size);

    if (output.auc != NULL || matrix.prefer_rows()) {
        internal::scan_matrix_by_row<false>(
            matrix, 
            ngroups,
            group,
            nblocks,
            block,
            ncombos,
            combinations.data(),
            combo_means,
            combo_vars,
            combo_detected,
            output.auc,
            combo_sizes,
            combo_weights, 
            options.threshold,
            options.num_threads
        );

    } else {
        internal::scan_matrix_by_column(
            matrix,
            ncombos,
            combinations.data(),
            combo_means,
            combo_vars,
            combo_detected,
            combo_sizes,
            options.num_threads
        );
    }
 
    internal::process_simple_pairwise_effects(
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
        options.num_threads);
}

/**
 * Overload of `score_markers_pairwise()` that allocates memory for the output statistics.
 *
 * @tparam Stat_ Floating-point type to store the statistics.
 * @tparam Value_ Matrix data type.
 * @tparam Index_ Matrix index type.
 * @tparam Group_ Integer type for the group assignments.
 *
 * @param matrix A **tatami** matrix instance.
 * @param[in] group Pointer to an array of length equal to the number of columns in `matrix`, containing the group assignments.
 * Group identifiers should be 0-based and should contain all integers in \f$[0, N)\f$ where \f$N\f$ is the number of unique groups.
 * @param options Further options.
 *
 * @return Object containing the pairwise effects, plus the mean expression and detected proportion in each group.
 */
template<typename Stat_ = double, typename Value_, typename Index_, typename Group_>
ScoreMarkersPairwiseResults<Stat_> score_markers_pairwise(const tatami::Matrix<Value_, Index_>& matrix, const Group_* group, const ScoreMarkersPairwiseOptions& options) {
    size_t ngroups = tatami_stats::total_groups(group, matrix.ncol());
    ScoreMarkersPairwiseResults<Stat_> res;
    auto buffers = internal::fill_pairwise_results(matrix.nrow(), ngroups, res, options);
    score_markers_pairwise(matrix, group, options, buffers);
    return res; 
}

/**
 * Overload of `score_markers_pairwise_blocked()` that allocates memory for the output statistics.
 *
 * @tparam Stat_ Floating-point type to store the statistics.
 * @tparam Value_ Matrix data type.
 * @tparam Index_ Matrix index type.
 * @tparam Group_ Integer type for the group assignments.
 * @tparam Block_ Integer type for the block assignments. 
 *
 * @param matrix A **tatami** matrix instance.
 * @param[in] group Pointer to an array of length equal to the number of columns in `matrix`, containing the group assignments.
 * Group identifiers should be 0-based and should contain all integers in \f$[0, N)\f$ where \f$N\f$ is the number of unique groups.
 * @param[in] block Pointer to an array of length equal to the number of columns in `matrix`, containing the blocking factor.
 * Block identifiers should be 0-based and should contain all integers in \f$[0, B)\f$ where \f$B\f$ is the number of unique blocking levels.
 * @param options Further options.
 *
 * @return Object containing the pairwise effects, plus the mean expression and detected proportion in each group.
 */
template<typename Stat_ = double, typename Value_, typename Index_, typename Group_, typename Block_>
ScoreMarkersPairwiseResults<Stat_> score_markers_pairwise_blocked(const tatami::Matrix<Value_, Index_>& matrix, const Group_* group, const Block_* block, const ScoreMarkersPairwiseOptions& options) {
    size_t ngroups = tatami_stats::total_groups(group, matrix.ncol());
    ScoreMarkersPairwiseResults<Stat_> res;
    auto buffers = internal::fill_pairwise_results(matrix.nrow(), ngroups, res, options);
    score_markers_pairwise_blocked(matrix, group, block, options, buffers);
    return res;
}

}

#endif
