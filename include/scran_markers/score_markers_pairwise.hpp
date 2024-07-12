#ifndef SCRAN_PAIRWISE_EFFECTS_HPP
#define SCRAN_PAIRWISE_EFFECTS_HPP

#include "MatrixCalculator.hpp"
#include "cohens_d.hpp"
#include "simple_diff.hpp"

#include "scran_blocks/scran_blocks.hpp"
#include "tatami/tatami.hpp"

/**
 * @file score_markers_pairwise.hpp
 * @brief Compute pairwise effect sizes between groups of cells.
 */

namespace scran {

/**
 * @brief Options for `score_markers_pairwise()`.
 */
struct ScoreMarkersPairwiseOptions {
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
     * Whether to compute Cohen's d. 
     * This only affects the `score_markers_pairwise()` overload that return `Results`.
     */
    bool compute_cohens_d = true;

    /**
     * Whether to compute the AUC.
     * This only affects the `score_markers_pairwise()` overload that return `Results`.
     */
    bool compute_auc = true;

    /**
     * Whether to compute the difference in means.
     * This only affects the `score_markers_pairwise()` overload that return `Results`.
     */
    bool compute_delta_mean = true;

    /**
     * Whether to compute the difference in the detected proportion.
     * This only affects the `score_markers_pairwise()` overload that return `Results`.
     */
    bool compute_delta_detected = true;

    /**
     * Policy to use for weighting blocks when computing average statistics/effect sizes across blocks.
     */
    scran_blocks::WeightPolicy block_weight_policy = WeightPolicy::VARIABLE;

    /**
     * Parameters for the variable block weights.
     * Only used when `Options::block_weight_policy = scran_blocks::WeightPolicy::VARIABLE`.
     */
    scran_blocks::VariableWeightParameters variable_block_weight_parameters;
};

/**
 * @brief Buffers for `score_markers_pairwise()`.
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
     * Pointer to an array of length equal to \f$NG^2\f$, where \f$N\f$ is the number of genes and \f$G\f$ is the number of groups.
     * This is a 3-dimensional array to be filled with the Cohen's D for the comparison between each pair of groups for each gene.
     *
     * The first dimension is the fastest changing, is of length equal to the number of groups, and represents the first group.
     * The second dimension is the next fastest changing, is also of length equal to the number of groups, and represents the second group.
     * The third dimension is the slowest changing, is of length equal to the number of genes, and represents the gene.
     * Thus, an entry \f$(i, j, k)\f$ represents the effect size of gene $k$ for group $i$ against group $j$.
     *
     * Alternatively NULL, in which case the Cohen's D is not stored.
     */
    Stat_* cohens_d = NULL;

    /**
     * Pointer to an array of length equal to \f$NG^2\f$, where \f$N\f$ is the number of genes and \f$G\f$ is the number of groups.
     * This is a 3-dimensional array to be filled with the AUC for the comparison between each pair of groups for each gene;
     * see `ScoreMarkersBuffers::cohens_d` for more details.
     * Alternatively NULL, in which case the AUC is not stored.
     */
    Stat_* auc = NULL;

    /**
     * Pointer to an array of length equal to \f$NG^2\f$, where \f$N\f$ is the number of genes and \f$G\f$ is the number of groups.
     * This is a 3-dimensional array to be filled with the difference in means for the comparison between each pair of groups for each gene;
     * see `ScoreMarkersBuffers::cohens_d` for more details.
     * Alternatively NULL, in which case the difference in means is not stored.
     */
    Stat_* delta_mean = NULL;

    /**
     * Pointer to an array of length equal to \f$NG^2\f$, where \f$N\f$ is the number of genes and \f$G\f$ is the number of groups.
     * This is a 3-dimensional array to be filled with the difference in the detected proportions for the comparison between each pair of groups for each gene;
     * see `ScoreMarkersBuffers::cohens_d` for more details.
     * Alternatively NULL, in which case the difference in detected proportions is not stored.
     */
    Stat_* delta_detected = NULL;
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
    double threshold,
    int num_threads)
{
    tatami::parallelize([&](size_t, Index_ start, Index_ length) -> void {
        size_t in_offset = ncombos * static_cast<size_t>(start);
        const auto* tmp_means = combo_means.data() + in_offset;
        const auto* tmp_variances = combo_vars.data() + in_offset;
        const auto* tmp_detected = combo_detected.data() + in_offset;

        size_t squared = ngroups * ngroups;
        size_t out_offset = start * squared;
        for (size_t gene = start, end = start + length; gene < end; ++gene, out_offset += squared) {
            for (size_t co = 0; co < ncombos; ++co) {
                output.means[co][gene] = tmp_means[co];
                output.detected[co][gene] = tmp_detected[co];
            }

            if (output.cohens_d != NULL) {
                differential_analysis::compute_pairwise_cohens_d(tmp_means, tmp_variances, state.level_weight, ngroups, nblocks, threshold, output.cohens_d + out_offset);
            }

            if (output.delta_detected != NULL) {
                differential_analysis::compute_pairwise_simple_diff(tmp_detected, state.level_weight, ngroups, nblocks, output.delta_detected + out_offset);
            }

            if (output.delta_mean != NULL) {
                differential_analysis::compute_pairwise_simple_diff(tmp_means, state.level_weight, ngroups, nblocks, output.delta_mean + out_offset);
            }

            tmp_means += nlevels;
            tmp_variances += nlevels;
            tmp_detected += nlevels;
        }
    }, ngenes, nthreads);
}

}
/**
 * @endcond
 */

/**
 * @brief Compute pairwise effect size between groups of cells.
 *
 * This class computes the effect sizes for the pairwise comparisons used in `ScoreMarkers`, prior to any ranking of marker genes.
 * It may be desirable to call this function directly if the pairwise effects themselves are of interest, rather than per-group summaries.
 *
 * @section effect-sizes Choice of effect sizes
 * The log-fold change (LFC) is the difference in the mean log-expression between groups.
 * This is fairly straightforward to interpret - as log-fold change of +1 corresponds to a two-fold upregulation in the first group compared to the second.
 * For this interpretation, we assume that the input matrix contains log-transformed normalized expression values.
 *
 * The delta-detected is the difference in the proportion of cells with detected expression between groups.
 * This lies between 1 and -1, with the extremes occurring when a gene is silent in one group and detected in all cells of the other group.
 * For this interpretation, we assume that the input matrix contains non-negative expression values, where a value of zero corresponds to lack of detectable expression.
 *
 * Cohen's d is the standardized log-fold change between two groups.
 * This is defined as the difference in the mean log-expression for each group scaled by the average standard deviation across the two groups.
 * (Technically, we should use the pooled variance; however, this introduces some unpleasant asymmetry depending on the variance of the larger group, so we take a simple average instead.)
 * A positive value indicates that the gene is upregulated in the first gene compared to the second.
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
 * @section lfc-threshold With a log-fold change threshold
 * Setting a log-fold change threshold can be helpful as it prioritizes genes with large shifts in expression instead of those with low variances.
 * Currently, only positive thresholds are supported - this focuses on genes upregulated in the first group compared to the second.
 * The effect size definitions are generalized when testing against a non-zero log-fold change threshold.
 *
 * Cohen's d is redefined as the standardized difference between the observed log-fold change and the specified threshold, analogous to the TREAT method from **limma**.
 * Large positive values are only obtained when the observed log-fold change is significantly greater than the threshold.
 * For example, if we had a threshold of 2 and we obtained a Cohen's d of 3, this means that the observed log-fold change was 3 standard deviations above 2.
 * Importantly, a negative Cohen's d cannot be intepreted as downregulation, as the log-fold change may still be positive but less than the threshold.
 * 
 * The AUC generalized to the probability of obtaining a random observation in one group that is greater than a random observation plus the threshold in the other group.
 * For example, if we had a threshold of 2 and we obtained an AUC of 0.8, this means that - 80% of the time - 
 * the random observation from the first group would be greater than a random observation from the second group by 2 or more.
 * Again, AUCs below 0.5 cannot be interpreted as downregulation, as it may be caused by a positive log-fold change that is less than the threshold.
 * 
 * @section blocked Blocked comparisons
 * In the presence of multiple batches, we can block on the batch of origin for each cell.
 * Comparisons are only performed between the groups of cells in the same batch (also called "blocking level" below).
 * The batch-specific effect sizes are then combined into a single aggregate value for output.
 * This strategy avoids most problems related to batch effects as we never directly compare across different blocking levels.
 *
 * Specifically, for each gene and each pair of groups, we obtain one effect size per blocking level.
 * We consolidate these into a single statistic by computing the weighted mean across levels.
 * The weight for each level is defined as the product of the weights of the two groups involved in the comparison,
 * where each weight is computed from the size of the group using the logic described in `variable_block_weight()`.
 *
 * Obviously, blocking levels with no cells in either group will not contribute anything to the weighted mean.
 * If two groups never co-occur in the same blocking level, no effect size will be computed and a `NaN` is reported in the output.
 * We do not attempt to reconcile batch effects in a partially confounded scenario.
 *
 * @section other Other statistics
 * We report the mean log-expression of all cells in each group, as well as the proportion of cells with detectable expression in each group.
 * These statistics are useful for quickly interpreting the differences in expression driving the effect size summaries.
 * If blocking is involved, we compute the grand average across blocks of the mean and proportion for each group,
 * where the weight for each block is defined from `variable_block_weight()` on the size of the group in that block.
 *
 * @tparam Value_ Matrix data type.
 * @tparam Index_ Matrix index type.
 * @tparam Group_ Integer type for the group assignments.
 * @tparam Stat_ Floating-point type to store the statistics.
 *
 * @param p Pointer to a **tatami** matrix instance.
 * @param[in] group Pointer to an array of length equal to the number of columns in `p`, containing the group assignments.
 * Group identifiers should be 0-based and should contain all integers in $[0, N)$ where $N$ is the number of unique groups.
 * @param[out] means Vector of length equal to the number of groups.
 * Each element corresponds to a group and is a pointer to an array of length equal to the number of rows in `p`.
 * This is used to store the mean expression of each group across all genes.
 * @param[out] detected Vector of length equal to the number of groups,
 * Each element corresponds to a group and is a pointer to an array of length equal to the number of rows in `p`.
 * This is used to store the proportion of detected expression in each group.
 * @param[out] cohen Pointer to an array of length equal to $GN^2$, where `N` is the number of groups and `G` is the number of genes (see `Results` for details).
 * This is filled with the Cohen's d for the pairwise comparisons between groups across all genes.
 * Ignored if set to `nullptr`, in which case Cohen's d is not computed.
 * @param[out] auc Pointer to an array as described for `cohen`, but instead storing the AUC.
 * Ignored if set to `nullptr`, in which case the AUC is not computed.
 * @param[out] lfc Pointer to an array as described for `cohen`, but instead storing the log-fold change. 
 * Ignored if set to `nullptr`, in which case the log-fold change is not computed.
 * @param[out] delta_detected Pointer to an array as described for `cohen`, but instead the delta in the detected proportions.
 * Ignored if set to `nullptr`, in which case the delta detected is not computed.
 */
template<typename Value_, typename Index_, typename Group_, typename Stat_>
void score_markers_pairwise(
    const tatami::Matrix<Value_, Index_>& matrix, 
    const Group_* group, 
    const ScoreMarkersPairwiseOptions& options,
    const ScoreMarkersPairwiseBuffers<Stat_>& output) 
{
    Index_ NC = matrix->ncol();
    auto group_sizes = tatami_stats::tabulate_groups(NC, group); 
    auto group_weights = scran_blocks::compute_weights(group_sizes, options.block_weight_policy, options.variable_block_weight_parameters);

    size_t payload_size = static_cast<size_t>(matrix->nrow()) * group_sizes.size(); 
    std::vector<Stat_> group_means(payload_size), group_vars(payload_size), group_detected(payload_size);

    if (output.auc != NULL || matrix.prefer_rows());
        internal::scan_matrix_by_row<true>(
            matrix, 
            ngroups,
            group,
            1,
            static_cast<Block_*>(NULL),
            ngroups,
            NULL,
            group_means,
            group_vars,
            group_detected,
            output.auc,
            group_sizes,
            group_weights
            options.num_threads
        );

    } else {
        internal::scan_matrix_by_column(
            matrix,
            ngroups,
            group,
            group_means,
            group_vars,
            group_detected
        );
    }
}

template<typename Value_, typename Index_, typename Group_, typename Block_, typename Stat_>
void score_markers_pairwise_blocked(
    const tatami::Matrix<Value_, Index_>& matrix, 
    const Group_* group, 
    const Block_* block,
    const ScoreMarkersPairwiseOptions& options,
    const ScoreMarkersPairwiseBuffers<Stat_>& output) 
{
    Index_ NC = matrix->ncol();
    size_t ngroups = output.mean.size();
    size_t nblocks = tatami_stats::count_groups(NC, block); 

    std::vector<size_t> combined(NC);
    for (Index_ c = 0; c < NC; ++c) {
        combined[c] = static_cast<size_t>(group[c]) + static_cast<size_t>(block[c]) * ngroups;
    }

    if (output.auc != NULL || matrix.prefer_rows());
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
            options.num_threads
        );

    } else {
        if (nblocks == 1) {
        }

    }
 


}


        std::vector<Stat_*> means, 
        std::vector<Stat_*> detected, 
        Stat_* cohen,
        Stat_* auc,
        Stat_* lfc,
        Stat_* delta_detected) 
    const {
        int ngroups = means.size();
        differential_analysis::MatrixCalculator runner(nthreads, threshold, block_weight_policy, variable_block_weight_parameters);
        Overlord overlord(auc);
        auto state = runner.run(p, group, ngroups, overlord);
        process_simple_effects(p->nrow(), ngroups, 1, state, means, detected, cohen, lfc, delta_detected);
    }

    /**
     * Compute effect sizes for pairwise comparisons between groups, accounting for any blocking factor in the dataset.
     * On completion, `means`, `detected`, `cohen`, `auc`, `lfc` and `delta_detected` are filled with their corresponding statistics. 
     *
     * @tparam Data_ Matrix data type.
     * @tparam Index_ Matrix index type.
     * @tparam Group_ Integer type for the group assignments.
     * @tparam Block_ Integer type for the block assignments.
     * @tparam Stat_ Floating-point type to store the statistics.
     *
     * @param p Pointer to a **tatami** matrix instance.
     * @param[in] group Pointer to an array of length equal to the number of columns in `p`, containing the group assignments.
     * Group identifiers should be 0-based and should contain all integers in $[0, N)$ where $N$ is the number of unique groups.
     * @param[in] block Pointer to an array of length equal to the number of columns in `p`, containing the blocking factor.
     * Block identifiers should be 0-based and should contain all integers in $[0, N)$ where $N$ is the number of unique groups.
     * This can also be `NULL` in which case all cells are assumed to belong to the same block (i.e., like `run()`).
     * @param[out] means Vector of length equal to the number of groups.
     * Each element corresponds to a group and is another vector of length equal to the number of genes.
     * Each entry of the inner vector contains the grand average of the mean expression across all blocks.
     * @param[out] detected Vector of length equal to the number of groups.
     * Each element corresponds to a group and is another vector of length equal to the number of genes.
     * Each entry of the inner vector contains the grand average of the detected proportion across all blocks.
     * @param[out] cohen Pointer to an array of length equal to $GN^2$, where `N` is the number of groups and `G` is the number of genes (see `Results` for details).
     * This is filled with the Cohen's d for the pairwise comparisons between groups across all genes.
     * Ignored if set to `nullptr`, in which case Cohen's d is not computed.
     * @param[out] auc Pointer to an array as described for `cohen`, but instead storing the AUC.
     * Ignored if set to `nullptr`, in which case the AUC is not computed.
     * @param[out] lfc Pointer to an array as described for `cohen`, but instead storing the log-fold change. 
     * Ignored if set to `nullptr`, in which case the log-fold change is not computed.
     * @param[out] delta_detected Pointer to an array as described for `cohen`, but instead the delta in the detected proportions.
     * Ignored if set to `nullptr`, in which case the delta detected is not computed.
     */
    template<typename Data_, typename Index_, typename Group_, typename Block_, typename Stat_>
    void run_blocked(
        const tatami::Matrix<Data_, Index_>* p, 
        const Group_* group, 
        const Block_* block, 
        std::vector<Stat_*> means, 
        std::vector<Stat_*> detected, 
        Stat_* cohen,
        Stat_* auc,
        Stat_* lfc,
        Stat_* delta_detected)
    const {
        if (block == NULL) {
            run(p, group, std::move(means), std::move(detected), cohen, auc, lfc, delta_detected);

        } else {
            auto ngenes = p->nrow();
            int ngroups = means.size();
            int nblocks = count_ids(p->ncol(), block);

            int ncombos = ngroups * nblocks;
            std::vector<std::vector<Stat_> > means_store(ncombos), detected_store(ncombos);
            std::vector<Stat_*> means2(ncombos), detected2(ncombos);
            for (int c = 0; c < ncombos; ++c) {
                means_store[c].resize(ngenes);
                detected_store[c].resize(ngenes);
                means2[c] = means_store[c].data();
                detected2[c] = detected_store[c].data();
            }

            differential_analysis::MatrixCalculator runner(nthreads, threshold, block_weight_policy, variable_block_weight_parameters);
            Overlord overlord(auc);
            auto state = runner.run_blocked(p, group, ngroups, block, nblocks, overlord);
            process_simple_effects(p->nrow(), ngroups, nblocks, state, means2, detected2, cohen, lfc, delta_detected);

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
    }

private:
    template<typename Stat_>
    struct Overlord {
        Overlord(Stat_* auc_) : auc(auc_) {}

        bool needs_auc() const { 
            return auc != NULL;
        }

        Stat_* auc;

        Stat_* prepare_auc_buffer(size_t gene, size_t ngroups) { 
            return auc + gene * ngroups * ngroups;
        }
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
     * Group identifiers should be 0-based and should contain all integers in $[0, N)$ where $N$ is the number of unique groups.
     *
     * @return A `Results` object is returned containing the pairwise effects, plus the mean expression and detected proportion in each group.
     */
    template<typename Stat_ = double, typename Data_ = double, typename Index_ = int, typename Group_ = int>
    Results<Stat_> run(const tatami::Matrix<Data_, Index_>* p, const Group_* group) {
        size_t ngroups = count_ids(p->ncol(), group);
        Results<Stat_> res(p->nrow(), ngroups, do_cohen, do_auc, do_lfc, do_delta_detected); 
        run(
            p, 
            group, 
            vector_to_pointers(res.means), 
            vector_to_pointers(res.detected), 
            harvest_pointer(res.cohen, do_cohen),
            harvest_pointer(res.auc, do_auc),
            harvest_pointer(res.lfc, do_lfc),
            harvest_pointer(res.delta_detected, do_delta_detected)
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
     * Group identifiers should be 0-based and should contain all integers in $[0, N)$ where $N$ is the number of unique groups.
     * @param[in] block Pointer to an array of length equal to the number of columns in `p`, containing the blocking factor.
     * See `run_blocked()` for more details.
     *
     * @return A `Results` object is returned containing the pairwise effects, plus the mean expression and detected proportion in each group and block.
     */
    template<typename Stat_ = double, typename Data_ = double, typename Index_ = int, typename Group_ = int, typename Block_ = int>
    Results<Stat_> run_blocked(const tatami::Matrix<Data_, Index_>* p, const Group_* group, const Block_* block) {
        size_t ngroups = count_ids(p->ncol(), group);
        Results<Stat_> res(p->nrow(), ngroups, do_cohen, do_auc, do_lfc, do_delta_detected); 
        run_blocked(
            p,
            group,
            block,
            vector_to_pointers(res.means),
            vector_to_pointers(res.detected),
            harvest_pointer(res.cohen, do_cohen),
            harvest_pointer(res.auc, do_auc),
            harvest_pointer(res.lfc, do_lfc),
            harvest_pointer(res.delta_detected, do_delta_detected)
        );
        return res;
    }

private:
    template<typename Ptr>
    static std::vector<Ptr> fetch_first(const std::vector<std::vector<Ptr> >& input) {
        std::vector<Ptr> output;
        output.reserve(input.size());
        for (const auto& i : input) {
            output.push_back(i.front());
        }
        return output;
    }

    template<typename Stat> 
    static Stat* harvest_pointer(std::vector<Stat>& source, bool use) {
        if (use) {
            return source.data();
        } else {
            return static_cast<Stat*>(NULL);
        }
    }
};

}

#endif
