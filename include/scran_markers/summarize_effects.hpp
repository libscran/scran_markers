#ifndef SCRAN_MARKERS_SUMMARIZE_EFFECTS_HPP
#define SCRAN_MARKERS_SUMMARIZE_EFFECTS_HPP

#include "summarize_comparisons.hpp"

#include <vector>

/**
 * @file summarize_effects.hpp
 * @brief Summarize the effect sizes from pairwise comparisons.
 */

namespace scran_markers {

/**
 * @brief Options for `summarize_effects()`.
 */
struct SummarizeEffectsOptions {
    /**
     * Number of threads to use. 
     */
    int num_threads = 1;

    /**
     * Whether to report the minimum of the pairwise effects.
     * Only affects the `summarize_effects()` overload that returns a vector of `SummaryResults`.
     */
    bool compute_min = true;

    /**
     * Whether to report the mean of the pairwise effects.
     * Only affects the `summarize_effects()` overload that returns a vector of `SummaryResults`.
     */
    bool compute_mean = true;

    /**
     * Whether to report the median of the pairwise effects.
     * Only affects the `summarize_effects()` overload that returns a vector of `SummaryResults`.
     */
    bool compute_median = true;

    /**
     * Whether to report the maximum of the pairwise effects.
     * Only affects the `summarize_effects()` overload that returns a vector of `SummaryResults`.
     */
    bool compute_max = true;

    /**
     * Whether to report the min-rank of the pairwise effects.
     * Only affects the `summarize_effects()` overload that returns a vector of `SummaryResults`.
     */
    bool compute_min_rank = true;
};

/**
 * Given `n` groups, each group is involved in `n - 1` pairwise comparisons and thus has `n - 1` effect sizes (e.g., as computed by `compute_pairwise_effects()`).
 * For each group, we compute summary statistics - e.g., the minimum, median, mean - of the effect sizes across all of that group's comparisons.
 * Users can then sort by any of these summaries to obtain a ranking of potential marker genes for each group.
 * 
 * The choice of summary statistic dictates the interpretation of the ranking.
 * Given a group \f$X\f$:
 * 
 * - A large mean effect size indicates that the gene is upregulated in \f$X\f$ compared to the average of the other groups.
 *   A small value indicates that the gene is downregulated in \f$X\f$ instead.
 *   This is a good general-purpose summary statistic for ranking, usually by decreasing size to obtain upregulated markers in \f$X\f$.
 * - A large median effect size indicates that the gene is upregulated in \f$X\f$ compared to most (>50%) other groups.
 *   A small value indicates that the gene is downregulated in \f$X\f$ instead.
 *   This is also a good general-purpose summary, with the advantage of being more robust to outlier effects compared to the mean.
 *   However, it also has the disadvantage of being less sensitive to strong effects in a minority of comparisons.
 * - A large minimum effect size indicates that the gene is upregulated in \f$X\f$ compared to all other groups.
 *   A small value indicates that the gene is downregulated in \f$X\f$ compared to at least one other group.
 *   For upregulation, this is the most stringent summary as markers will only have extreme values if they are _uniquely_ upregulated in \f$X\f$ compared to every other group.
 *   However, it may not be effective if \f$X\f$ is closely related to any of the groups.
 * - A large maximum effect size indicates that the gene is upregulated in \f$X\f$ compared to at least one other group.
 *   A small value indicates that the gene is downregulated in \f$X\f$ compared to all other groups.
 *   For downregulation, this is the most stringent summary as markers will only have extreme values if they are _uniquely_ downregulated in \f$X\f$ compared to every other group.
 *   However, it may not be effective if \f$X\f$ is closely related to any of the groups.
 * - The "minimum rank" (a.k.a. min-rank) is defined by ranking genes based on decreasing effect size _within_ each comparison, and then taking the smallest rank _across_ comparisons.
 *   A minimum rank of 1 means that the gene is the top upregulated gene in at least one comparison to another group.
 *   More generally, a minimum rank of \f$T\f$ indicates that the gene is the \f$T\f$-th upregulated gene in at least one comparison. 
 *   Applying a threshold on the minimum rank is useful for obtaining a set of genes that, in combination, are guaranteed to distinguish \f$X\f$ from every other group.
 *
 * The exact definition of "large" and "small" depends on the choice of effect size. 
 * For signed effects like Cohen's d, log-fold change and delta-detected, the value must be positive to be considered "large", and negative to be considered "small".
 * For the AUC, a value greater than 0.5 is considered "large" and less than 0.5 is considered "small".
 *
 * The interpretation above is also contingent on the log-fold change threshold used (see `compute_pairwise_effects()` or `score_markers()` for details).
 * For positive thresholds, small effects cannot be unambiguously interpreted as downregulation, as the effect is already adjusted to account for the threshold.
 * As a result, only large effects can be interpreted as evidence for upregulation.
 *
 * NaN effect sizes are allowed, e.g., if two groups do not exist in the same block for a blocked analysis in `compute_pairwise_effects_blocked()`.
 * This class will ignore NaN values when computing each summary.
 * If all effects are NaN for a particular group, the summary statistic will also be `NaN`.
 *
 * All choices of summary statistics are enumerated by `Summary`.
 *
 * @tparam Index_ Integer type for the number of genes.
 * @tparam Stat_ Floating-point type for the statistics.
 * @tparam Rank_ Numeric type for the minimum rank.
 *
 * @param ngenes Number of genes.
 * @param ngroups Number of groups.
 * @param[in] effects Pointer to a 3-dimensional array containing the pairwise statistics.
 * The first dimension is the fastest changing, is of length equal to the number of groups, and represents the first group.
 * The second dimension is the next fastest changing, is also of length equal to the number of groups, and represents the second group.
 * The third dimension is the slowest changing, is of length equal to the number of genes, and represents the gene.
 * Thus, an entry \f$(i, j, k)\f$ represents the effect size of gene $k$ for group $i$ against group $j$.
 * @param[out] summaries Vector of length equal to the number of groups.
 * Each entry corresponds to a group and is used to store the summary statistics for that group.
 * Each pointer in any given `SummaryBuffers` should either point to an array of length equal to the number of genes, 
 * or be NULL to indicate that the corresponding summary statistic should not be computed for that group.
 * @param options Further options.
 */
template<typename Index_, typename Stat_, typename Rank_>
void summarize_effects(Index_ ngenes, size_t ngroups, const Stat_* effects, const std::vector<SummaryBuffers<Stat_, Rank_> >& summaries, const SummarizeEffectsOptions& options) {
    internal::compute_min_rank_pairwise(ngenes, ngroups, effects, summaries, options.num_threads);
    internal::summarize_comparisons(ngenes, ngroups, effects, summaries, options.num_threads); 
}

/**
 * Overload of `summarize_effects()` that allocates memory for the output summary statistics.
 *
 * @tparam Index_ Integer type for the number of genes.
 * @tparam Stat Floating point type for the statistics.
 * @tparam Rank_ Numeric type for the minimum rank.
 *
 * @param ngenes Number of genes.
 * @param ngroups Number of groups.
 * @param[in] effects Pointer to a 3-dimensional array containing the pairwise statistics.
 * The first dimension is the fastest changing, is of length equal to the number of groups, and represents the first group.
 * The second dimension is the next fastest changing, is also of length equal to the number of groups, and represents the second group.
 * The third dimension is the slowest changing, is of length equal to the number of genes, and represents the gene.
 * Thus, an entry \f$(i, j, k)\f$ represents the effect size of gene $k$ for group $i$ against group $j$.
 * @param options Further options.
 *
 * @return A vector of length equal to the number of groups.
 * Each `SummaryResults` corresponds to a group and contains the summary statistics (depending on `options`) for that group.
 */
template<typename Stat_ = double, typename Rank_ = int, typename Index_>
std::vector<SummaryResults<Stat_, Rank_> > summarize_effects(Index_ ngenes, size_t ngroups, const Stat_* effects, const SummarizeEffectsOptions& options) {
    std::vector<SummaryResults<Stat_, Rank_> > output(ngroups);
    std::vector<SummaryBuffers<Stat_, Rank_> > ptrs(ngroups);

    for (size_t g = 0; g < ngroups; ++g) {
        auto& out = output[g];
        auto& ptr = ptrs[g];

        if (options.compute_min) {
            out.min.resize(ngenes);
            ptr.min = out.min.data();
        }
        if (options.compute_mean) {
            out.mean.resize(ngenes);
            ptr.mean = out.mean.data();
        }
        if (options.compute_median) {
            out.median.resize(ngenes);
            ptr.median = out.median.data();
        }
        if (options.compute_max) {
            out.max.resize(ngenes);
            ptr.max = out.max.data();
        }
        if (options.compute_min_rank) {
            out.min_rank.resize(ngenes);
            ptr.min_rank = out.min_rank.data();
        }
    }

    summarize_effects(ngenes, ngroups, effects, ptrs, options);
    return output;
}

}

#endif
