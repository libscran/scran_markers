#include "scran_tests/scran_tests.hpp"

#include <cmath>
#include <vector>
#include "scran_markers/simple_diff.hpp"

TEST(SimpleDiff, Unblocked) {
    std::vector<double> means{0.1, 0.2, 0.3, 0.4};
    std::vector<double> combo_weights{ 10, 5, 12, 34 }; // don't really matter.

    const auto ngroups = means.size();
    std::vector<double> output(ngroups * ngroups
#ifdef SCRAN_MARKERS_TEST_INIT
        , SCRAN_MARKERS_TEST_INIT
#endif
    );
    {
        scran_markers::internal::PrecomputedPairwiseWeights preweights(combo_weights.size(), 1, combo_weights.data());
        scran_markers::internal::compute_pairwise_simple_diff_blockmean(means.data(), combo_weights.size(), 1, preweights, output.data());

        for (size_t g = 0; g < means.size(); ++g) {
            for (size_t g2 = 0; g2 < means.size(); ++g2) {
                EXPECT_FLOAT_EQ(output[g * means.size() + g2], means[g] - means[g2]);
            }
        }
    }

    {
        std::vector<double> q_output(ngroups * ngroups
#ifdef SCRAN_MARKERS_TEST_INIT
            , SCRAN_MARKERS_TEST_INIT
#endif
        );
        std::vector<double> buffer;
        scran_blocks::SingleQuantileVariable<double, typename std::vector<double>::iterator> qcalc(1, 0.5);
        scran_markers::internal::compute_pairwise_simple_diff_blockquantile(means.data(), ngroups, 1, buffer, qcalc, q_output.data());
        scran_tests::compare_almost_equal_containers(output, q_output, {});
    }
}

TEST(SimpleDiff, ZeroSize) {
    double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> means{ nan, 0.2, 0.3, 0.4, 0.1 };
    std::vector<double> combo_weights{ 0, 5, 1, 1, 5 }; 

    const auto ngroups = means.size();
    std::vector<double> output(ngroups * ngroups
#ifdef SCRAN_MARKERS_TEST_INIT
        , SCRAN_MARKERS_TEST_INIT
#endif
    );
    {
        scran_markers::internal::PrecomputedPairwiseWeights preweights(combo_weights.size(), 1, combo_weights.data());
        scran_markers::internal::compute_pairwise_simple_diff_blockmean(means.data(), combo_weights.size(), 1, preweights, output.data());

        for (size_t g = 0; g < means.size(); ++g) {
            for (size_t g2 = 0; g2 < means.size(); ++g2) {
                if (g == g2) {
                    continue;
                } 
                
                double x = output[g * means.size() + g2];
                if (combo_weights[g] == 0 || combo_weights[g2] == 0) {
                    EXPECT_TRUE(std::isnan(x));
                } else {
                    double delta = means[g] - means[g2];
                    EXPECT_FLOAT_EQ(x, delta);
                }
            }
        }
    }

    {
        std::vector<double> q_output(ngroups * ngroups
#ifdef SCRAN_MARKERS_TEST_INIT
            , SCRAN_MARKERS_TEST_INIT
#endif
        );
        std::vector<double> buffer;
        scran_blocks::SingleQuantileVariable<double, typename std::vector<double>::iterator> qcalc(1, 0.5);
        scran_markers::internal::compute_pairwise_simple_diff_blockquantile(means.data(), ngroups, 1, buffer, qcalc, q_output.data());
        scran_tests::compare_almost_equal_containers(output, q_output, {});
    }
}

TEST(SimpleDiff, Blocked) {
    int nblocks = 2, ngroups = 4;
    std::vector<double> means{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
    std::vector<double> combo_weights{ 10, 5, 12, 34, 15, 2, 3, 6 }; 

    std::vector<double> output(ngroups * ngroups
#ifdef SCRAN_MARKERS_TEST_INIT
        , SCRAN_MARKERS_TEST_INIT
#endif
    );
    {
        scran_markers::internal::PrecomputedPairwiseWeights preweights(ngroups, nblocks, combo_weights.data());
        scran_markers::internal::compute_pairwise_simple_diff_blockmean(means.data(), ngroups, nblocks, preweights, output.data());

        for (int g1 = 0; g1 < ngroups; ++g1) {
            for (int g2 = 0; g2 < ngroups; ++g2) {
                double totalnum = 0, totaldenom = 0;

                for (int b = 0; b < nblocks; ++b) {
                    int offset1 = g1 + b * ngroups;
                    int offset2 = g2 + b * ngroups;
                    double d = means[offset1] - means[offset2];
                    double w = combo_weights[offset1] * combo_weights[offset2];
                    totalnum += d * w;
                    totaldenom += w;
                }

                EXPECT_FLOAT_EQ(output[g1 * ngroups + g2], totalnum/totaldenom);
            }
        }
    }

    {
        std::vector<double> q_output(ngroups * ngroups
#ifdef SCRAN_MARKERS_TEST_INIT
            , SCRAN_MARKERS_TEST_INIT
#endif
        );
        std::vector<double> buffer;
        scran_blocks::SingleQuantileVariable<double, typename std::vector<double>::iterator> qcalc(nblocks, 0.5);
        scran_markers::internal::compute_pairwise_simple_diff_blockquantile(means.data(), ngroups, nblocks, buffer, qcalc, q_output.data());

        for (int g1 = 0; g1 < ngroups; ++g1) {
            for (int g2 = 0; g2 < ngroups; ++g2) {
                double totalnum = 0;

                for (int b = 0; b < nblocks; ++b) {
                    int offset1 = g1 + b * ngroups;
                    int offset2 = g2 + b * ngroups;
                    totalnum += means[offset1] - means[offset2];
                }

                // Default quantile of 0.5 is the median, and with just 2 blocks, this is the mean.
                EXPECT_FLOAT_EQ(output[g1 * ngroups + g2], totalnum/2);
            }
        }
    }
}

TEST(SimpleDiff, BlockedMissing) {
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    int nblocks = 2, ngroups = 4;
    std::vector<double> means{ 0.0, nan, nan, nan, 0.2, 0.4, 0.6, 0.8 };
    std::vector<double> combo_weights{ 0, 1, 0, 0, 34, 23, 6, 55 }; // exactly one group in the first block with non-zero size.

    std::vector<double> output(ngroups * ngroups
#ifdef SCRAN_MARKERS_TEST_INIT
        , SCRAN_MARKERS_TEST_INIT
#endif
    );
    {
        scran_markers::internal::PrecomputedPairwiseWeights preweights(ngroups, nblocks, combo_weights.data());
        scran_markers::internal::compute_pairwise_simple_diff_blockmean(means.data(), ngroups, nblocks, preweights, output.data());

        // Effectively excising the first block, as there is no comparison with a
        // non-zero product for the combined weight. We shouldn't see any NaNs
        // bleeding through, as the 0-size groups should be skipped.

        std::vector<double> output2(ngroups * ngroups);
        scran_markers::internal::PrecomputedPairwiseWeights preweights2(ngroups, nblocks - 1, combo_weights.data() + ngroups);
        scran_markers::internal::compute_pairwise_simple_diff_blockmean(means.data() + ngroups, ngroups, nblocks - 1, preweights2, output2.data());
        EXPECT_EQ(output, output2);
    }

    {
        std::vector<double> q_output(ngroups * ngroups
#ifdef SCRAN_MARKERS_TEST_INIT
            , SCRAN_MARKERS_TEST_INIT
#endif
        );
        std::vector<double> buffer;
        scran_blocks::SingleQuantileVariable<double, typename std::vector<double>::iterator> qcalc(nblocks, 0.5);
        scran_markers::internal::compute_pairwise_simple_diff_blockquantile(means.data(), ngroups, nblocks, buffer, qcalc, q_output.data());
        scran_tests::compare_almost_equal_containers(output, q_output, {});
    }
}
