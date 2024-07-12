#include <gtest/gtest.h>

#include "scran_markers/cohens_d.hpp"

TEST(CohensD, RawEdgeCases) {
    EXPECT_EQ(scran_markers::internal::compute_cohens_d<double>(1, 0, 0.5, 0), 2.0);
    EXPECT_EQ(scran_markers::internal::compute_cohens_d<double>(1, 0, 0.5, 1), 0.0);
    EXPECT_EQ(scran_markers::internal::compute_cohens_d<double>(1, 0, 0, 0), std::numeric_limits<double>::infinity());
    EXPECT_EQ(scran_markers::internal::compute_cohens_d<double>(0, 1, 0, 0), -std::numeric_limits<double>::infinity());
    EXPECT_EQ(scran_markers::internal::compute_cohens_d<double>(0, 0, 0, 0), 0);
    EXPECT_TRUE(std::isnan(scran_markers::internal::compute_cohens_d<double>(0, 0, std::numeric_limits<double>::quiet_NaN(), 0)));
}

TEST(CohensD, Denominator) {
    EXPECT_FLOAT_EQ(scran_markers::internal::cohen_denominator<double>(4.0, 0), std::sqrt(2.0));
    EXPECT_FLOAT_EQ(scran_markers::internal::cohen_denominator<double>(5.0, 3.0), 2.0);

    double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FLOAT_EQ(scran_markers::internal::cohen_denominator<double>(4.0, nan), 2.0);
    EXPECT_FLOAT_EQ(scran_markers::internal::cohen_denominator<double>(nan, 4.0), 2.0);
    EXPECT_TRUE(std::isnan(scran_markers::internal::cohen_denominator<double>(nan, nan)));
}

TEST(CohensD, Unblocked) {
    std::vector<double> means{0.1, 0.2, 0.3, 0.4};
    std::vector<double> variances{1.5, 2.3, 0.5, 1.2};
    std::vector<int> group_sizes{ 10, 5, 12, 34 }; // don't really matter.

    size_t ngroups = means.size();
    std::vector<double> output(ngroups * ngroups);
    scran_markers::internal::PrecomputedPairwiseWeights preweights(group_sizes.size(), 1, group_sizes.data());
    scran_markers::internal::compute_pairwise_cohens_d(means.data(), variances.data(), ngroups, 1, preweights, 0.0, output.data());

    for (size_t g = 0; g < ngroups; ++g) {
        for (size_t g2 = 0; g2 < ngroups; ++g2) {
            EXPECT_FLOAT_EQ(output[g * means.size() + g2], (means[g] - means[g2]) / std::sqrt((variances[g] + variances[g2])/2.0));
        }
    }
}

TEST(CohensD, Thresholded) {
    std::vector<double> means{0.1, 0.2, 0.3, 0.4};
    std::vector<double> variances{1.5, 2.3, 0.5, 1.2};
    std::vector<int> group_sizes{ 10, 5, 12, 34 }; // don't really matter.

    size_t ngroups = means.size();
    std::vector<double> output(ngroups * ngroups);
    scran_markers::internal::PrecomputedPairwiseWeights preweights(group_sizes.size(), 1, group_sizes.data());
    scran_markers::internal::compute_pairwise_cohens_d(means.data(), variances.data(), ngroups, 1, preweights, 1.0, output.data());

    for (size_t g = 0; g < ngroups; ++g) {
        for (size_t g2 = 0; g2 < ngroups; ++g2) {
            if (g != g2) {
                EXPECT_FLOAT_EQ(output[g * ngroups + g2], (means[g] - means[g2] - 1.0) / std::sqrt((variances[g] + variances[g2])/2.0));
            }
        }
    }
}

TEST(CohensD, MissingValues) {
    double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> means{ nan, 0.2, 0.3, 0.4, 0.1 };
    std::vector<double> variances{ nan, 2.3, nan, nan, 1.2 };
    std::vector<int> group_sizes{ 0, 5, 1, 1, 5 }; 

    size_t ngroups = means.size();
    std::vector<double> output(means.size() * means.size());
    scran_markers::internal::PrecomputedPairwiseWeights preweights(ngroups, 1, group_sizes.data());
    scran_markers::internal::compute_pairwise_cohens_d(means.data(), variances.data(), ngroups, 1, preweights, 0.0, output.data());

    for (size_t g = 0; g < ngroups; ++g) {
        for (size_t g2 = 0; g2 < ngroups; ++g2) {
            if (g == g2) {
                continue;
            } 
            
            double x = output[g * ngroups + g2];
            if (std::isnan(means[g]) || std::isnan(means[g2]) || (std::isnan(variances[g]) && std::isnan(variances[g2]))) {
                EXPECT_TRUE(std::isnan(x));
            } else {
                double delta = means[g] - means[g2];
                if (!std::isnan(variances[g]) && !std::isnan(variances[g2])) {
                    EXPECT_FLOAT_EQ(x, delta / std::sqrt((variances[g] + variances[g2])/2.0));
                } else if (!std::isnan(variances[g])) {
                    EXPECT_FLOAT_EQ(x, delta / std::sqrt(variances[g]));
                } else {
                    EXPECT_FLOAT_EQ(x, delta / std::sqrt(variances[g2]));
                }
            }
        }
    }
}

TEST(CohensD, Blocked) {
    int nblocks = 2, ngroups = 4;
    std::vector<double> means{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
    std::vector<double> variances{1.5, 2.3, 0.5, 1.2, 0.1, 1.2, 0.4, 0.5 };
    std::vector<int> combo_sizes{ 10, 5, 12, 34, 15, 2, 3, 6 }; 

    std::vector<double> output(ngroups * ngroups);
    scran_markers::internal::PrecomputedPairwiseWeights preweights(ngroups, nblocks, combo_sizes.data());
    scran_markers::internal::compute_pairwise_cohens_d(means.data(), variances.data(), ngroups, nblocks, preweights, 0.0, output.data());

    for (int g1 = 0; g1 < ngroups; ++g1) {
        for (int g2 = 0; g2 < ngroups; ++g2) {
            double totalnum = 0, totaldenom = 0;

            for (int b = 0; b < nblocks; ++b) {
                int offset1 = g1 + b * ngroups;
                int offset2 = g2 + b * ngroups;
                double d = (means[offset1] - means[offset2])/std::sqrt((variances[offset1] + variances[offset2])/2);
                double w = combo_sizes[offset1] * combo_sizes[offset2];
                totalnum += d * w;
                totaldenom += w;
            }

            EXPECT_FLOAT_EQ(output[g1 * ngroups + g2], totalnum/totaldenom);
        }
    }
}

TEST(CohensD, BlockedMissing) {
    double nan = std::numeric_limits<double>::quiet_NaN();
    int nblocks = 2, ngroups = 4;
    std::vector<double> means{nan, 0.2, nan, nan, 0.5, 0.6, 0.7, 0.8 };
    std::vector<double> variances{nan, nan, nan, nan, 0.1, 1.2, 0.4, 0.5 };
    std::vector<int> combo_sizes{ 0, 1, 0, 0, 15, 2, 3, 6 }; 

    std::vector<double> output(ngroups * ngroups);
    scran_markers::internal::PrecomputedPairwiseWeights preweights(ngroups, nblocks, combo_sizes.data());
    scran_markers::internal::compute_pairwise_cohens_d(means.data(), variances.data(), ngroups, nblocks, preweights, 0.0, output.data());

    // Effectively excising the first block, as there is no comparison with a
    // non-zero product for the combined weight. We shouldn't see any NaNs
    // bleeding through, as the 0-size groups should be skipped.

    std::vector<double> output2(ngroups * ngroups);
    scran_markers::internal::PrecomputedPairwiseWeights preweights2(ngroups, nblocks - 1, combo_sizes.data() + ngroups);
    scran_markers::internal::compute_pairwise_cohens_d(means.data() + ngroups, variances.data() + ngroups, ngroups, nblocks - 1, preweights2, 0.0, output2.data());

    EXPECT_EQ(output, output2);
}
