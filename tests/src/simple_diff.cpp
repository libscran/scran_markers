#include <gtest/gtest.h>

#include <cmath>
#include <vector>
#include "scran_markers/simple_diff.hpp"

TEST(SimpleDiff, Unblocked) {
    std::vector<double> means{0.1, 0.2, 0.3, 0.4};
    std::vector<int> group_sizes{ 10, 5, 12, 34 }; // don't really matter.

    std::vector<double> output(means.size() * means.size());
    scran_markers::internal::compute_pairwise_simple_diff(means.data(), group_sizes.size(), 1, group_sizes.data(), output.data());

    for (size_t g = 0; g < means.size(); ++g) {
        for (size_t g2 = 0; g2 < means.size(); ++g2) {
            EXPECT_FLOAT_EQ(output[g * means.size() + g2], means[g] - means[g2]);
        }
    }
}

TEST(SimpleDiff, ZeroSize) {
    double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> means{ nan, 0.2, 0.3, 0.4, 0.1 };
    std::vector<int> group_sizes{ 0, 5, 1, 1, 5 }; 

    std::vector<double> output(means.size() * means.size());
    scran_markers::internal::compute_pairwise_simple_diff(means.data(), group_sizes.size(), 1, group_sizes.data(), output.data());

    for (size_t g = 0; g < means.size(); ++g) {
        for (size_t g2 = 0; g2 < means.size(); ++g2) {
            if (g == g2) {
                continue;
            } 
            
            double x = output[g * means.size() + g2];
            if (group_sizes[g] == 0 || group_sizes[g2] == 0) {
                EXPECT_TRUE(std::isnan(x));
            } else {
                double delta = means[g] - means[g2];
                EXPECT_FLOAT_EQ(x, delta);
            }
        }
    }
}

TEST(SimpleDiff, Blocked) {
    int nblocks = 2, ngroups = 4;
    std::vector<double> means{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
    std::vector<int> combo_sizes{ 10, 5, 12, 34, 15, 2, 3, 6 }; 

    std::vector<double> output(ngroups * ngroups);
    scran_markers::internal::compute_pairwise_simple_diff(means.data(), ngroups, nblocks, combo_sizes.data(), output.data());

    for (int g1 = 0; g1 < ngroups; ++g1) {
        for (int g2 = 0; g2 < ngroups; ++g2) {
            double totalnum = 0, totaldenom = 0;

            for (int b = 0; b < nblocks; ++b) {
                int offset1 = g1 + b * ngroups;
                int offset2 = g2 + b * ngroups;
                double d = means[offset1] - means[offset2];
                double w = combo_sizes[offset1] * combo_sizes[offset2];
                totalnum += d * w;
                totaldenom += w;
            }

            EXPECT_FLOAT_EQ(output[g1 * ngroups + g2], totalnum/totaldenom);
        }
    }
}

TEST(SimpleDiff, BlockedMissing) {
    double nan = std::numeric_limits<double>::quiet_NaN();
    int nblocks = 2, ngroups = 4;
    std::vector<double> means{ 0.0, nan, nan, nan, 0.2, 0.4, 0.6, 0.8 };
    std::vector<int> combo_sizes{ 0, 1, 0, 0, 34, 23, 6, 55 }; // exactly one group in the first block with non-zero size.

    std::vector<double> output(ngroups * ngroups);
    scran_markers::internal::compute_pairwise_simple_diff(means.data(), ngroups, nblocks, combo_sizes.data(), output.data());

    // Effectively excising the first block, as there is no comparison with a
    // non-zero product for the combined weight. We shouldn't see any NaNs
    // bleeding through, as the 0-size groups should be skipped.

    std::vector<double> output2(ngroups * ngroups);
    scran_markers::internal::compute_pairwise_simple_diff(means.data() + ngroups, ngroups, nblocks - 1, combo_sizes.data() + ngroups, output2.data());
    EXPECT_EQ(output, output2);
}
