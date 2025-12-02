#include "scran_tests/scran_tests.hpp"

#include <vector>
#include <algorithm>

#include "scran_markers/block_averages.hpp"

TEST(PrecomputedPairwiseWeights, SingleBlock) {
    // Single block falls back to all-unity weights to ensure that the weighted mean is a no-op.
    std::vector<double> weights{ 1, 2, 3, 4, 5};
    const int ngroups = weights.size();
    scran_markers::internal::PrecomputedPairwiseWeights<double> preweights(ngroups, 1, weights.data());

    for (int g1 = 0; g1 < ngroups; ++g1) {
        for (int g2 = 0; g2 < g1; ++g2) {
            const auto info = preweights.get(g1, g2);
            EXPECT_EQ(info.first[0], 1);
            EXPECT_EQ(info.second, 1);

            const auto info2 = preweights.get(g2, g1);
            EXPECT_EQ(info2.first[0], 1);
            EXPECT_EQ(info2.second, 1);
        }
    }
}

TEST(PrecomputedPairwiseWeights, TwoBlocks) {
    std::vector<double> weights{ 1, 2, 3, 4, 5, 6 };
    const int ngroups = 3;
    const int nblocks = 2;
    scran_markers::internal::PrecomputedPairwiseWeights<double> preweights(ngroups, nblocks, weights.data());

    for (int g1 = 0; g1 < ngroups; ++g1) {
        for (int g2 = 0; g2 < ngroups; ++g2) {
            if (g1 == g2) {
                continue;
            }

            const auto info = preweights.get(g1, g2);
            if ((g1 == 0 && g2 == 1) || (g1 == 1 && g2 == 0)) {
                EXPECT_EQ(info.first[0], 2);
                EXPECT_EQ(info.first[1], 20);
            } else if ((g1 == 0 && g2 == 2) || (g1 == 2 && g2 == 0)) {
                EXPECT_EQ(info.first[0], 3);
                EXPECT_EQ(info.first[1], 24);
            } else { // if ((g1 == 1 && g2 == 2) || (g1 == 2 && g2 == 1)) {
                EXPECT_EQ(info.first[0], 6);
                EXPECT_EQ(info.first[1], 30);
            }

            double combined = 0; 
            for (int b = 0; b < nblocks; ++b) {
                combined += info.first[b];
            }
            EXPECT_EQ(info.second, combined);
        }
    }
}

TEST(PrecomputedPairwiseWeights, StrippedBlock) {
    // Everything devolves to a weight of 1 as no pair is present in more than one block.
    std::vector<double> weights{ 1, 0, 0, 0, 5, 6 };
    const int ngroups = 3;
    const int nblocks = 2;
    scran_markers::internal::PrecomputedPairwiseWeights<double> preweights(ngroups, nblocks, weights.data());

    for (int g1 = 0; g1 < ngroups; ++g1) {
        for (int g2 = 0; g2 < ngroups; ++g2) {
            if (g1 == g2) {
                continue;
            }

            const auto info = preweights.get(g1, g2);
            if ((g1 == 0 && g2 == 1) || (g1 == 1 && g2 == 0)) {
                EXPECT_EQ(info.first[0], 0);
                EXPECT_EQ(info.first[1], 0);
                EXPECT_EQ(info.second, 0);
            } else if ((g1 == 0 && g2 == 2) || (g1 == 2 && g2 == 0)) {
                EXPECT_EQ(info.first[0], 0);
                EXPECT_EQ(info.first[1], 0);
                EXPECT_EQ(info.second, 0);
            } else { // if ((g1 == 1 && g2 == 2) || (g1 == 2 && g2 == 1)) {
                EXPECT_EQ(info.first[0], 0);
                EXPECT_EQ(info.first[1], 1);
                EXPECT_EQ(info.second, 1);
            }
        }
    }
}

TEST(BlockAverageInfo, Basic) {
    std::vector<double> weights{ 1, 2, 3, 4, 5, 6 };

    scran_markers::internal::BlockAverageInfo<double> info(weights);
    EXPECT_TRUE(info.use_mean());
    EXPECT_EQ(info.combo_weights(), weights);

    scran_markers::internal::BlockAverageInfo<double> qinfo(0.5);
    EXPECT_FALSE(qinfo.use_mean());
    EXPECT_EQ(qinfo.quantile(), 0.5);
}
