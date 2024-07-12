#include <gtest/gtest.h>

#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"

#include "scran_tests/scran_tests.hpp"
#include "scran_markers/score_markers_pairwise.hpp"

#include "utils.h"

/*********************************************/

// We compare against the reference to check that the account-keeping
// with respect to threading and threshold specification is correct.

class ScoreMarkersPairwiseUnblockedTest : public ::testing::TestWithParam<std::tuple<int, double, bool, int> > {
protected:
    inline static std::shared_ptr<tatami::Matrix<double, int> > dense_row, dense_column, sparse_row, sparse_column;

    static void SetUpTestSuite() {
        size_t nr = 898, nc = 76;
        dense_row.reset(
            new tatami::DenseRowMatrix<double, int>(
                nr,
                nc,
                scran_tests::simulate_vector(
                    nr * nc, 
                    []{
                        scran_tests::SimulationParameters sparam;
                        sparam.density = 0.1;
                        sparam.seed = 99998;
                        return sparam;
                    }()
                )
            )
        );

        dense_column = tatami::convert_to_dense(dense_row.get(), false);
        sparse_row = tatami::convert_to_compressed_sparse(dense_row.get(), true);
        sparse_column = tatami::convert_to_compressed_sparse(dense_row.get(), false);
    }

    static auto simple_reference(const tatami::Matrix<double, int>& mat, const int* group, double threshold) {
        size_t ngenes = mat.nrow();
        auto group_sizes = tatami_stats::tabulate_groups(group, mat.ncol());
        size_t ngroups = group_sizes.size();

        scran_markers::ScoreMarkersPairwiseResults<double> output;
        output.cohens_d.resize(ngroups * ngroups * ngenes);
        output.delta_mean = output.cohens_d;
        output.delta_detected = output.cohens_d;

        output.mean = tatami_stats::grouped_sums::by_row(&mat, group, tatami_stats::grouped_sums::Options());
        auto all_variances = tatami_stats::grouped_variances::by_row(&mat, group, tatami_stats::grouped_variances::Options());

        auto nonzero = tatami::make_DelayedUnaryIsometricOperation(
            tatami::wrap_shared_ptr(&mat), 
            tatami::DelayedUnaryIsometricCompareScalar<tatami::CompareOperation::NOT_EQUAL, double>(0)
        );
        output.detected = tatami_stats::grouped_sums::by_row(nonzero.get(), group, tatami_stats::grouped_sums::Options());

        for (size_t g = 0; g < ngroups; ++g) {
            double current = group_sizes[g];
            for (auto& r : output.mean[g]) {
                r /= current;
            }
            for (auto& r : output.detected[g]) {
                r /= current;
            }
        }

        scran_markers::internal::PrecomputedPairwiseWeights preweights(ngroups, 1, group_sizes.data());
        std::vector<double> means(ngroups), variances(ngroups), detected(ngroups);
        for (size_t r = 0; r < ngenes; ++r) {
            for (size_t g = 0; g < ngroups; ++g) {
                means[g] = output.mean[g][r];
                variances[g] = all_variances[g][r];
                detected[g] = output.detected[g][r];
            }

            size_t out_offset = r * ngroups * ngroups;
            scran_markers::internal::compute_pairwise_cohens_d(means.data(), variances.data(), ngroups, 1, preweights, threshold, output.cohens_d.data() + out_offset);
            scran_markers::internal::compute_pairwise_simple_diff(means.data(), ngroups, 1, preweights, output.delta_mean.data() + out_offset);
            scran_markers::internal::compute_pairwise_simple_diff(detected.data(), ngroups, 1, preweights, output.delta_detected.data() + out_offset);
        }

        return output;
    }

    static void compare_results(
        const scran_markers::ScoreMarkersPairwiseResults<double>& expected, 
        const scran_markers::ScoreMarkersPairwiseResults<double>& observed,
        bool include_auc) 
    {
        ASSERT_EQ(expected.mean.size(), observed.mean.size());
        for (size_t g = 0; g < expected.mean.size(); ++g) {
            scran_tests::compare_almost_equal(expected.mean[g], observed.mean[g]);
        }

        ASSERT_EQ(expected.detected.size(), observed.detected.size());
        for (size_t g = 0; g < expected.detected.size(); ++g) {
            scran_tests::compare_almost_equal(expected.detected[g], observed.detected[g]);
        }

        scran_tests::compare_almost_equal(expected.cohens_d, observed.cohens_d);
        scran_tests::compare_almost_equal(expected.delta_mean, observed.delta_mean);
        scran_tests::compare_almost_equal(expected.delta_detected, observed.delta_detected);

        if (include_auc) {
            scran_tests::compare_almost_equal(expected.auc, observed.auc);
        }
    }
};

TEST_P(ScoreMarkersPairwiseUnblockedTest, Reference) {
    auto param = GetParam();
    auto ngroups = std::get<0>(param);
    auto threshold = std::get<1>(param);
    auto auc = std::get<2>(param);
    auto nthreads = std::get<3>(param);

    std::vector<int> groupings = create_groupings(dense_row->ncol(), ngroups);

    scran_markers::ScoreMarkersPairwiseOptions opt;
    opt.threshold = threshold;
    opt.compute_auc = auc;
    auto ref = scran_markers::score_markers_pairwise(*dense_row, groupings.data(), opt);

    // Checking that all the values match up to the reference.
    if (nthreads == 1) {
        auto simple = simple_reference(*dense_row, groupings.data(), threshold);
        compare_results(ref, simple, /* include_auc = */ false);
    } else {
        opt.num_threads = nthreads;
        auto drres = scran_markers::score_markers_pairwise(*dense_row, groupings.data(), opt);
        compare_results(ref, drres, auc);
    }

    auto dcres = scran_markers::score_markers_pairwise(*dense_column, groupings.data(), opt);
    compare_results(ref, dcres, auc);

    auto srres = scran_markers::score_markers_pairwise(*sparse_row, groupings.data(), opt);
    compare_results(ref, srres, auc);

    auto scres = scran_markers::score_markers_pairwise(*sparse_column, groupings.data(), opt);
    compare_results(ref, scres, auc);
}

INSTANTIATE_TEST_SUITE_P(
    ScoreMarkersPairwiseUnblocked,
    ScoreMarkersPairwiseUnblockedTest,
    ::testing::Combine(
        ::testing::Values(2, 3, 4, 5), // number of clusters
        ::testing::Values(0, 0.5), // threshold to use
        ::testing::Values(false, true), // whether to compute AUCs.
        ::testing::Values(1, 3) // number of threads
    )
);

///*********************************************/
//
//class ScoreMarkersPairwiseBlockedTest : public ::testing::TestWithParam<std::tuple<int, int, scran::WeightPolicy, int> >, public DifferentialAnalysisTestCore {
//protected:
//    static constexpr size_t nrows = 100;
//    static constexpr size_t ncols = 50;
//
//    void SetUp() {
//        assemble(nrows, ncols);
//    }
//};
//
//TEST_P(ScoreMarkersPairwiseBlockedTest, VersusReference) {
//    auto param = GetParam();
//    auto ngroups = std::get<0>(param);
//    auto nblocks = std::get<1>(param);
//    auto policy = std::get<2>(param);
//    auto nthreads = std::get<3>(param);
//
//    scran::ScoreMarkersPairwise chd;
//    chd.set_num_threads(nthreads);
//    chd.set_block_weight_policy(policy);
//    scran::VariableBlockWeightParameters vparams{0, 10};
//    chd.set_variable_block_weight_parameters(vparams); // get some interesting variable weights, if we can.
//
//    auto groups = create_groupings(ncols, ngroups);
//    auto blocks = create_blocks(ncols, nblocks);
//
//    std::vector<double> ref_cohen(ngroups * ngroups * nrows);
//    auto ref_auc = ref_cohen;
//    auto ref_lfc = ref_cohen;
//    auto ref_delta_detected = ref_cohen;
//    std::vector<int> total_product_weights(ngroups * ngroups);
//
//    std::vector<std::vector<double> > ref_means(ngroups, std::vector<double>(nrows));
//    auto ref_detected = ref_means;
//    std::vector<int> total_group_weights(ngroups);
//
//    for (int b = 0; b < nblocks; ++b) {
//        std::vector<int> subset;
//        std::vector<int> subgroups;
//        for (int i = 0; i < ncols; ++i) {
//            if (blocks[i] == b) {
//                subset.push_back(i);
//                subgroups.push_back(groups[i]);
//            }
//        }
//
//        auto sub = tatami::make_DelayedSubset<1>(dense_row, std::move(subset));
//        auto res = chd.run(sub.get(), subgroups.data());
//
//        auto subcount = scran::tabulate_ids(subgroups.size(), subgroups.data());
//        std::vector<double> subweights = scran::compute_block_weights(subcount, policy, vparams);
//
//        for (size_t i = 0; i < nrows; ++i) {
//            for (int g1 = 0; g1 < ngroups; ++g1) {
//                for (int g2 = 0; g2 < ngroups; ++g2) {
//                    size_t offset = i * ngroups * ngroups + g1 * ngroups + g2;
//                    double weight = subweights[g1] * subweights[g2];
//                    ref_cohen[offset] += weight * res.cohen[offset];
//                    ref_lfc[offset] += weight * res.lfc[offset];
//                    ref_delta_detected[offset] += weight * res.delta_detected[offset];
//                    ref_auc[offset] += weight * res.auc[offset];
//                }
//            }
//
//            for (int g = 0; g < ngroups; ++g) {
//                ref_means[g][i] += res.means[g][i] * subweights[g];
//                ref_detected[g][i] += res.detected[g][i] * subweights[g];
//            }
//        }
//
//        for (int g1 = 0; g1 < ngroups; ++g1) {
//            total_group_weights[g1] += subweights[g1];
//            for (int g2 = 0; g2 < ngroups; ++g2) {
//                total_product_weights[g1 * ngroups + g2] += subweights[g1] * subweights[g2];
//            }
//        }
//    }
//
//    for (size_t i = 0; i < nrows; ++i) {
//        auto offset = i * ngroups * ngroups;
//
//        for (int g = 0; g < ngroups * ngroups; ++g) {
//            ref_cohen[offset + g] /= total_product_weights[g];
//            ref_lfc[offset + g] /= total_product_weights[g];
//            ref_auc[offset + g] /= total_product_weights[g];
//            ref_delta_detected[offset + g] /= total_product_weights[g];
//        }
//
//        for (int g = 0; g < ngroups; ++g) {
//            ref_means[g][i] /= total_group_weights[g];
//            ref_detected[g][i] /= total_group_weights[g];
//        }
//    }
//
//    // Alright, running all the tests.
//    auto res = chd.run_blocked(sparse_row.get(), groups.data(), blocks.data());
//    compare_almost_equal(ref_cohen, res.cohen);
//    compare_almost_equal(ref_lfc, res.lfc);
//    compare_almost_equal(ref_delta_detected, res.delta_detected);
//    compare_almost_equal(ref_auc, res.auc);
//
//    for (int g = 0; g < ngroups; ++g) {
//        compare_almost_equal(ref_means[g], res.means[g]);
//        compare_almost_equal(ref_detected[g], res.detected[g]);
//    }
//}
//
//INSTANTIATE_TEST_SUITE_P(
//    ScoreMarkersPairwiseBlocked,
//    ScoreMarkersPairwiseBlockedTest,
//    ::testing::Combine(
//        ::testing::Values(2, 3, 4, 5), // number of clusters
//        ::testing::Values(1, 2, 3), // number of blocks
//        ::testing::Values(scran::WeightPolicy::NONE, scran::WeightPolicy::EQUAL, scran::WeightPolicy::VARIABLE), // block weighting method.
//        ::testing::Values(1, 3) // number of threads
//    )
//);
//
///*********************************************/
//
//class ScoreMarkersPairwiseScenarioTest : public ::testing::Test, public DifferentialAnalysisTestCore {
//protected:
//    static constexpr size_t nrows = 100;
//    static constexpr size_t ncols = 50;
//
//    void SetUp() {
//        assemble(nrows, ncols);
//    }
//};
//
//TEST_F(ScoreMarkersPairwiseScenarioTest, Self) {
//    int copies = 3;
//
//    // Replicating the same matrix 'ngroups' times.
//    std::vector<std::shared_ptr<tatami::NumericMatrix> > stuff;
//    for (int i = 0; i < copies; ++i) {
//        stuff.push_back(dense_row);
//    }
//    auto combined = tatami::make_DelayedBind<1>(std::move(stuff));
//
//    // Creating two groups; second group can be larger than the first, to check
//    // for correct behavior w.r.t. imbalanced groups.
//    std::vector<int> groupings(ncols * copies);
//    std::fill(groupings.begin(), groupings.begin() + ncols, 0);
//    std::fill(groupings.begin() + ncols, groupings.end(), 1); 
//
//    scran::ScoreMarkersPairwise chd;
//    auto res = chd.run(combined.get(), groupings.data());
//
//    // All AUCs should be 0.5, all Cohen/LFC/delta-d's should be 0.
//    int ngroups = 2;
//    std::vector<double> cohen(ngroups * ngroups * nrows);
//    auto lfc = cohen, delta_detected = cohen;
//    std::vector<double> auc(cohen.size(), 0.5);
//
//    for (size_t g = 0; g < nrows; ++g) {
//        for (int l = 0; l < ngroups; ++l) {
//            size_t offset = g * ngroups * ngroups + l * ngroups + l;  
//            auc[offset] = 0;
//        }
//    }
//
//    compare_almost_equal(cohen, res.cohen);
//    compare_almost_equal(auc, res.auc);
//    compare_almost_equal(lfc, res.lfc);
//    compare_almost_equal(delta_detected, res.delta_detected);
//}
//
//TEST_F(ScoreMarkersPairwiseScenarioTest, Perfect) {
//    int ngroups = 5;
//    std::vector<int> groupings = create_groupings(ncols, ngroups);
//
//    int nrows = 33;
//    std::vector<double> pretend;
//    for (int r = 0; r < nrows; ++r) {
//        pretend.insert(pretend.end(), groupings.begin(), groupings.end());
//    }
//
//    tatami::DenseRowMatrix<double, int> mat(nrows, groupings.size(), std::move(pretend));
//    scran::ScoreMarkersPairwise chd;
//    auto res = chd.run(&mat, groupings.data());
//
//    for (size_t g = 0; g < nrows; ++g) {
//        for (int l = 0; l < ngroups; ++l) {
//            for (int l2 = 0; l2 < ngroups; ++l2) {
//                if (l == l2) {
//                    continue;
//                }
//
//                size_t offset = g * ngroups * ngroups + l * ngroups + l2;  
//                EXPECT_EQ(res.lfc[offset], l - l2);
//                EXPECT_EQ(res.delta_detected[offset], (l > 0) - (l2 > 0));
//                EXPECT_EQ(res.auc[offset], static_cast<double>(l > l2));
//                EXPECT_TRUE(std::isinf(res.cohen[offset]));
//                EXPECT_EQ(res.cohen[offset] > 0, l > l2);
//            }
//        }
//    }
//}
//
//TEST_F(ScoreMarkersPairwiseScenarioTest, Thresholds) {
//    int ngroups = 3;
//    std::vector<int> groupings = create_groupings(ncols, ngroups);
//
//    scran::ScoreMarkersPairwise chd;
//    auto ref = chd.run(dense_row.get(), groupings.data());
//    auto out = chd.set_threshold(1).run(dense_row.get(), groupings.data());
//    EXPECT_EQ(ref.lfc, out.lfc);
//    EXPECT_EQ(ref.delta_detected, out.delta_detected);
//
//    bool some_diff = false;
//    for (size_t g = 0; g < nrows; ++g) {
//        for (int l = 0; l < ngroups; ++l) {
//            for (int l2 = 0; l2 < ngroups; ++l2) {
//                if (l == l2) {
//                    continue;
//                }
//
//                // Threshold should have some effect for cohen.
//                size_t offset = g * ngroups * ngroups + l * ngroups + l2;  
//                EXPECT_TRUE(ref.cohen[offset] > out.cohen[offset]);
//
//                // '>' is not guaranteed due to imprecision with ranks... but (see below).
//                EXPECT_TRUE(ref.auc[offset] >= out.auc[offset]); 
//            }
//        }
//    }
//
//    // There should be at least some difference here.
//    EXPECT_NE(ref.auc, out.auc);
//}
//
//TEST_F(ScoreMarkersPairwiseScenarioTest, Missing) {
//    int ngroups = 3;
//    std::vector<int> groupings = create_groupings(ncols, ngroups);
//
//    scran::ScoreMarkersPairwise chd;
//    auto ref = chd.run(dense_row.get(), groupings.data());
//
//    // Zero is effectively the missing group here.
//    for (auto& g : groupings) {
//        ++g;
//    }
//    auto lost = chd.run(dense_row.get(), groupings.data());
//
//    // Everything should be NaN.
//    int ngroups_p1 = ngroups + 1;
//    for (size_t g = 0; g < nrows; ++g) {
//        for (int l2 = 1; l2 < ngroups_p1; ++l2) {
//            // For the comparisons from group 0 to the others.
//            size_t offset = g * ngroups_p1 * ngroups_p1 + l2;  
//            EXPECT_TRUE(std::isnan(lost.cohen[offset]));
//            EXPECT_TRUE(std::isnan(lost.lfc[offset]));
//            EXPECT_TRUE(std::isnan(lost.delta_detected[offset]));
//            EXPECT_TRUE(std::isnan(lost.auc[offset]));
//
//            // For the comparisons in the other direction.
//            offset = g * ngroups_p1 * ngroups_p1 + l2 * ngroups_p1;  
//            EXPECT_TRUE(std::isnan(lost.cohen[offset]));
//            EXPECT_TRUE(std::isnan(lost.lfc[offset]));
//            EXPECT_TRUE(std::isnan(lost.delta_detected[offset]));
//            EXPECT_TRUE(std::isnan(lost.auc[offset]));
//        }
//    }
//
//    // Other metrics should be the same as usual.
//    for (size_t g = 0; g < nrows; ++g) {
//        for (int l = 0; l < ngroups; ++l) {
//            size_t ref_offset = g * ngroups * ngroups + l * ngroups;  
//            size_t lost_offset = g * ngroups_p1 * ngroups_p1 + (l + 1) * ngroups_p1 + 1; // skip group 0, and also the NaN in the comparison against group 0.
//
//            EXPECT_EQ(
//                std::vector<double>(ref.cohen.begin() + ref_offset, ref.cohen.begin() + ref_offset + ngroups), 
//                std::vector<double>(lost.cohen.begin() + lost_offset, lost.cohen.begin() + lost_offset + ngroups)
//            );
//
//            EXPECT_EQ(
//                std::vector<double>(ref.lfc.begin() + ref_offset, ref.lfc.begin() + ref_offset + ngroups), 
//                std::vector<double>(lost.lfc.begin() + lost_offset, lost.lfc.begin() + lost_offset + ngroups)
//            );
//
//            EXPECT_EQ(
//                std::vector<double>(ref.delta_detected.begin() + ref_offset, ref.delta_detected.begin() + ref_offset + ngroups), 
//                std::vector<double>(lost.delta_detected.begin() + lost_offset, lost.delta_detected.begin() + lost_offset + ngroups)
//            );
//
//            EXPECT_EQ(
//                std::vector<double>(ref.cohen.begin() + ref_offset, ref.cohen.begin() + ref_offset + ngroups), 
//                std::vector<double>(lost.cohen.begin() + lost_offset, lost.cohen.begin() + lost_offset + ngroups)
//            );
//        }
//    }
//}
//
//TEST_F(ScoreMarkersPairwiseScenarioTest, BlockConfounded) {
//    auto ngroups = 4;
//    std::vector<int> groupings = create_groupings(ncols, ngroups);
//
//    // Block is fully confounded with one group.
//    std::vector<int> blocks(ncols);
//    for (size_t c = 0; c < ncols; ++c) {
//        blocks[c] = groupings[c] == 0;
//    }
//
//    scran::ScoreMarkersPairwise chd;
//    auto comres = chd.run_blocked(dense_row.get(), groupings.data(), blocks.data());
//
//    // First group should only be NaN's.
//    size_t ngenes = dense_row->nrow();
//    for (size_t g = 0; g < nrows; ++g) {
//        for (int l2 = 1; l2 < ngroups; ++l2) {
//            // For the comparisons from group 0 to the others.
//            size_t offset = g * ngroups * ngroups + l2;  
//            EXPECT_TRUE(std::isnan(comres.cohen[offset]));
//            EXPECT_TRUE(std::isnan(comres.lfc[offset]));
//            EXPECT_TRUE(std::isnan(comres.delta_detected[offset]));
//            EXPECT_TRUE(std::isnan(comres.auc[offset]));
//
//            // For the comparisons in the other direction.
//            offset = g * ngroups * ngroups + l2 * ngroups;  
//            EXPECT_TRUE(std::isnan(comres.cohen[offset]));
//            EXPECT_TRUE(std::isnan(comres.lfc[offset]));
//            EXPECT_TRUE(std::isnan(comres.delta_detected[offset]));
//            EXPECT_TRUE(std::isnan(comres.auc[offset]));
//        }
//    }
//
//    // Excluding the confounded group and running on the remaining samples.
//    std::vector<int> subgroups;
//    std::vector<int> keep;
//    for (size_t c = 0; c < ncols; ++c) {
//        auto g = groupings[c];
//        if (g != 0) {
//            subgroups.push_back(g - 1);
//            keep.push_back(c);
//        }
//    }
//
//    auto sub = tatami::make_DelayedSubset<1>(dense_row, std::move(keep));
//    auto ref = chd.run(sub.get(), subgroups.data()); 
//    int ngroups_m1 = ngroups - 1;
//
//    for (size_t g = 0; g < nrows; ++g) {
//        for (int l = 0; l < ngroups_m1; ++l) {
//            size_t ref_offset = g * ngroups_m1 * ngroups_m1 + l * ngroups_m1;  
//            size_t comres_offset = g * ngroups * ngroups + (l + 1) * ngroups + 1; // skip group 0 as well as the NaN in the comparison against group 0.
//
//            EXPECT_EQ(
//                std::vector<double>(ref.cohen.begin() + ref_offset, ref.cohen.begin() + ref_offset + ngroups_m1), 
//                std::vector<double>(comres.cohen.begin() + comres_offset, comres.cohen.begin() + comres_offset + ngroups_m1)
//            );
//
//            EXPECT_EQ(
//                std::vector<double>(ref.lfc.begin() + ref_offset, ref.lfc.begin() + ref_offset + ngroups_m1), 
//                std::vector<double>(comres.lfc.begin() + comres_offset, comres.lfc.begin() + comres_offset + ngroups_m1)
//            );
//
//            EXPECT_EQ(
//                std::vector<double>(ref.delta_detected.begin() + ref_offset, ref.delta_detected.begin() + ref_offset + ngroups_m1), 
//                std::vector<double>(comres.delta_detected.begin() + comres_offset, comres.delta_detected.begin() + comres_offset + ngroups_m1)
//            );
//
//            EXPECT_EQ(
//                std::vector<double>(ref.cohen.begin() + ref_offset, ref.cohen.begin() + ref_offset + ngroups_m1), 
//                std::vector<double>(comres.cohen.begin() + comres_offset, comres.cohen.begin() + comres_offset + ngroups_m1)
//            );
//        }
//    }
//}
