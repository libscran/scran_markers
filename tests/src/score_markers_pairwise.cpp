#include <gtest/gtest.h>

#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"

#include "scran_tests/scran_tests.hpp"
#include "scran_markers/score_markers_pairwise.hpp"

#include "utils.h"

void compare_results(
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

/*********************************************/

// We compare against the reference to check that the account-keeping
// with respect to threading and threshold specification is correct.

class ScoreMarkersPairwiseUnblockedTest : public ::testing::TestWithParam<std::tuple<int, double, bool, int> > {
protected:
    inline static std::shared_ptr<tatami::Matrix<double, int> > dense_row, dense_column, sparse_row, sparse_column;

    static void SetUpTestSuite() {
        size_t nr = 898, nc = 176;
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

        if (auc) {
            for (int r = 0; r < dense_row->nrow(); ++r) {
                for (int g1 = 0; g1 < ngroups; ++g1) {
                    for (int g2 = 0; g2 < ngroups; ++g2) {
                        auto val = ref.auc[r * ngroups * ngroups + g1 * ngroups + g2];
                        if (g1 == g2) {
                            EXPECT_EQ(val, 0);
                        } else {
                            EXPECT_GE(val, 0); // checking correct bounds.
                            EXPECT_LE(val, 1);
                        }
                    }
                }
            }
        }

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

/*********************************************/

class ScoreMarkersPairwiseBlockedTest : public ::testing::TestWithParam<std::tuple<int, int, bool, scran_blocks::WeightPolicy, int> > {
protected:
    inline static std::shared_ptr<tatami::Matrix<double, int> > dense_row, dense_column, sparse_row, sparse_column;

    static void SetUpTestSuite() {
        size_t nr = 398, nc = 157; // use a prime number of columns to check for non-equal weights.
        dense_row.reset(
            new tatami::DenseRowMatrix<double, int>(
                nr,
                nc,
                scran_tests::simulate_vector(
                    nr * nc, 
                    []{
                        scran_tests::SimulationParameters sparam;
                        sparam.density = 0.1;
                        sparam.seed = 999998;
                        return sparam;
                    }()
                )
            )
        );

        dense_column = tatami::convert_to_dense(dense_row.get(), false);
        sparse_row = tatami::convert_to_compressed_sparse(dense_row.get(), true);
        sparse_column = tatami::convert_to_compressed_sparse(dense_row.get(), false);
    }

    static auto blocked_reference(const tatami::Matrix<double, int>& mat, const int* group, const int* blocks, const scran_markers::ScoreMarkersPairwiseOptions& opt) {
        scran_markers::ScoreMarkersPairwiseResults<double> output;

        size_t ngenes = mat.nrow();
        size_t ngroups = tatami_stats::total_groups(group, mat.ncol());
        output.mean.reserve(ngroups);
        output.detected.reserve(ngroups);
        for (size_t g = 0; g < ngroups; ++g) {
            output.mean.emplace_back(ngenes);
            output.detected.emplace_back(ngenes);
        }

        output.cohens_d.resize(ngroups * ngroups * ngenes);
        output.delta_mean = output.cohens_d;
        output.delta_detected = output.cohens_d;
        if (opt.compute_auc) {
            output.auc = output.cohens_d;
        }

        int nblocks = tatami_stats::total_groups(blocks, mat.ncol());
        std::vector<double> total_group_weights(ngroups);
        std::vector<double> total_product_weights(ngroups * ngroups);

        for (int b = 0; b < nblocks; ++b) {
            // Slicing the matrix.
            std::vector<int> subset;
            std::vector<int> subgroups;
            int ncols = mat.ncol();
            for (int i = 0; i < ncols; ++i) {
                if (blocks[i] == b) {
                    subset.push_back(i);
                    subgroups.push_back(group[i]);
                }
            }

            auto sub = tatami::make_DelayedSubset(dense_row, std::move(subset), false);
            auto res = scran_markers::score_markers_pairwise(*sub, subgroups.data(), opt);
            auto subcount = tatami_stats::tabulate_groups(subgroups.data(), subgroups.size());
            std::vector<double> subweights = scran_blocks::compute_weights(subcount, opt.block_weight_policy, opt.variable_block_weight_parameters);

            for (size_t i = 0; i < ngenes; ++i) {
                for (size_t g1 = 0; g1 < ngroups; ++g1) {
                    for (size_t g2 = 0; g2 < ngroups; ++g2) {
                        size_t offset = i * ngroups * ngroups + g1 * ngroups + g2;
                        double weight = subweights[g1] * subweights[g2];
                        output.cohens_d[offset] += weight * res.cohens_d[offset];
                        output.delta_mean[offset] += weight * res.delta_mean[offset];
                        output.delta_detected[offset] += weight * res.delta_detected[offset];
                        if (opt.compute_auc) {
                            output.auc[offset] += weight * res.auc[offset];
                        }
                    }
                }

                for (size_t g = 0; g < ngroups; ++g) {
                    output.mean[g][i] += res.mean[g][i] * subweights[g];
                    output.detected[g][i] += res.detected[g][i] * subweights[g];
                }
            }

            for (size_t g1 = 0; g1 < ngroups; ++g1) {
                total_group_weights[g1] += subweights[g1];
                for (size_t g2 = 0; g2 < ngroups; ++g2) {
                    total_product_weights[g1 * ngroups + g2] += subweights[g1] * subweights[g2];
                }
            }
        }

        for (size_t i = 0; i < ngenes; ++i) {
            auto offset = i * ngroups * ngroups;

            for (size_t g1 = 0; g1 < ngroups; ++g1) {
                for (size_t g2 = 0; g2 < ngroups; ++g2) {
                    size_t from = g1 * ngroups + g2;
                    size_t to = offset + from;
                    output.cohens_d[to] /= total_product_weights[from];
                    output.delta_mean[to] /= total_product_weights[from];
                    output.delta_detected[to] /= total_product_weights[from];
                    if (opt.compute_auc) {
                        output.auc[to] /= total_product_weights[from];
                    }
                }
            }

            for (size_t g = 0; g < ngroups; ++g) {
                output.mean[g][i] /= total_group_weights[g];
                output.detected[g][i] /= total_group_weights[g];
            }
        }

        return output;
    }
};

TEST_P(ScoreMarkersPairwiseBlockedTest, VersusReference) {
    auto param = GetParam();
    auto ngroups = std::get<0>(param);
    auto nblocks = std::get<1>(param);
    auto auc = std::get<2>(param);
    auto policy = std::get<3>(param);
    auto nthreads = std::get<4>(param);

    scran_markers::ScoreMarkersPairwiseOptions opt;
    opt.block_weight_policy = policy;
    opt.compute_auc = auc;

    auto ncols = dense_row->ncol();
    auto groups = create_groupings(ncols, ngroups);
    auto blocks = create_blocks(ncols, nblocks);
    auto ref = scran_markers::score_markers_pairwise_blocked(*dense_row, groups.data(), blocks.data(), opt);

    if (nthreads == 1) {
        auto simple = blocked_reference(*dense_row, groups.data(), blocks.data(), opt);
        compare_results(ref, simple, auc);
    } else {
        opt.num_threads = nthreads;
        auto drres = scran_markers::score_markers_pairwise_blocked(*dense_row, groups.data(), blocks.data(), opt);
        compare_results(ref, drres, auc);
    }

    auto dcres = scran_markers::score_markers_pairwise_blocked(*dense_column, groups.data(), blocks.data(), opt);
    compare_results(ref, dcres, auc);

    auto srres = scran_markers::score_markers_pairwise_blocked(*sparse_row, groups.data(), blocks.data(), opt);
    compare_results(ref, srres, auc);

    auto scres = scran_markers::score_markers_pairwise_blocked(*sparse_column, groups.data(), blocks.data(), opt);
    compare_results(ref, scres, auc);
}

INSTANTIATE_TEST_SUITE_P(
    ScoreMarkersPairwiseBlocked,
    ScoreMarkersPairwiseBlockedTest,
    ::testing::Combine(
        ::testing::Values(2, 3, 4, 5), // number of clusters
        ::testing::Values(1, 2, 3), // number of blocks
        ::testing::Values(false, true), // whether to compute AUC or not.
        ::testing::Values(scran_blocks::WeightPolicy::NONE, scran_blocks::WeightPolicy::EQUAL), // block weighting method.
        ::testing::Values(1, 3) // number of threads
    )
);

/*********************************************/

TEST(ScoreMarkersPairwiseScenarios, Self) {
    int nrows = 132, ncols = 97;
    std::shared_ptr<tatami::NumericMatrix> mat(
        new tatami::DenseRowMatrix<double, int>(
            nrows,
            ncols,
            scran_tests::simulate_vector(
                nrows * ncols,
                []{
                    scran_tests::SimulationParameters sparam;
                    sparam.seed = 69;
                    return sparam;
                }()
            )
        )
    );

    // Replicating the same matrix 3 times.
    const int copies = 3;
    std::vector<std::shared_ptr<tatami::NumericMatrix> > stuff;
    for (int i = 0; i < copies; ++i) {
        stuff.push_back(mat);
    }
    auto combined = tatami::make_DelayedBind(std::move(stuff), false);

    // Creating two groups; second group can be larger than the first, to check
    // for correct behavior w.r.t. imbalanced groups.
    std::vector<int> groupings(ncols * copies);
    std::fill(groupings.begin(), groupings.begin() + ncols, 0);
    std::fill(groupings.begin() + ncols, groupings.end(), 1); 

    scran_markers::ScoreMarkersPairwiseOptions opt;
    auto res = scran_markers::score_markers_pairwise(*combined, groupings.data(), opt);

    // All AUCs should be 0.5, all Cohen/LFC/delta-d's should be 0.
    int ngroups = 2;
    std::vector<double> cohen(ngroups * ngroups * nrows);
    auto lfc = cohen, delta_detected = cohen;
    std::vector<double> auc(cohen.size(), 0.5);

    for (int g = 0; g < nrows; ++g) {
        for (int l = 0; l < ngroups; ++l) {
            size_t offset = g * ngroups * ngroups + l * ngroups + l;  
            auc[offset] = 0;
        }
    }

    scran_tests::compare_almost_equal(cohen, res.cohens_d);
    scran_tests::compare_almost_equal(auc, res.auc);
    scran_tests::compare_almost_equal(lfc, res.delta_mean);
    scran_tests::compare_almost_equal(delta_detected, res.delta_detected);
}

TEST(ScoreMarkersPairwiseScenarios, Perfect) {
    int ngroups = 5;
    int ncols = 71;
    std::vector<int> groupings = create_groupings(ncols, ngroups);

    int nrows = 33;
    std::vector<double> pretend;
    for (int r = 0; r < nrows; ++r) {
        pretend.insert(pretend.end(), groupings.begin(), groupings.end());
    }

    tatami::DenseRowMatrix<double, int> mat(nrows, groupings.size(), std::move(pretend));
    scran_markers::ScoreMarkersPairwiseOptions opt;
    auto res = scran_markers::score_markers_pairwise(mat, groupings.data(), opt);

    for (int g = 0; g < nrows; ++g) {
        for (int l = 0; l < ngroups; ++l) {
            for (int l2 = 0; l2 < ngroups; ++l2) {
                if (l == l2) {
                    continue;
                }

                size_t offset = g * ngroups * ngroups + l * ngroups + l2;  
                EXPECT_EQ(res.delta_mean[offset], l - l2);
                EXPECT_EQ(res.delta_detected[offset], (l > 0) - (l2 > 0));
                EXPECT_EQ(res.auc[offset], static_cast<double>(l > l2));
                EXPECT_TRUE(std::isinf(res.cohens_d[offset]));
                EXPECT_EQ(res.cohens_d[offset] > 0, l > l2);
            }
        }
    }
}

TEST(ScoreMarkersPairwiseScenarios, Thresholds) {
    int nrows = 67, ncols = 91;
    tatami::DenseRowMatrix<double, int> mat(
        nrows,
        ncols,
        scran_tests::simulate_vector(
            nrows * ncols,
            []{
                scran_tests::SimulationParameters sparam;
                sparam.seed = 696969;
                return sparam;
            }()
        )
    );

    int ngroups = 3;
    std::vector<int> groupings = create_groupings(ncols, ngroups);
    scran_markers::ScoreMarkersPairwiseOptions opt;
    auto ref = scran_markers::score_markers_pairwise(mat, groupings.data(), opt);

    opt.threshold = 1;
    auto out = scran_markers::score_markers_pairwise(mat, groupings.data(), opt);
    EXPECT_EQ(ref.delta_mean, out.delta_mean);
    EXPECT_EQ(ref.delta_detected, out.delta_detected);

    for (int g = 0; g < nrows; ++g) {
        for (int l = 0; l < ngroups; ++l) {
            for (int l2 = 0; l2 < ngroups; ++l2) {
                if (l == l2) {
                    continue;
                }

                // Threshold should have some effect for cohen.
                size_t offset = g * ngroups * ngroups + l * ngroups + l2;  
                EXPECT_TRUE(ref.cohens_d[offset] > out.cohens_d[offset]);

                // '>' is not guaranteed due to imprecision with ranks... but (see below).
                EXPECT_TRUE(ref.auc[offset] >= out.auc[offset]); 
            }
        }
    }

    // There should be at least some difference here.
    EXPECT_NE(ref.auc, out.auc);
}

TEST(ScoreMarkersPairwiseScenarios, Missing) {
    int nrows = 144, ncols = 109;
    tatami::DenseRowMatrix<double, int> mat(
        nrows,
        ncols,
        scran_tests::simulate_vector(
            nrows * ncols,
            []{
                scran_tests::SimulationParameters sparam;
                sparam.seed = 696969;
                return sparam;
            }()
        )
    );

    int ngroups = 4;
    std::vector<int> groupings = create_groupings(ncols, ngroups);
    scran_markers::ScoreMarkersPairwiseOptions opt;
    auto ref = scran_markers::score_markers_pairwise(mat, groupings.data(), opt);

    // Zero is effectively the missing group here.
    for (auto& g : groupings) {
        ++g;
    }
    auto lost = scran_markers::score_markers_pairwise(mat, groupings.data(), opt);

    // Everything should be NaN.
    int ngroups_p1 = ngroups + 1;
    for (int g = 0; g < nrows; ++g) {
        for (int l2 = 1; l2 < ngroups_p1; ++l2) {
            // For the comparisons from group 0 to the others.
            size_t offset = g * ngroups_p1 * ngroups_p1 + l2;  
            EXPECT_TRUE(std::isnan(lost.cohens_d[offset]));
            EXPECT_TRUE(std::isnan(lost.delta_mean[offset]));
            EXPECT_TRUE(std::isnan(lost.delta_detected[offset]));
            EXPECT_TRUE(std::isnan(lost.auc[offset]));

            // For the comparisons in the other direction.
            offset = g * ngroups_p1 * ngroups_p1 + l2 * ngroups_p1;  
            EXPECT_TRUE(std::isnan(lost.cohens_d[offset]));
            EXPECT_TRUE(std::isnan(lost.delta_mean[offset]));
            EXPECT_TRUE(std::isnan(lost.delta_detected[offset]));
            EXPECT_TRUE(std::isnan(lost.auc[offset]));
        }
    }

    // Other metrics should be the same as usual.
    for (int g = 0; g < nrows; ++g) {
        for (int l = 0; l < ngroups; ++l) {
            size_t ref_offset = g * ngroups * ngroups + l * ngroups;  
            size_t lost_offset = g * ngroups_p1 * ngroups_p1 + (l + 1) * ngroups_p1 + 1; // skip group 0, and also the NaN in the comparison against group 0.

            EXPECT_EQ(scran_tests::vector_n(ref.cohens_d.data() + ref_offset, ngroups), scran_tests::vector_n(lost.cohens_d.data() + lost_offset, ngroups));
            EXPECT_EQ(scran_tests::vector_n(ref.auc.data() + ref_offset, ngroups), scran_tests::vector_n(lost.auc.data() + lost_offset, ngroups));
            EXPECT_EQ(scran_tests::vector_n(ref.delta_mean.data() + ref_offset, ngroups), scran_tests::vector_n(lost.delta_mean.data() + lost_offset, ngroups));
            EXPECT_EQ(scran_tests::vector_n(ref.delta_detected.data() + ref_offset, ngroups), scran_tests::vector_n(lost.delta_detected.data() + lost_offset, ngroups));
        }
    }
}

TEST(ScoreMarkersPairwiseScenarios, BlockConfounded) {
    int nrows = 198, ncols = 99;
    std::shared_ptr<tatami::Matrix<double, int> > mat(
        new tatami::DenseRowMatrix<double, int>(
            nrows,
            ncols,
            scran_tests::simulate_vector(
                nrows * ncols,
                []{
                    scran_tests::SimulationParameters sparam;
                    sparam.seed = 69696969;
                    return sparam;
                }()
            )
        )
    );

    int ngroups = 4;
    std::vector<int> groupings = create_groupings(ncols, ngroups);

    // Block is fully confounded with one group.
    std::vector<int> blocks(ncols);
    for (int c = 0; c < ncols; ++c) {
        blocks[c] = groupings[c] == 0;
    }

    scran_markers::ScoreMarkersPairwiseOptions opt;
    auto comres = scran_markers::score_markers_pairwise_blocked(*mat, groupings.data(), blocks.data(), opt);

    // First group should only be NaN's.
    for (int g = 0; g < nrows; ++g) {
        for (int l2 = 1; l2 < ngroups; ++l2) {
            // For the comparisons from group 0 to the others.
            size_t offset = g * ngroups * ngroups + l2;  
            EXPECT_TRUE(std::isnan(comres.cohens_d[offset]));
            EXPECT_TRUE(std::isnan(comres.delta_mean[offset]));
            EXPECT_TRUE(std::isnan(comres.delta_detected[offset]));
            EXPECT_TRUE(std::isnan(comres.auc[offset]));

            // For the comparisons in the other direction.
            offset = g * ngroups * ngroups + l2 * ngroups;  
            EXPECT_TRUE(std::isnan(comres.cohens_d[offset]));
            EXPECT_TRUE(std::isnan(comres.delta_mean[offset]));
            EXPECT_TRUE(std::isnan(comres.delta_detected[offset]));
            EXPECT_TRUE(std::isnan(comres.auc[offset]));
        }
    }

    // Excluding the confounded group and running on the remaining samples.
    std::vector<int> subgroups;
    std::vector<int> keep;
    for (int c = 0; c < ncols; ++c) {
        auto g = groupings[c];
        if (g != 0) {
            subgroups.push_back(g - 1);
            keep.push_back(c);
        }
    }

    auto sub = tatami::make_DelayedSubset(mat, std::move(keep), false);
    auto ref = scran_markers::score_markers_pairwise(*sub, subgroups.data(), opt);
    int ngroups_m1 = ngroups - 1;

    for (int g = 0; g < nrows; ++g) {
        for (int l = 0; l < ngroups_m1; ++l) {
            size_t ref_offset = g * ngroups_m1 * ngroups_m1 + l * ngroups_m1;  
            size_t comres_offset = g * ngroups * ngroups + (l + 1) * ngroups + 1; // skip group 0 as well as the NaN in the comparison against group 0.

            EXPECT_EQ(scran_tests::vector_n(ref.cohens_d.data() + ref_offset, ngroups_m1), scran_tests::vector_n(comres.cohens_d.data() + comres_offset, ngroups_m1));
            EXPECT_EQ(scran_tests::vector_n(ref.delta_mean.data() + ref_offset, ngroups_m1), scran_tests::vector_n(comres.delta_mean.data() + comres_offset, ngroups_m1));
            EXPECT_EQ(scran_tests::vector_n(ref.delta_detected.data() + ref_offset, ngroups_m1), scran_tests::vector_n(comres.delta_detected.data() + comres_offset, ngroups_m1));
            EXPECT_EQ(scran_tests::vector_n(ref.auc.data() + ref_offset, ngroups_m1), scran_tests::vector_n(comres.auc.data() + comres_offset, ngroups_m1));
        }
    }
}
