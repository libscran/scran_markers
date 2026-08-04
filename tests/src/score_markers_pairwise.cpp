#include "scran_tests/scran_tests.hpp"
#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"

#include "scran_markers/score_markers_pairwise.hpp"

#include "utils.h"

static void compare_averages(const std::vector<std::vector<double> >& res, const std::vector<std::vector<double> >& other) {
    const int ngroups = res.size();
    ASSERT_EQ(ngroups, other.size());
    for (int l = 0; l < ngroups; ++l) {
        scran_tests::compare_almost_equal(res[l], other[l]);
    }
}

static void compare_results(
    const scran_markers::ScoreMarkersPairwiseResults<double>& expected, 
    const scran_markers::ScoreMarkersPairwiseResults<double>& observed,
    bool include_auc) 
{
    compare_averages(expected.mean, observed.mean);
    compare_averages(expected.detected, observed.detected);

    scran_tests::compare_almost_equal(expected.cohens_d, observed.cohens_d);
    scran_tests::compare_almost_equal(expected.delta_mean, observed.delta_mean);
    scran_tests::compare_almost_equal(expected.delta_detected, observed.delta_detected);

    if (include_auc) {
        scran_tests::compare_almost_equal(expected.auc, observed.auc);
    }
}

/*********************************************/

class ScoreMarkersPairwiseTest : public ::testing::TestWithParam<std::tuple<int, bool, int> > {
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
                        scran_tests::SimulateVectorParameters sparam;
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
};

static auto simple_reference(const tatami::Matrix<double, int>& mat, const int* group, std::size_t num_groups, double threshold) {
    const int num_genes = mat.nrow();
    auto group_sizes = scran_markers::tabulate_groups(mat.ncol(), group, num_groups);

    scran_markers::ScoreMarkersPairwiseResults<double> output;
    output.cohens_d.resize(num_groups * num_groups * num_genes);
    output.delta_mean = output.cohens_d;
    output.delta_detected = output.cohens_d;

    auto var_out = tatami_stats::group_variance(true, mat, group, num_groups, {});
    output.mean = var_out.mean;
    const auto& all_variances = var_out.variance;

    auto nonzero = tatami::DelayedUnaryIsometricOperation<double, double, int>(
        tatami::wrap_shared_ptr(&mat), 
        std::make_shared<tatami::DelayedUnaryIsometricCompareScalarHelper<tatami::CompareOperation::NOT_EQUAL, double, double, int, int> >(0)
    );
    output.detected = tatami_stats::group_sum(true, nonzero, group, num_groups, {});

    for (std::size_t g = 0; g < num_groups; ++g) {
        double current = group_sizes[g];
        for (auto& r : output.detected[g]) {
            r /= current;
        }
    }

    std::vector<double> combo_weights(group_sizes.begin(), group_sizes.end());
    scran_markers::internal::PrecomputedPairwiseWeights preweights(num_groups, 1, combo_weights.data());

    std::vector<double> means(num_groups), variances(num_groups), detected(num_groups);
    for (int r = 0; r < num_genes; ++r) {
        for (std::size_t g = 0; g < num_groups; ++g) {
            means[g] = output.mean[g][r];
            variances[g] = all_variances[g][r];
            detected[g] = output.detected[g][r];
        }

        size_t out_offset = r * num_groups * num_groups;
        scran_markers::internal::compute_pairwise_cohens_d_blockmean(means.data(), variances.data(), num_groups, 1, threshold, preweights, output.cohens_d.data() + out_offset);
        scran_markers::internal::compute_pairwise_simple_diff_blockmean(means.data(), num_groups, 1, preweights, output.delta_mean.data() + out_offset);
        scran_markers::internal::compute_pairwise_simple_diff_blockmean(detected.data(), num_groups, 1, preweights, output.delta_detected.data() + out_offset);
    }

    return output;
}

TEST_P(ScoreMarkersPairwiseTest, Reference) {
    auto param = GetParam();
    const auto ngroups = std::get<0>(param);
    const auto auc = std::get<1>(param);
    const auto nthreads = std::get<2>(param);

    auto groupings = create_interleaved_factor(dense_row->ncol(), ngroups);

    scran_markers::ScoreMarkersPairwiseOptions opt;
    opt.compute_auc = auc;
    auto ref = scran_markers::score_markers_pairwise(*dense_row, groupings.data(), ngroups, opt);

    // Checking that all the values match up to the reference.
    if (nthreads == 1) {
        auto simple = simple_reference(*dense_row, groupings.data(), ngroups, 0.0);
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
        auto drres = scran_markers::score_markers_pairwise(*dense_row, groupings.data(), ngroups, opt);
        compare_results(ref, drres, auc);
    }

    // Testing the other matrix representations.
    {
        auto dcres = scran_markers::score_markers_pairwise(*dense_column, groupings.data(), ngroups, opt);
        compare_results(ref, dcres, auc);

        auto srres = scran_markers::score_markers_pairwise(*sparse_row, groupings.data(), ngroups, opt);
        compare_results(ref, srres, auc);

        auto scres = scran_markers::score_markers_pairwise(*sparse_column, groupings.data(), ngroups, opt);
        compare_results(ref, scres, auc);
    }
}

INSTANTIATE_TEST_SUITE_P(
    ScoreMarkersPairwise,
    ScoreMarkersPairwiseTest,
    ::testing::Combine(
        ::testing::Values(2, 5), // number of groups
        ::testing::Values(false, true), // whether to compute AUCs.
        ::testing::Values(1, 3) // number of threads
    )
);

/*********************************************/

class ScoreMarkersPairwiseBlockedTest : public ::testing::TestWithParam<std::tuple<int, bool, scran_blocks::WeightPolicy, int> > {
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
                        scran_tests::SimulateVectorParameters sparam;
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
};

TEST_P(ScoreMarkersPairwiseBlockedTest, SingleBlock) {
    auto param = GetParam();
    auto ngroups = std::get<0>(param);
    auto auc = std::get<1>(param);
    auto policy = std::get<2>(param);
    auto nthreads = std::get<3>(param);

    scran_markers::ScoreMarkersPairwiseOptions opt;
    opt.block_weight_policy = policy;
    opt.compute_auc = auc;
    opt.num_threads = nthreads;

    const auto ncols = dense_row->ncol();
    auto groups = create_interleaved_factor(ncols, ngroups);
    std::vector<int> blocks(ncols);

    auto ref = scran_markers::score_markers_pairwise(*dense_row, groups.data(), ngroups, opt);

    // Should get the same result with 1 block.
    {
        auto drres = scran_markers::score_markers_pairwise_blocked(*dense_row, groups.data(), ngroups, blocks.data(), 1, opt);
        compare_results(ref, drres, auc);

        auto dcres = scran_markers::score_markers_pairwise_blocked(*dense_column, groups.data(), ngroups, blocks.data(), 1, opt);
        compare_results(ref, dcres, auc);

        auto srres = scran_markers::score_markers_pairwise_blocked(*sparse_row, groups.data(), ngroups, blocks.data(), 1, opt);
        compare_results(ref, srres, auc);

        auto scres = scran_markers::score_markers_pairwise_blocked(*sparse_column, groups.data(), ngroups, blocks.data(), 1, opt);
        compare_results(ref, scres, auc);
    }

    // Same results for quantile.
    {
        auto qopt = opt;
        qopt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;

        auto drres = scran_markers::score_markers_pairwise_blocked(*dense_row, groups.data(), ngroups, blocks.data(), 1, qopt);
        compare_results(ref, drres, auc);

        auto dcres = scran_markers::score_markers_pairwise_blocked(*dense_column, groups.data(), ngroups, blocks.data(), 1, qopt);
        compare_results(ref, dcres, auc);

        auto srres = scran_markers::score_markers_pairwise_blocked(*sparse_row, groups.data(), ngroups, blocks.data(), 1, qopt);
        compare_results(ref, srres, auc);

        auto scres = scran_markers::score_markers_pairwise_blocked(*sparse_column, groups.data(), ngroups, blocks.data(), 1, qopt);
        compare_results(ref, scres, auc);
    }
}

static scran_markers::ScoreMarkersPairwiseResults<double> allocate_output(std::size_t num_genes, std::size_t ngroups, bool do_auc) {
    scran_markers::ScoreMarkersPairwiseResults<double> output;
    output.mean.reserve(ngroups);
    output.detected.reserve(ngroups);
    for (std::size_t g = 0; g < ngroups; ++g) {
        output.mean.emplace_back(num_genes);
        output.detected.emplace_back(num_genes);
    }

    const std::size_t full_size = ngroups * ngroups * num_genes;
    output.cohens_d.resize(full_size);
    output.delta_mean.resize(full_size);
    output.delta_detected.resize(full_size);
    if (do_auc) {
        output.auc.resize(full_size);
    }

    return output;
}

static auto blocked_reference_mean(
    const tatami::Matrix<double, int>& mat,
    const int* group,
    const std::size_t num_groups,
    const int* blocks,
    const std::size_t num_blocks,
    const scran_markers::ScoreMarkersPairwiseOptions& opt
) {
    const int num_genes = mat.nrow();
    auto output = allocate_output(num_genes, num_groups, opt.compute_auc);

    std::vector<double> total_group_weights(num_groups);
    std::vector<double> total_product_weights(num_groups * num_groups);

    for (std::size_t b = 0; b < num_blocks; ++b) {
        std::vector<int> subset, subgroups;
        int ncols = mat.ncol();
        for (int i = 0; i < ncols; ++i) {
            if (sanisizer::is_equal(blocks[i], b)) {
                subset.push_back(i);
                subgroups.push_back(group[i]);
            }
        }

        auto sub = tatami::make_DelayedSubset(tatami::wrap_shared_ptr(&mat), std::move(subset), false);
        auto res = scran_markers::score_markers_pairwise(*sub, subgroups.data(), num_groups, opt);
        auto subcount = scran_markers::tabulate_groups(subgroups.size(), subgroups.data(), num_groups);
        auto subweights = scran_blocks::compute_weights(subcount, opt.block_weight_policy, opt.variable_block_weight_parameters);

        for (int i = 0; i < num_genes; ++i) {
            for (std::size_t g1 = 0; g1 < num_groups; ++g1) {
                for (std::size_t g2 = 0; g2 < num_groups; ++g2) {
                    size_t offset = i * num_groups * num_groups + g1 * num_groups + g2;
                    double weight = subweights[g1] * subweights[g2];
                    output.cohens_d[offset] += weight * res.cohens_d[offset];
                    output.delta_mean[offset] += weight * res.delta_mean[offset];
                    output.delta_detected[offset] += weight * res.delta_detected[offset];
                    if (opt.compute_auc) {
                        output.auc[offset] += weight * res.auc[offset];
                    }
                }
            }

            for (std::size_t g = 0; g < num_groups; ++g) {
                output.mean[g][i] += res.mean[g][i] * subweights[g];
                output.detected[g][i] += res.detected[g][i] * subweights[g];
            }
        }

        for (std::size_t g1 = 0; g1 < num_groups; ++g1) {
            total_group_weights[g1] += subweights[g1];
            for (size_t g2 = 0; g2 < num_groups; ++g2) {
                total_product_weights[g1 * num_groups + g2] += subweights[g1] * subweights[g2];
            }
        }
    }

    for (int i = 0; i < num_genes; ++i) {
        auto offset = i * num_groups * num_groups;
        for (std::size_t g1 = 0; g1 < num_groups; ++g1) {
            for (std::size_t g2 = 0; g2 < num_groups; ++g2) {
                std::size_t from = g1 * num_groups + g2;
                std::size_t to = offset + from;
                output.cohens_d[to] /= total_product_weights[from];
                output.delta_mean[to] /= total_product_weights[from];
                output.delta_detected[to] /= total_product_weights[from];
                if (opt.compute_auc) {
                    output.auc[to] /= total_product_weights[from];
                }
            }
        }

        for (std::size_t g = 0; g < num_groups; ++g) {
            output.mean[g][i] /= total_group_weights[g];
            output.detected[g][i] /= total_group_weights[g];
        }
    }

    return output;
}

TEST_P(ScoreMarkersPairwiseBlockedTest, ReferenceMean) {
    auto param = GetParam();
    auto ngroups = std::get<0>(param);
    auto auc = std::get<1>(param);
    auto policy = std::get<2>(param);
    auto nthreads = std::get<3>(param);

    scran_markers::ScoreMarkersPairwiseOptions opt;
    opt.block_weight_policy = policy;
    opt.compute_auc = auc;

    const auto ncols = dense_row->ncol();
    auto groups = create_interleaved_factor(ncols, ngroups);
    const int num_blocks = 3;
    auto blocks = create_contiguous_factor(ncols, num_blocks);

    auto ref = scran_markers::score_markers_pairwise_blocked(*dense_row, groups.data(), ngroups, blocks.data(), num_blocks, opt);

    if (nthreads == 1) {
        auto simple = blocked_reference_mean(*dense_row, groups.data(), ngroups, blocks.data(), num_blocks, opt);
        compare_results(ref, simple, auc);
    } else {
        opt.num_threads = nthreads;
        auto drres = scran_markers::score_markers_pairwise_blocked(*dense_row, groups.data(), ngroups, blocks.data(), num_blocks, opt);
        compare_results(ref, drres, auc);
    }

    auto dcres = scran_markers::score_markers_pairwise_blocked(*dense_column, groups.data(), ngroups, blocks.data(), num_blocks, opt);
    compare_results(ref, dcres, auc);

    auto srres = scran_markers::score_markers_pairwise_blocked(*sparse_row, groups.data(), ngroups, blocks.data(), num_blocks, opt);
    compare_results(ref, srres, auc);

    auto scres = scran_markers::score_markers_pairwise_blocked(*sparse_column, groups.data(), ngroups, blocks.data(), num_blocks, opt);
    compare_results(ref, scres, auc);
}

static auto blocked_reference_quantile(
    const tatami::Matrix<double, int>& mat,
    const int* group,
    const std::size_t num_groups,
    const int* blocks,
    const std::size_t num_blocks,
    const scran_markers::ScoreMarkersPairwiseOptions& opt
) {
    const int num_genes = mat.nrow();
    auto output = allocate_output(num_genes, num_groups, opt.compute_auc);

    // Indexing goes: group, gene, blocks.
    std::vector<std::vector<std::vector<double> > > qbuffers_mean(num_groups), qbuffers_det(num_groups);
    for (std::size_t g = 0; g < num_groups; ++g) {
        qbuffers_mean[g].resize(num_genes);
        qbuffers_det[g].resize(num_genes);
    }

    // Indexing goes: group 1, group 2, gene, blocks.
    std::vector<std::vector<std::vector<std::vector<double> > > > qbuffers_cohen(num_groups), qbuffers_dmean(num_groups), qbuffers_ddet(num_groups), qbuffers_auc;
    if (opt.compute_auc) {
        qbuffers_auc.resize(num_groups);
    }

    for (std::size_t g1 = 0; g1 < num_groups; ++g1) {
        qbuffers_cohen[g1].resize(num_groups);
        qbuffers_dmean[g1].resize(num_groups);
        qbuffers_ddet[g1].resize(num_groups);
        if (opt.compute_auc) {
            qbuffers_auc[g1].resize(num_groups);
        }

        for (std::size_t g2 = 0; g2 < num_groups; ++g2) {
            qbuffers_cohen[g1][g2].resize(num_genes);
            qbuffers_dmean[g1][g2].resize(num_genes);
            qbuffers_ddet[g1][g2].resize(num_genes);
            if (opt.compute_auc) {
                qbuffers_auc[g1][g2].resize(num_genes);
            }
        }
    }

    for (std::size_t b = 0; b < num_blocks; ++b) {
        std::vector<int> subset, subgroups;
        const int ncols = mat.ncol();
        for (int i = 0; i < ncols; ++i) {
            if (sanisizer::is_equal(blocks[i], b)) {
                subset.push_back(i);
                subgroups.push_back(group[i]);
            }
        }

        auto sub = tatami::make_DelayedSubset(tatami::wrap_shared_ptr(&mat), std::move(subset), false);
        auto res = scran_markers::score_markers_pairwise(*sub, subgroups.data(), num_groups, opt);

        for (int i = 0; i < num_genes; ++i) {
            for (std::size_t g1 = 0; g1 < num_groups; ++g1) {
                for (std::size_t g2 = 0; g2 < num_groups; ++g2) {
                    std::size_t offset = i * num_groups * num_groups + g1 * num_groups + g2;
                    qbuffers_cohen[g1][g2][i].push_back(res.cohens_d[offset]);
                    qbuffers_dmean[g1][g2][i].push_back(res.delta_mean[offset]);
                    qbuffers_ddet[g1][g2][i].push_back(res.delta_detected[offset]);
                    if (opt.compute_auc) {
                        qbuffers_auc[g1][g2][i].push_back(res.auc[offset]);
                    }
                }
            }

            for (std::size_t g = 0; g < num_groups; ++g) {
                qbuffers_mean[g][i].push_back(res.mean[g][i]);
                qbuffers_det[g][i].push_back(res.detected[g][i]);
            }
        }
    }

    quickstats::SingleQuantileFixedNumber<double> qcalc(num_blocks, opt.block_quantile);
    for (int i = 0; i < num_genes; ++i) {
        std::size_t offset = i * num_groups * num_groups;

        for (std::size_t g1 = 0; g1 < num_groups; ++g1) {
            for (std::size_t g2 = 0; g2 < num_groups; ++g2) {
                std::size_t from = g1 * num_groups + g2;
                std::size_t to = offset + from;
                output.cohens_d[to] = qcalc(qbuffers_cohen[g1][g2][i].data());
                output.delta_mean[to] = qcalc(qbuffers_dmean[g1][g2][i].data());
                output.delta_detected[to] = qcalc(qbuffers_ddet[g1][g2][i].data());
                if (opt.compute_auc) {
                    output.auc[to] = qcalc(qbuffers_auc[g1][g2][i].data());
                }
            }
        }

        for (std::size_t g = 0; g < num_groups; ++g) {
            output.mean[g][i] = qcalc(qbuffers_mean[g][i].data());
            output.detected[g][i] = qcalc(qbuffers_det[g][i].data());
        }
    }

    return output;
}

TEST_P(ScoreMarkersPairwiseBlockedTest, ReferenceQuantile) {
    auto param = GetParam();
    auto ngroups = std::get<0>(param);
    auto auc = std::get<1>(param);
    auto policy = std::get<2>(param);
    auto nthreads = std::get<3>(param);

    const auto ncols = dense_row->ncol();
    const auto num_blocks = 3;
    auto groups = create_contiguous_factor(ncols, ngroups);
    auto blocks = create_interleaved_factor(ncols, num_blocks);

    scran_markers::ScoreMarkersPairwiseOptions opt;
    opt.compute_auc = auc;
    opt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
    opt.block_weight_policy = policy; // shouldn't really matter as quantiles are unweighted, but whatever.
    auto ref = scran_markers::score_markers_pairwise_blocked(*dense_row, groups.data(), ngroups, blocks.data(), num_blocks, opt);

    if (nthreads == 1) {
        auto simple = blocked_reference_quantile(*dense_row, groups.data(), ngroups, blocks.data(), num_blocks, opt);
        compare_results(ref, simple, auc);
    } else {
        opt.num_threads = nthreads;
        auto rres = scran_markers::score_markers_pairwise_blocked(*dense_row, groups.data(), ngroups, blocks.data(), num_blocks, opt);
        compare_results(ref, rres, auc);
    }

    auto dcres = scran_markers::score_markers_pairwise_blocked(*dense_column, groups.data(), ngroups, blocks.data(), num_blocks, opt);
    compare_results(ref, dcres, auc);

    auto srres = scran_markers::score_markers_pairwise_blocked(*sparse_row, groups.data(), ngroups, blocks.data(), num_blocks, opt);
    compare_results(ref, srres, auc);

    auto scres = scran_markers::score_markers_pairwise_blocked(*sparse_column, groups.data(), ngroups, blocks.data(), num_blocks, opt);
    compare_results(ref, scres, auc);
}

INSTANTIATE_TEST_SUITE_P(
    ScoreMarkersPairwiseBlocked,
    ScoreMarkersPairwiseBlockedTest,
    ::testing::Combine(
        ::testing::Values(2, 5), // number of groups
        ::testing::Values(false, true), // whether to compute AUC or not.
        ::testing::Values(scran_blocks::WeightPolicy::NONE, scran_blocks::WeightPolicy::EQUAL), // block weighting method.
        ::testing::Values(1, 3) // number of threads
    )
);

/*********************************************/

TEST(ScoreMarkersPairwise, VersusSelf) {
    int nrows = 132, ncols = 97;
    std::shared_ptr<tatami::NumericMatrix> mat(
        new tatami::DenseRowMatrix<double, int>(
            nrows,
            ncols,
            scran_tests::simulate_vector(
                nrows * ncols,
                []{
                    scran_tests::SimulateVectorParameters sparam;
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

    // Creating two groups where the second group is larger than the first.
    // This aims to check correct behavior w.r.t. imbalanced groups.
    std::vector<int> groupings(ncols * copies);
    std::fill(groupings.begin(), groupings.begin() + ncols, 0);
    std::fill(groupings.begin() + ncols, groupings.end(), 1); 

    scran_markers::ScoreMarkersPairwiseOptions opt;
    auto res = scran_markers::score_markers_pairwise(*combined, groupings.data(), 2, opt);

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

    // Quantile should give the same results for a single block.
    auto qopt = opt;
    qopt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
    auto qres = scran_markers::score_markers_pairwise(*combined, groupings.data(), 2, qopt);
    compare_results(res, qres, true);
}

TEST(ScoreMarkersPairwise, PerfectSeparation) {
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
    auto res = scran_markers::score_markers_pairwise(mat, groupings.data(), ngroups, opt);

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

    // Quantile should give the same results for a single block.
    auto qopt = opt;
    qopt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
    auto qres = scran_markers::score_markers_pairwise(mat, groupings.data(), ngroups, qopt);
    compare_results(res, qres, true);
}

TEST(ScoreMarkersPairwise, Thresholds) {
    int nrows = 67, ncols = 91;
    tatami::DenseRowMatrix<double, int> mat(
        nrows,
        ncols,
        scran_tests::simulate_vector(
            nrows * ncols,
            []{
                scran_tests::SimulateVectorParameters sparam;
                sparam.seed = 696969;
                return sparam;
            }()
        )
    );

    const double threshold = 0.3;
    int ngroups = 3;
    std::vector<int> groupings = create_groupings(ncols, ngroups);
    auto simple = simple_reference(mat, groupings.data(), ngroups, threshold);

    scran_markers::ScoreMarkersPairwiseOptions opt;
    opt.threshold = threshold;
    auto out = scran_markers::score_markers_pairwise(mat, groupings.data(), ngroups, opt);
    compare_results(out, simple, /* include_auc = */ false);

    auto nothresh = scran_markers::score_markers_pairwise(mat, groupings.data(), ngroups, {});
    EXPECT_EQ(nothresh.delta_mean, out.delta_mean);
    EXPECT_EQ(nothresh.delta_detected, out.delta_detected);

    for (int g = 0; g < nrows; ++g) {
        for (int l = 0; l < ngroups; ++l) {
            for (int l2 = 0; l2 < ngroups; ++l2) {
                if (l == l2) {
                    continue;
                }

                // Threshold should have some effect for cohen.
                size_t offset = g * ngroups * ngroups + l * ngroups + l2;  
                EXPECT_TRUE(nothresh.cohens_d[offset] > out.cohens_d[offset]);

                // '>' is not guaranteed due to imprecision with ranks... but (see below).
                EXPECT_TRUE(nothresh.auc[offset] >= out.auc[offset]); 
            }
        }
    }

    // There should be at least some difference here.
    EXPECT_NE(nothresh.auc, out.auc);

    // Quantile should give the same results for a single block.
    auto qopt = opt;
    qopt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
    auto qout = scran_markers::score_markers_pairwise(mat, groupings.data(), ngroups, qopt);
    compare_results(out, qout, true);
}

static std::vector<double> populate_expected_effect_array_with_lost_groups(
    const int NR,
    const std::vector<int>& present,
    const std::vector<int>& lost,
    const std::vector<double>& source
) {
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    const std::size_t old_ngroups = present.size();
    const std::size_t new_ngroups = old_ngroups + lost.size();
    std::vector<double> expected(new_ngroups * new_ngroups * NR, nan);

    for (int r = 0; r < NR; ++r) {
        for (std::size_t i1 = 0; i1 < old_ngroups; ++i1) {
            for (std::size_t i2 = 0; i2 < old_ngroups; ++i2) {
                const auto in_pos = sanisizer::nd_offset<std::size_t>(i2, old_ngroups, i1, old_ngroups, r);
                const auto out_pos = sanisizer::nd_offset<std::size_t>(present[i2], new_ngroups, present[i1], new_ngroups, r);
                expected[out_pos] = source[in_pos];
            }
        }
        // Set self-comparisons to zero.
        for (const auto l : lost) {
            expected[sanisizer::nd_offset<std::size_t>(l, new_ngroups, l, new_ngroups, r)] = 0;
        }
    }

    return expected;
}

TEST(ScoreMarkersPairwise, EmptyGroups) {
    int nrows = 144, ncols = 109;
    tatami::DenseRowMatrix<double, int> mat(
        nrows,
        ncols,
        scran_tests::simulate_vector(
            nrows * ncols,
            []{
                scran_tests::SimulateVectorParameters sparam;
                sparam.seed = 696969;
                return sparam;
            }()
        )
    );

    int ngroups = 4;
    std::vector<int> groupings = create_groupings(ncols, ngroups);
    scran_markers::ScoreMarkersPairwiseOptions opt;
    auto ref = scran_markers::score_markers_pairwise(mat, groupings.data(), ngroups, opt);

    // First and last groups are empty.
    for (auto& g : groupings) {
        ++g;
    }
    const int ngroups_p2 = ngroups + 2;
    auto lost = scran_markers::score_markers_pairwise(mat, groupings.data(), ngroups_p2, opt);

    for (int g = 0; g < nrows; ++g) {
        EXPECT_TRUE(std::isnan(lost.mean[0][g]));
        EXPECT_TRUE(std::isnan(lost.detected[0][g]));
        EXPECT_TRUE(std::isnan(lost.mean[ngroups_p2 - 1][g]));
        EXPECT_TRUE(std::isnan(lost.detected[ngroups_p2 - 1][g]));
    }

    std::vector<int> present(ngroups);
    std::iota(present.begin(), present.end(), 1);
    std::vector<int> lost_groups{ 0, ngroups_p2 - 1 };

    auto expected_delta_mean = populate_expected_effect_array_with_lost_groups(nrows, present, lost_groups, ref.delta_mean);
    scran_tests::compare_almost_equal(expected_delta_mean, lost.delta_mean);

    auto expected_delta_detected = populate_expected_effect_array_with_lost_groups(nrows, present, lost_groups, ref.delta_detected);
    scran_tests::compare_almost_equal(expected_delta_detected, lost.delta_detected);

    auto expected_cohens_d = populate_expected_effect_array_with_lost_groups(nrows, present, lost_groups, ref.cohens_d);
    scran_tests::compare_almost_equal(expected_cohens_d, lost.cohens_d);

    auto expected_auc = populate_expected_effect_array_with_lost_groups(nrows, present, lost_groups, ref.auc);
    scran_tests::compare_almost_equal(expected_auc, lost.auc);

    // Quantile should give the same results for a single block.
    auto qopt = opt;
    qopt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
    auto qlost = scran_markers::score_markers_pairwise(mat, groupings.data(), ngroups_p2, qopt);
    compare_results(lost, qlost, true);
}

TEST(ScoreMarkersPairwise, BlockConfounded) {
    int nrows = 198, ncols = 99;
    std::shared_ptr<tatami::Matrix<double, int> > mat(
        new tatami::DenseRowMatrix<double, int>(
            nrows,
            ncols,
            scran_tests::simulate_vector(
                nrows * ncols,
                []{
                    scran_tests::SimulateVectorParameters sparam;
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
        blocks[c] = (groupings[c] == 0);
    }

    scran_markers::ScoreMarkersPairwiseOptions opt;
    auto comres = scran_markers::score_markers_pairwise_blocked(*mat, groupings.data(), ngroups, blocks.data(), 2, opt);

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
    int ngroups_m1 = ngroups - 1;
    auto ref = scran_markers::score_markers_pairwise(*sub, subgroups.data(), ngroups_m1, opt);

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

    // Quantile should give the same results, as there's basically only one block;
    // the second block is fully confounded.
    auto qopt = opt;
    qopt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
    auto qcomres = scran_markers::score_markers_pairwise_blocked(*mat, groupings.data(), ngroups, blocks.data(), 2, qopt);
    compare_results(comres, qcomres, true);
}

/*********************************************/

static void check_one_at_a_time(
    const tatami::Matrix<double, int>& mat,
    const std::vector<int>& groupings,
    const std::size_t ngroups,
    const scran_markers::ScoreMarkersPairwiseResults<double>& ref
) {
    // Only the group mean.
    {
        scran_markers::ScoreMarkersPairwiseOptions opt;
        opt.compute_group_detected = false;
        opt.compute_cohens_d = false;
        opt.compute_auc = false;
        opt.compute_delta_mean = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_pairwise<double>(mat, groupings.data(), ngroups, opt);
        compare_averages(ref.mean, alt.mean);
        EXPECT_TRUE(alt.detected.empty());
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only the group detected proportions.
    {
        scran_markers::ScoreMarkersPairwiseOptions opt;
        opt.compute_group_mean = false;
        opt.compute_cohens_d = false;
        opt.compute_auc = false;
        opt.compute_delta_mean = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_pairwise<double>(mat, groupings.data(), ngroups, opt);
        compare_averages(ref.detected, alt.detected);
        EXPECT_TRUE(alt.mean.empty());
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only Cohen's d.
    {
        scran_markers::ScoreMarkersPairwiseOptions opt;
        opt.compute_group_mean = false;
        opt.compute_group_detected = false;
        opt.compute_auc = false;
        opt.compute_delta_mean = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_pairwise<double>(mat, groupings.data(), ngroups, opt);
        scran_tests::compare_almost_equal(alt.cohens_d, ref.cohens_d);
        EXPECT_TRUE(alt.mean.empty());
        EXPECT_TRUE(alt.detected.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only AUC.
    {
        scran_markers::ScoreMarkersPairwiseOptions opt;
        opt.compute_group_mean = false;
        opt.compute_group_detected = false;
        opt.compute_cohens_d = false;
        opt.compute_delta_mean = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_pairwise<double>(mat, groupings.data(), ngroups, opt);
        scran_tests::compare_almost_equal(alt.auc, ref.auc);
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only delta-mean.
    {
        scran_markers::ScoreMarkersPairwiseOptions opt;
        opt.compute_group_mean = false;
        opt.compute_group_detected = false;
        opt.compute_cohens_d = false;
        opt.compute_auc = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_pairwise<double>(mat, groupings.data(), ngroups, opt);
        scran_tests::compare_almost_equal(alt.delta_mean, ref.delta_mean);
        EXPECT_TRUE(alt.mean.empty());
        EXPECT_TRUE(alt.detected.empty());
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only delta-detected.
    {
        scran_markers::ScoreMarkersPairwiseOptions opt;
        opt.compute_group_mean = false;
        opt.compute_group_detected = false;
        opt.compute_cohens_d = false;
        opt.compute_auc = false;
        opt.compute_delta_mean = false;

        auto alt = scran_markers::score_markers_pairwise<double>(mat, groupings.data(), ngroups, opt);
        scran_tests::compare_almost_equal(alt.delta_detected, ref.delta_detected);
        EXPECT_TRUE(alt.mean.empty());
        EXPECT_TRUE(alt.detected.empty());
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
    }
}

TEST(ScoreMarkersPairwise, OneAtATime) {
    int nr = 128, nc = 302;
    auto dense_row = std::make_unique<tatami::DenseRowMatrix<double, int> >(
        nr,
        nc,
        scran_tests::simulate_vector(
            nr * nc, 
            []{
                scran_tests::SimulateVectorParameters sparam;
                sparam.density = 0.2;
                sparam.seed = 96;
                return sparam;
            }()
        )
    );

    auto dense_column = tatami::convert_to_dense(dense_row.get(), false);
    auto sparse_row = tatami::convert_to_compressed_sparse(dense_row.get(), true);
    auto sparse_column = tatami::convert_to_compressed_sparse(dense_row.get(), false);

    const int ngroups = 3;
    std::vector<int> groupings = create_groupings(nc, ngroups);
    auto ref = scran_markers::score_markers_pairwise<double>(*dense_row, groupings.data(), ngroups, {});

    check_one_at_a_time(*dense_row, groupings, ngroups, ref);
    check_one_at_a_time(*sparse_row, groupings, ngroups, ref);
    check_one_at_a_time(*dense_column, groupings, ngroups, ref);
    check_one_at_a_time(*sparse_column, groupings, ngroups, ref);
}
