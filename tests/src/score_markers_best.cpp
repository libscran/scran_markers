#include "scran_tests/scran_tests.hpp"
#include "tatami/tatami.hpp"
#include "topicks/topicks.hpp"

#include "scran_markers/score_markers_best.hpp"
#include "scran_markers/score_markers_pairwise.hpp"
#include "scran_markers/summarize_effects.hpp"

#include "utils.h"

class ScoreMarkersBestTestCore {
protected:
    static void compare_averages(const std::vector<std::vector<double> >& res, const std::vector<std::vector<double> >& other) {
        const int ngroups = res.size();
        ASSERT_EQ(ngroups, other.size());
        for (int l = 0; l < ngroups; ++l) {
            scran_tests::compare_almost_equal(res[l], other[l]);
        }
    }

    void compare_best(
        const std::vector<std::vector<std::vector<std::pair<int, double> > > >& left,
        const std::vector<std::vector<std::vector<std::pair<int, double> > > >& right
    ) {
        ASSERT_EQ(left.size(), right.size());
        const int ngroups = left.size();
        for (int g1 = 0; g1 < ngroups; ++g1) {
            ASSERT_EQ(left[g1].size(), right[g1].size());
            for (int g2 = 0; g2 < ngroups; ++g2) {
                const auto n = left[g1][g2].size();
                ASSERT_EQ(n, right[g1][g2].size());

                for (std::size_t i = 0; i < n; ++i) {
                    const auto& lval = left[g1][g2][i];
                    const auto& rval = right[g1][g2][i];
                    EXPECT_EQ(lval.first, rval.first);
                    scran_tests::compare_almost_equal(lval.second, rval.second);
                }
            }
        }
    }

    void compare_best(
        const scran_markers::ScoreMarkersBestResults<double, int>& left,
        const scran_markers::ScoreMarkersBestResults<double, int>& right
    ) {
        compare_best(left.cohens_d, right.cohens_d);
        compare_best(left.auc, right.auc);
        compare_best(left.delta_mean, right.delta_mean);
        compare_best(left.delta_detected, right.delta_detected);
    }

    void compare_best_to_pairwise(
        const std::vector<std::vector<std::vector<std::pair<int, double> > > >& best,
        const std::vector<double>& effects,
        int ngenes,
        int ngroups,
        int top,
        bool larger,
        bool keep_ties,
        const std::optional<double>& bound
    ) {
        std::vector<double> buffer(ngenes);

        topicks::PickTopGenesOptions<double> opt;
        opt.check_nan = true;
        opt.keep_ties = keep_ties;
        if (bound.has_value()) {
            opt.bound = *bound;
        }

        for (int g1 = 0; g1 < ngroups; ++g1) {
            for (int g2 = 0; g2 < ngroups; ++g2) {
                const auto& curbest = best[g1][g2];
                if (g1 == g2) {
                    EXPECT_TRUE(curbest.empty());
                    continue;
                }

                for (int r = 0; r < ngenes; ++r) {
                    buffer[r] = effects[sanisizer::nd_offset<std::size_t>(g2, ngroups, g1, ngroups, r)];
                }

                const auto curtop = topicks::pick_top_genes_index(ngenes, buffer.data(), top, larger, opt);
                std::vector<std::pair<int, double> > expected;
                expected.reserve(curtop.size());
                for (auto t : curtop) {
                    expected.emplace_back(t, buffer[t]);
                }
                std::sort(expected.begin(), expected.end(), [&](const std::pair<int, double>& left, const std::pair<int, double>& right) -> bool {
                    if (left.second == right.second) {
                        return left.first < right.first;
                    } else if (larger) {
                        return left.second > right.second;
                    } else {
                        return left.second < right.second;
                    }
                });

                EXPECT_EQ(expected, curbest);
            }
        }
    }
};

/*********************************************/

class ScoreMarkersBestTest : public ScoreMarkersBestTestCore, public ::testing::TestWithParam<std::tuple<int, bool, int, bool, bool, bool, int> > {
protected:
    inline static std::shared_ptr<tatami::Matrix<double, int> > dense_row, dense_column, sparse_row, sparse_column;

    static void SetUpTestSuite() {
        size_t nr = 428, nc = 177;
        dense_row.reset(
            new tatami::DenseRowMatrix<double, int>(
                nr,
                nc,
                scran_tests::simulate_vector(
                    nr * nc, 
                    []{
                        scran_tests::SimulateVectorParameters sparam;
                        sparam.density = 0.2;
                        sparam.seed = 6900;
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

TEST_P(ScoreMarkersBestTest, Basic) {
    auto param = GetParam();
    auto ngroups = std::get<0>(param);
    bool do_auc = std::get<1>(param);
    int top = std::get<2>(param);
    bool larger = std::get<3>(param);
    bool use_bound = std::get<4>(param);
    bool keep_ties = std::get<5>(param);
    auto nthreads = std::get<6>(param);

    std::vector<int> groupings = create_groupings(dense_row->ncol(), ngroups);
    size_t ngenes = dense_row->nrow();

    scran_markers::ScoreMarkersBestOptions opt;
    opt.largest_cohens_d = larger;
    opt.largest_delta_detected = larger;
    opt.largest_delta_mean = larger;
    opt.largest_auc = larger;

    if (!use_bound) {
        opt.threshold_cohens_d.reset();
        opt.threshold_delta_detected.reset();
        opt.threshold_delta_mean.reset();
        opt.threshold_auc.reset();
    }

    opt.keep_ties = keep_ties;
    opt.compute_auc = do_auc; // false, if we want to check the running implementations.
    auto ref = scran_markers::score_markers_best<double>(*dense_row, groupings.data(), top, opt);

    if (nthreads == 1) {
        // Comparing score_markers_best against score_markers_pairwise + topicks::pick_top_genes.
        // The latter is less mind-bending but requires holding a large 3D matrix in memory.
        scran_markers::ScoreMarkersPairwiseOptions popt;
        popt.compute_auc = do_auc;
        auto pairres = scran_markers::score_markers_pairwise(*dense_row, groupings.data(), popt);
        compare_averages(ref.mean, pairres.mean);
        compare_averages(ref.detected, pairres.detected);

        compare_best_to_pairwise(ref.cohens_d, pairres.cohens_d, ngenes, ngroups, top, opt.largest_cohens_d, keep_ties, opt.threshold_cohens_d);
        compare_best_to_pairwise(ref.delta_mean, pairres.delta_mean, ngenes, ngroups, top, opt.largest_delta_mean, keep_ties, opt.threshold_delta_mean);
        compare_best_to_pairwise(ref.delta_detected, pairres.delta_detected, ngenes, ngroups, top, opt.largest_delta_detected, keep_ties, opt.threshold_delta_detected);
        if (do_auc) {
            compare_best_to_pairwise(ref.auc, pairres.auc, ngenes, ngroups, top, opt.largest_auc, keep_ties, opt.threshold_auc);
        }

    } else {
        opt.num_threads = nthreads;
        auto dr = scran_markers::score_markers_best<double>(*dense_row, groupings.data(), top, opt);
        compare_averages(ref.mean, dr.mean);
        compare_averages(ref.detected, dr.detected);
        compare_best(ref, dr);
    }

    // Comparing to all of the other matrix representations.
    {
        auto dc = scran_markers::score_markers_best<double>(*dense_column, groupings.data(), top, opt);
        compare_averages(ref.mean, dc.mean);
        compare_averages(ref.detected, dc.detected);
        compare_best(ref, dc);

        auto sr = scran_markers::score_markers_best<double>(*sparse_row, groupings.data(), top, opt);
        compare_averages(ref.mean, sr.mean);
        compare_averages(ref.detected, sr.detected);
        compare_best(ref, sr);

        auto sc = scran_markers::score_markers_best<double>(*sparse_column, groupings.data(), top, opt);
        compare_averages(ref.mean, sc.mean);
        compare_averages(ref.detected, sc.detected);
        compare_best(ref, sc);
    }

    // Comparing to quantile best; should be the same for 1 block.
    {
        auto qopt = opt;
        qopt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;

        auto qdr = scran_markers::score_markers_best<double>(*dense_row, groupings.data(), top, qopt);
        compare_averages(ref.mean, qdr.mean);
        compare_averages(ref.detected, qdr.detected);
        compare_best(ref, qdr);

        auto qdc = scran_markers::score_markers_best<double>(*dense_column, groupings.data(), top, qopt);
        compare_averages(ref.mean, qdc.mean);
        compare_averages(ref.detected, qdc.detected);
        compare_best(ref, qdc);

        auto qsr = scran_markers::score_markers_best<double>(*sparse_row, groupings.data(), top, qopt);
        compare_averages(ref.mean, qsr.mean);
        compare_averages(ref.detected, qsr.detected);
        compare_best(ref, qsr);

        auto qsc = scran_markers::score_markers_best<double>(*sparse_column, groupings.data(), top, qopt);
        compare_averages(ref.mean, qsc.mean);
        compare_averages(ref.detected, qsc.detected);
        compare_best(ref, qsc);
    }
}

INSTANTIATE_TEST_SUITE_P(
    ScoreMarkersBest,
    ScoreMarkersBestTest,
    ::testing::Combine(
        ::testing::Values(2, 5), // number of clusters
        ::testing::Values(false, true), // with or without the AUC?
        ::testing::Values(10, 50), // number of top markers
        ::testing::Values(true, false), // use larger or smaller effects
        ::testing::Values(true, false), // use bounds or not
        ::testing::Values(true, false), // keep ties or not
        ::testing::Values(1, 3) // number of threads
    )
);

/*********************************************/

class ScoreMarkersBestBlockedTest : public ScoreMarkersBestTestCore, public ::testing::TestWithParam<std::tuple<int, bool, scran_blocks::WeightPolicy, int, bool, bool, bool> > {
protected:
    inline static std::shared_ptr<tatami::Matrix<double, int> > dense_row, dense_column, sparse_row, sparse_column;

    static void SetUpTestSuite() {
        size_t nr = 198, nc = 201;
        dense_row.reset(
            new tatami::DenseRowMatrix<double, int>(
                nr,
                nc,
                scran_tests::simulate_vector(
                    nr * nc, 
                    []{
                        scran_tests::SimulateVectorParameters sparam;
                        sparam.density = 0.3;
                        sparam.seed = 4200;
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

TEST_P(ScoreMarkersBestBlockedTest, AgainstPairwiseMean) {
    auto param = GetParam();
    auto ngroups = std::get<0>(param);
    bool do_auc = std::get<1>(param);
    auto policy = std::get<2>(param);
    int top = std::get<3>(param);
    bool larger = std::get<4>(param);
    bool use_bound = std::get<5>(param);
    bool keep_ties = std::get<6>(param);

    auto NC = dense_row->ncol();
    std::vector<int> groupings = create_groupings(NC, ngroups);
    std::vector<int> blocks = create_blocks(NC, 3);
    auto ngenes = dense_row->nrow();

    scran_markers::ScoreMarkersBestOptions opt;
    opt.largest_cohens_d = larger;
    opt.largest_delta_detected = larger;
    opt.largest_delta_mean = larger;
    opt.largest_auc = larger;

    if (!use_bound) {
        opt.threshold_cohens_d.reset();
        opt.threshold_delta_detected.reset();
        opt.threshold_delta_mean.reset();
        opt.threshold_auc.reset();
    }

    opt.keep_ties = keep_ties;
    opt.compute_auc = do_auc; // false, if we want to check the running implementations.
    opt.block_weight_policy = policy;
    auto ref = scran_markers::score_markers_best_blocked<double>(*dense_row, groupings.data(), blocks.data(), top, opt);

    // Comparing score_markers_best_blocked against score_markers_pairwise_blocked + topicks::pick_top_genes.
    // The latter is less mind-bending but requires holding a large 3D matrix in memory.
    {
        scran_markers::ScoreMarkersPairwiseOptions popt;
        popt.compute_auc = do_auc;
        popt.block_weight_policy = policy;
        auto pairres = scran_markers::score_markers_pairwise_blocked(*dense_row, groupings.data(), blocks.data(), popt);
        compare_averages(ref.mean, pairres.mean);
        compare_averages(ref.detected, pairres.detected);

        compare_best_to_pairwise(ref.cohens_d, pairres.cohens_d, ngenes, ngroups, top, opt.largest_cohens_d, keep_ties, opt.threshold_cohens_d);
        compare_best_to_pairwise(ref.delta_mean, pairres.delta_mean, ngenes, ngroups, top, opt.largest_delta_mean, keep_ties, opt.threshold_delta_mean);
        compare_best_to_pairwise(ref.delta_detected, pairres.delta_detected, ngenes, ngroups, top, opt.largest_delta_detected, keep_ties, opt.threshold_delta_detected);
        if (do_auc) {
            compare_best_to_pairwise(ref.auc, pairres.auc, ngenes, ngroups, top, opt.largest_auc, keep_ties, opt.threshold_auc);
        }
    }

    // Note: skipping the multi-threaded checks as we've already got too many test cases.
    // Besides, the multi-threaded code is the same as the unblocked case.

    // Comparing to all of the other matrix representations.
    {
        auto dc = scran_markers::score_markers_best_blocked<double>(*dense_column, groupings.data(), blocks.data(), top, opt);
        compare_averages(ref.mean, dc.mean);
        compare_averages(ref.detected, dc.detected);
        compare_best(ref, dc);

        auto sr = scran_markers::score_markers_best_blocked<double>(*sparse_row, groupings.data(), blocks.data(), top, opt);
        compare_averages(ref.mean, sr.mean);
        compare_averages(ref.detected, sr.detected);
        compare_best(ref, sr);

        auto sc = scran_markers::score_markers_best_blocked<double>(*sparse_column, groupings.data(), blocks.data(), top, opt);
        compare_averages(ref.mean, sc.mean);
        compare_averages(ref.detected, sc.detected);
        compare_best(ref, sc);
    }
}

TEST_P(ScoreMarkersBestBlockedTest, AgainstPairwiseQuantile) {
    auto param = GetParam();
    auto ngroups = std::get<0>(param);
    bool do_auc = std::get<1>(param);
    auto policy = std::get<2>(param);
    int top = std::get<3>(param);
    bool larger = std::get<4>(param);
    bool use_bound = std::get<5>(param);
    bool keep_ties = std::get<6>(param);

    // Block weighting has no effect here, so we'll just short-circuit.
    if (policy == scran_blocks::WeightPolicy::EQUAL) {
        return;
    }

    auto NC = dense_row->ncol();
    std::vector<int> groupings = create_groupings(NC, ngroups);
    std::vector<int> blocks = create_blocks(NC, 3);
    auto ngenes = dense_row->nrow();

    scran_markers::ScoreMarkersBestOptions opt;
    opt.largest_cohens_d = larger;
    opt.largest_delta_detected = larger;
    opt.largest_delta_mean = larger;
    opt.largest_auc = larger;

    if (!use_bound) {
        opt.threshold_cohens_d.reset();
        opt.threshold_delta_detected.reset();
        opt.threshold_delta_mean.reset();
        opt.threshold_auc.reset();
    }

    opt.keep_ties = keep_ties;
    opt.compute_auc = do_auc; // false, if we want to check the running implementations.
    opt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
    auto ref = scran_markers::score_markers_best_blocked<double>(*dense_row, groupings.data(), blocks.data(), top, opt);

    // Comparing score_markers_best_blocked against score_markers_pairwise_blocked + topicks::pick_top_genes.
    // The latter is less mind-bending but requires holding a large 3D matrix in memory.
    {
        scran_markers::ScoreMarkersPairwiseOptions popt;
        popt.compute_auc = do_auc;
        popt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
        auto pairres = scran_markers::score_markers_pairwise_blocked(*dense_row, groupings.data(), blocks.data(), popt);
        compare_averages(ref.mean, pairres.mean);
        compare_averages(ref.detected, pairres.detected);

        compare_best_to_pairwise(ref.cohens_d, pairres.cohens_d, ngenes, ngroups, top, opt.largest_cohens_d, keep_ties, opt.threshold_cohens_d);
        compare_best_to_pairwise(ref.delta_mean, pairres.delta_mean, ngenes, ngroups, top, opt.largest_delta_mean, keep_ties, opt.threshold_delta_mean);
        compare_best_to_pairwise(ref.delta_detected, pairres.delta_detected, ngenes, ngroups, top, opt.largest_delta_detected, keep_ties, opt.threshold_delta_detected);
        if (do_auc) {
            compare_best_to_pairwise(ref.auc, pairres.auc, ngenes, ngroups, top, opt.largest_auc, keep_ties, opt.threshold_auc);
        }
    }

    // Note: skipping the multi-threaded checks as we've already got too many test cases.
    // Besides, the multi-threaded code is the same as the unblocked case.

    // Comparing to all of the other matrix representations.
    {
        auto dc = scran_markers::score_markers_best_blocked<double>(*dense_column, groupings.data(), blocks.data(), top, opt);
        compare_averages(ref.mean, dc.mean);
        compare_averages(ref.detected, dc.detected);
        compare_best(ref, dc);

        auto sr = scran_markers::score_markers_best_blocked<double>(*sparse_row, groupings.data(), blocks.data(), top, opt);
        compare_averages(ref.mean, sr.mean);
        compare_averages(ref.detected, sr.detected);
        compare_best(ref, sr);

        auto sc = scran_markers::score_markers_best_blocked<double>(*sparse_column, groupings.data(), blocks.data(), top, opt);
        compare_averages(ref.mean, sc.mean);
        compare_averages(ref.detected, sc.detected);
        compare_best(ref, sc);
    }
}

INSTANTIATE_TEST_SUITE_P(
    ScoreMarkersBest,
    ScoreMarkersBestBlockedTest,
    ::testing::Combine(
        ::testing::Values(2, 5), // number of clusters
        ::testing::Values(false, true), // with or without the AUC?
        ::testing::Values(scran_blocks::WeightPolicy::NONE, scran_blocks::WeightPolicy::EQUAL), // block weighting method.
        ::testing::Values(10, 50), // number of top markers
        ::testing::Values(true, false), // use larger or smaller effects
        ::testing::Values(true, false), // use bounds or not
        ::testing::Values(true, false)  // keep ties or not
    )
);

/*********************************************/

class ScoreMarkersBestScenariosTest : public ScoreMarkersBestTestCore, public ::testing::Test {};

TEST_F(ScoreMarkersBestScenariosTest, Thresholds) {
    int nrows = 291, ncols = 91;
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

    int ngroups = 3;
    std::vector<int> groupings = create_groupings(mat.ncol(), ngroups);

    int top = 10;
    scran_markers::ScoreMarkersBestOptions sopt;
    auto ref = scran_markers::score_markers_best<double>(mat, groupings.data(), top, sopt);
    sopt.threshold = 1;
    auto out = scran_markers::score_markers_best<double>(mat, groupings.data(), top, sopt);

    int some_diff = 0;
    for (int l = 0; l < ngroups; ++l) {
        // Not affected.
        EXPECT_EQ(ref.mean[l], out.mean[l]);
        EXPECT_EQ(ref.detected[l], out.detected[l]);

        for (int k = 0; k < ngroups; ++k) {
            if (k == l) {
                continue;
            }

            EXPECT_EQ(ref.delta_mean[l][k], out.delta_mean[l][k]);
            EXPECT_EQ(ref.delta_detected[l][k], out.delta_detected[l][k]);

            const auto& rcohen = ref.cohens_d[l][k];
            const auto& ocohen = out.cohens_d[l][k];
            const std::size_t ncohen = std::min(rcohen.size(), ocohen.size());
            for (std::size_t i = 0; i < ncohen; ++i) {
                EXPECT_GT(rcohen[i].second, ocohen[i].second);
            }

            const auto& rauc = ref.auc[l][k];
            const auto& oauc = out.auc[l][k];
            const std::size_t nauc = std::min(rauc.size(), oauc.size());
            for (std::size_t i = 0; i < nauc; ++i) {
                // Can't use GT as the ranks might not be affected by the threshold.
                EXPECT_GE(rauc[i].second, oauc[i].second);
                some_diff += (rauc[i].second != oauc[i].second);
            }
        }
    }

    EXPECT_GT(some_diff, 0); // though hopefully at least one is changed.

    // Quantile should give the same results for a single block.
    auto qopt = sopt;
    qopt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
    auto qout = scran_markers::score_markers_best<double>(mat, groupings.data(), top, qopt);
    compare_averages(out.mean, qout.mean);
    compare_averages(out.detected, qout.detected);
    compare_best(out, qout);
}

TEST_F(ScoreMarkersBestScenariosTest, Missing) {
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

    int top = 10;
    scran_markers::ScoreMarkersBestOptions opt;
    auto ref = scran_markers::score_markers_best<double>(mat, groupings.data(), top, opt);

    // Zero is effectively the missing group here.
    for (auto& g : groupings) {
        ++g;
    }
    auto lost = scran_markers::score_markers_best<double>(mat, groupings.data(), top, opt);

    // Everything should be empty.
    for (int g = 1; g <= ngroups; ++g) {
        EXPECT_TRUE(lost.cohens_d[0][g].empty());
        EXPECT_TRUE(lost.cohens_d[g][0].empty());
        EXPECT_TRUE(lost.auc[0][g].empty());
        EXPECT_TRUE(lost.auc[g][0].empty());
        EXPECT_TRUE(lost.delta_mean[0][g].empty());
        EXPECT_TRUE(lost.delta_mean[g][0].empty());
        EXPECT_TRUE(lost.delta_detected[0][g].empty());
        EXPECT_TRUE(lost.delta_detected[g][0].empty());
    }
    for (int r = 0; r < nrows; ++r) {
        EXPECT_TRUE(std::isnan(lost.mean[0][r]));
        EXPECT_TRUE(std::isnan(lost.detected[0][r]));
    }

    // Other metrics should be the same as usual.
    for (int g1 = 0; g1 < ngroups; ++g1) {
        EXPECT_EQ(lost.mean[g1+1], ref.mean[g1]);
        EXPECT_EQ(lost.detected[g1+1], ref.detected[g1]);

        for (int g2 = 0; g2 < ngroups; ++g2) {
            EXPECT_EQ(lost.cohens_d[g1 + 1][g2 + 1], ref.cohens_d[g1][g2]);
            EXPECT_EQ(lost.auc[g1 + 1][g2 + 1], ref.auc[g1][g2]);
            EXPECT_EQ(lost.delta_mean[g1 + 1][g2 + 1], ref.delta_mean[g1][g2]);
            EXPECT_EQ(lost.delta_detected[g1 + 1][g2 + 1], ref.delta_detected[g1][g2]);
        }
    }

    // Quantile should give the same results for a single block.
    auto qopt = opt;
    qopt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
    auto qlost = scran_markers::score_markers_best<double>(mat, groupings.data(), top, qopt);
    compare_averages(lost.mean, qlost.mean);
    compare_averages(lost.detected, qlost.detected);
    compare_best(lost, qlost);
}

TEST_F(ScoreMarkersBestScenariosTest, BlockConfounded) {
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

    int top = 10;
    scran_markers::ScoreMarkersBestOptions opt;
    auto comres = scran_markers::score_markers_best_blocked<double>(*mat, groupings.data(), blocks.data(), top, opt);

    // First group should only be NaN's.
    for (int g = 1; g < ngroups; ++g) {
        EXPECT_TRUE(comres.cohens_d[0][g].empty());
        EXPECT_TRUE(comres.cohens_d[g][0].empty());
        EXPECT_TRUE(comres.auc[0][g].empty());
        EXPECT_TRUE(comres.auc[g][0].empty());
        EXPECT_TRUE(comres.delta_mean[0][g].empty());
        EXPECT_TRUE(comres.delta_mean[g][0].empty());
        EXPECT_TRUE(comres.delta_detected[0][g].empty());
        EXPECT_TRUE(comres.delta_detected[g][0].empty());
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
    auto ref = scran_markers::score_markers_best<double>(*sub, subgroups.data(), top, opt);

    for (int g1 = 1; g1 < ngroups; ++g1) {
        EXPECT_EQ(comres.mean[g1], ref.mean[g1 - 1]);
        EXPECT_EQ(comres.detected[g1], ref.detected[g1 - 1]);

        for (int g2 = 1; g2 < ngroups; ++g2) {
            EXPECT_EQ(comres.cohens_d[g1][g2], ref.cohens_d[g1 - 1][g2 - 1]);
            EXPECT_EQ(comres.auc[g1][g2], ref.auc[g1 - 1][g2 - 1]);
            EXPECT_EQ(comres.delta_mean[g1][g2], ref.delta_mean[g1 - 1][g2 - 1]);
            EXPECT_EQ(comres.delta_detected[g1][g2], ref.delta_detected[g1 - 1][g2 - 1]);
        }
    }

    // Quantile should give the same results, as there's basically only one block;
    // the second block is fully confounded.
    auto qopt = opt;
    qopt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
    auto qcomres = scran_markers::score_markers_best_blocked<double>(*mat, groupings.data(), blocks.data(), top, qopt);
    compare_averages(comres.mean, qcomres.mean);
    compare_averages(comres.detected, qcomres.detected);
    compare_best(comres, qcomres);
}

TEST_F(ScoreMarkersBestScenariosTest, Empty) {
    int nrows = 0, ncols = 66;
    tatami::DenseMatrix<double, int, std::vector<double> > mat(nrows, ncols, std::vector<double>(), true);

    int ngroups = 4;
    std::vector<int> groupings = create_groupings(ncols, ngroups);

    scran_markers::ScoreMarkersBestOptions opts;
    int top = 10;
    auto out = scran_markers::score_markers_best<double>(mat, groupings.data(), top, opts);

    for (int g = 0; g < ngroups; ++g) {
        EXPECT_TRUE(out.mean[g].empty());
        EXPECT_TRUE(out.detected[g].empty());

        EXPECT_EQ(out.cohens_d[g].size(), ngroups);
        EXPECT_EQ(out.auc[g].size(), ngroups);
        EXPECT_EQ(out.delta_mean[g].size(), ngroups);
        EXPECT_EQ(out.delta_detected[g].size(), ngroups);

        for (int g2 = 0; g2 < ngroups; ++g2) {
            EXPECT_TRUE(out.cohens_d[g][g2].empty());
            EXPECT_TRUE(out.auc[g][g2].empty());
            EXPECT_TRUE(out.delta_mean[g][g2].empty());
            EXPECT_TRUE(out.delta_detected[g][g2].empty());
        }
    }
}

/*********************************************/

class ScoreMarkersBestOneAtATimeTest : public ScoreMarkersBestTestCore, public ::testing::TestWithParam<int> {
protected:
    inline static std::shared_ptr<tatami::Matrix<double, int> > dense_row, dense_column, sparse_row, sparse_column;

    static void SetUpTestSuite() {
        size_t nr = 128, nc = 302;
        dense_row.reset(
            new tatami::DenseRowMatrix<double, int>(
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
            )
        );

        dense_column = tatami::convert_to_dense(dense_row.get(), false);
        sparse_row = tatami::convert_to_compressed_sparse(dense_row.get(), true);
        sparse_column = tatami::convert_to_compressed_sparse(dense_row.get(), false);
    }
};

TEST_P(ScoreMarkersBestOneAtATimeTest, Basic) {
    auto NC = dense_row->ncol();
    std::vector<int> groupings = create_groupings(NC, 3);

    const tatami::Matrix<double, int>* mat;
    switch (GetParam()) {
        case 0:
            mat = dense_row.get(); break;
        case 1:
            mat = dense_column.get(); break;
        case 2:
            mat = sparse_row.get(); break;
        case 3:
            mat = sparse_column.get(); break;
    }

    scran_markers::ScoreMarkersBestOptions opt;
    int top = 15;
    auto ref = scran_markers::score_markers_best<double>(*mat, groupings.data(), top, opt);

    // Only the group mean.
    {
        scran_markers::ScoreMarkersBestOptions opt;
        opt.compute_group_detected = false;
        opt.compute_cohens_d = false;
        opt.compute_auc = false;
        opt.compute_delta_mean = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_best<double>(*mat, groupings.data(), top, opt);
        compare_averages(ref.mean, alt.mean);
        EXPECT_TRUE(alt.detected.empty());
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only the group detected proportions.
    {
        scran_markers::ScoreMarkersBestOptions opt;
        opt.compute_group_mean = false;
        opt.compute_cohens_d = false;
        opt.compute_auc = false;
        opt.compute_delta_mean = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_best<double>(*mat, groupings.data(), top, opt);
        compare_averages(ref.detected, alt.detected);
        EXPECT_TRUE(alt.mean.empty());
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only Cohen's d.
    {
        scran_markers::ScoreMarkersBestOptions opt;
        opt.compute_group_mean = false;
        opt.compute_group_detected = false;
        opt.compute_auc = false;
        opt.compute_delta_mean = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_best<double>(*mat, groupings.data(), top, opt);
        compare_best(alt.cohens_d, ref.cohens_d);
        EXPECT_TRUE(alt.mean.empty());
        EXPECT_TRUE(alt.detected.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only AUC.
    {
        scran_markers::ScoreMarkersBestOptions opt;
        opt.compute_group_mean = false;
        opt.compute_group_detected = false;
        opt.compute_cohens_d = false;
        opt.compute_delta_mean = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_best<double>(*mat, groupings.data(), top, opt);
        compare_best(alt.auc, ref.auc);
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only delta-mean.
    {
        scran_markers::ScoreMarkersBestOptions opt;
        opt.compute_group_mean = false;
        opt.compute_group_detected = false;
        opt.compute_cohens_d = false;
        opt.compute_auc = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_best<double>(*mat, groupings.data(), top, opt);
        compare_best(alt.delta_mean, ref.delta_mean);
        EXPECT_TRUE(alt.mean.empty());
        EXPECT_TRUE(alt.detected.empty());
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only delta-mean.
    {
        scran_markers::ScoreMarkersBestOptions opt;
        opt.compute_group_mean = false;
        opt.compute_group_detected = false;
        opt.compute_cohens_d = false;
        opt.compute_auc = false;
        opt.compute_delta_mean = false;

        auto alt = scran_markers::score_markers_best<double>(*mat, groupings.data(), top, opt);
        compare_best(alt.delta_detected, ref.delta_detected);
        EXPECT_TRUE(alt.mean.empty());
        EXPECT_TRUE(alt.detected.empty());
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
    }
}

INSTANTIATE_TEST_SUITE_P(
    ScoreMarkersBest,
    ScoreMarkersBestOneAtATimeTest,
    ::testing::Values(0, 1, 2, 3)
);
