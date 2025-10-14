#include "scran_tests/scran_tests.hpp"
#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"

#include "scran_markers/score_markers_summary.hpp"
#include "scran_markers/score_markers_pairwise.hpp"
#include "scran_markers/summarize_effects.hpp"

#include "utils.h"

class ScoreMarkersSummaryTestCore {
protected:
    static void compare_averages(const std::vector<std::vector<double> >& res, const std::vector<std::vector<double> >& other) {
        const int ngroups = res.size();
        ASSERT_EQ(ngroups, other.size());
        for (int l = 0; l < ngroups; ++l) {
            scran_tests::compare_almost_equal(res[l], other[l]);
        }
    }

    template<class Summaries_>
    static void compare_summaries_for_effect(int ngroups, const std::vector<Summaries_>& res, const std::vector<Summaries_>& other) {
        EXPECT_EQ(res.size(), ngroups);
        EXPECT_EQ(other.size(), ngroups);

        for (int l = 0; l < ngroups; ++l) {
            scran_tests::compare_almost_equal(res[l].min, other[l].min);
            scran_tests::compare_almost_equal(res[l].mean, other[l].mean);
            scran_tests::compare_almost_equal(res[l].median, other[l].median);
            scran_tests::compare_almost_equal(res[l].max, other[l].max);

            // Don't compare min-rank values here, as minor numerical differences
            // can change the ranks by a small amount when effects are tied.
            EXPECT_EQ(res[l].min_rank.size(), other[l].min_rank.size());
        }
    }

    template<class Results_>
    static void compare_effects(int ngroups, const Results_& res, const Results_& other, bool do_auc) {
        compare_summaries_for_effect(ngroups, res.cohens_d, other.cohens_d);
        compare_summaries_for_effect(ngroups, res.delta_mean, other.delta_mean);
        compare_summaries_for_effect(ngroups, res.delta_detected, other.delta_detected);
        if (do_auc) {
            compare_summaries_for_effect(ngroups, res.auc, other.auc);
        }
    }

    template<class Summaries_>
    static void check_effects(size_t ngenes, size_t group, const std::vector<Summaries_>& effects) {
        for (size_t g = 0; g < ngenes; ++g) {
            double curmin = effects[group].min[g];
            double curmean = effects[group].mean[g];
            double curmed = effects[group].median[g];
            double curmax = effects[group].max[g];
            double currank = effects[group].min_rank[g];

            EXPECT_LE(curmin, curmean);
            EXPECT_LE(curmin, curmed);
            EXPECT_GE(curmax, curmean);
            EXPECT_GE(curmax, curmed);

            if (curmed > curmin && std::isfinite(curmin)) { // i.e., not -Inf.
                EXPECT_GT(curmean, curmin);
            }
            if (curmed < curmax && std::isfinite(curmax)) { // i.e., not Inf.
                EXPECT_LT(curmean, curmax);
            }

            EXPECT_GE(currank, 1);
            EXPECT_LE(currank, ngenes);
        }
    }

    static void wipe_min_rank(scran_markers::ScoreMarkersSummaryResults<double, int>& results) { 
        for (auto& x : results.auc) {
            x.min_rank.clear();
        }
        for (auto& x : results.cohens_d) {
            x.min_rank.clear();
        }
        for (auto& x : results.delta_detected) {
            x.min_rank.clear();
        }
        for (auto& x : results.delta_mean) {
            x.min_rank.clear();
        }
    }
};

/*********************************************/

class ScoreMarkersSummaryTest : public ScoreMarkersSummaryTestCore, public ::testing::TestWithParam<std::tuple<int, bool, int> > {
protected:
    inline static std::shared_ptr<tatami::Matrix<double, int> > dense_row, dense_column, sparse_row, sparse_column;

    static void SetUpTestSuite() {
        size_t nr = 548, nc = 142;
        dense_row.reset(
            new tatami::DenseRowMatrix<double, int>(
                nr,
                nc,
                scran_tests::simulate_vector(
                    nr * nc, 
                    []{
                        scran_tests::SimulationParameters sparam;
                        sparam.density = 0.1;
                        sparam.seed = 4242;
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

TEST_P(ScoreMarkersSummaryTest, Basic) {
    auto param = GetParam();
    auto ngroups = std::get<0>(param);
    bool do_auc = std::get<1>(param);
    auto nthreads = std::get<2>(param);

    std::vector<int> groupings = create_groupings(dense_row->ncol(), ngroups);
    size_t ngenes = dense_row->nrow();

    scran_markers::ScoreMarkersSummaryOptions opt;
    opt.compute_auc = do_auc; // false, if we want to check the running implementations.
    auto ref = scran_markers::score_markers_summary(*dense_row, groupings.data(), opt);

    if (nthreads == 1) {
        // Running some further checks on the effects.
        for (int l = 0; l < ngroups; ++l) {
            check_effects(ngenes, l, ref.cohens_d);
            check_effects(ngenes, l, ref.delta_mean);

            check_effects(ngenes, l, ref.delta_detected);
            for (size_t g = 0; g < ngenes; ++g) {
                double curmin = ref.delta_detected[l].min[g];
                double curmax = ref.delta_detected[l].max[g];
                EXPECT_GE(curmin, -1);
                EXPECT_LE(curmax, 1);
            }

            if (do_auc) {
                check_effects(ngenes, l, ref.auc);
                for (size_t g = 0; g < ngenes; ++g) {
                    double curmin = ref.auc[l].min[g];
                    double curmax = ref.auc[l].max[g];
                    EXPECT_GE(curmin, 0);
                    EXPECT_LE(curmax, 1);
                }
            }
        }

        // Comparing the efficient ScoreMarkers implementation against the
        // PairwiseEffects + SummarizeEffects combination, which is less mind-bending
        // but requires holding a large 3D matrix in memory.
        scran_markers::ScoreMarkersPairwiseOptions popt;
        popt.compute_auc = do_auc;
        auto pairres = scran_markers::score_markers_pairwise(*dense_row, groupings.data(), popt);
        compare_averages(ref.mean, pairres.mean);
        compare_averages(ref.detected, pairres.detected);

        scran_markers::SummarizeEffectsOptions sopt;
        auto cohen_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.cohens_d.data(), sopt);
        compare_summaries_for_effect(ngroups, cohen_summ, ref.cohens_d);

        auto dm_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.delta_mean.data(), sopt);
        compare_summaries_for_effect(ngroups, dm_summ, ref.delta_mean);

        auto dd_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.delta_detected.data(), sopt);
        compare_summaries_for_effect(ngroups, dd_summ, ref.delta_detected);

        if (do_auc) {
            auto auc_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.auc.data(), sopt);
            compare_summaries_for_effect(ngroups, auc_summ, ref.auc);
        }

    } else {
        opt.num_threads = nthreads;
        auto dr = scran_markers::score_markers_summary(*dense_row, groupings.data(), opt);
        compare_averages(ref.mean, dr.mean);
        compare_averages(ref.detected, dr.detected);
        compare_effects(ngroups, ref, dr, do_auc);
    }

    // Comparing to all of the other matrix representations.
    {
        auto dc = scran_markers::score_markers_summary(*dense_column, groupings.data(), opt);
        compare_averages(ref.mean, dc.mean);
        compare_averages(ref.detected, dc.detected);
        compare_effects(ngroups, ref, dc, do_auc);

        auto sr = scran_markers::score_markers_summary(*sparse_row, groupings.data(), opt);
        compare_averages(ref.mean, sr.mean);
        compare_averages(ref.detected, sr.detected);
        compare_effects(ngroups, ref, sr, do_auc);

        auto sc = scran_markers::score_markers_summary(*sparse_column, groupings.data(), opt);
        compare_averages(ref.mean, sc.mean);
        compare_averages(ref.detected, sc.detected);
        compare_effects(ngroups, ref, sc, do_auc);
    }

    // Comparing to what happens if we disable minrank, which allows for more efficient AUC calculations.
    if (do_auc) {
        scran_markers::ScoreMarkersSummaryOptions opt_no_mr = opt;
        opt_no_mr.compute_min_rank = false;
        wipe_min_rank(ref);

        auto dr_nmr = scran_markers::score_markers_summary(*dense_row, groupings.data(), opt_no_mr);
        compare_effects(ngroups, ref, dr_nmr, true);

        auto dc_nmr = scran_markers::score_markers_summary(*dense_column, groupings.data(), opt_no_mr);
        compare_effects(ngroups, ref, dc_nmr, true);

        auto sr_nmr = scran_markers::score_markers_summary(*sparse_row, groupings.data(), opt_no_mr);
        compare_effects(ngroups, ref, sr_nmr, true);

        auto sc_nmr = scran_markers::score_markers_summary(*sparse_column, groupings.data(), opt_no_mr);
        compare_effects(ngroups, ref, sc_nmr, true);
    }
}

INSTANTIATE_TEST_SUITE_P(
    ScoreMarkers,
    ScoreMarkersSummaryTest,
    ::testing::Combine(
        ::testing::Values(2, 3, 4, 5), // number of clusters
        ::testing::Values(false, true), // with or without the AUC?
        ::testing::Values(1, 3) // number of threads
    )
);

/*********************************************/

class ScoreMarkersSummaryBlockedTest : public ScoreMarkersSummaryTestCore, public ::testing::TestWithParam<std::tuple<int, bool, scran_blocks::WeightPolicy, int> > {
protected:
    inline static std::shared_ptr<tatami::Matrix<double, int> > dense_row, dense_column, sparse_row, sparse_column;

    static void SetUpTestSuite() {
        size_t nr = 248, nc = 287;
        dense_row.reset(
            new tatami::DenseRowMatrix<double, int>(
                nr,
                nc,
                scran_tests::simulate_vector(
                    nr * nc, 
                    []{
                        scran_tests::SimulationParameters sparam;
                        sparam.density = 0.1;
                        sparam.seed = 666 ;
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

TEST_P(ScoreMarkersSummaryBlockedTest, AgainstPairwise) {
    auto param = GetParam();
    auto ngroups = std::get<0>(param);
    bool do_auc = std::get<1>(param);
    auto policy = std::get<2>(param);
    auto nthreads = std::get<3>(param);

    auto NC = dense_row->ncol();
    std::vector<int> groupings = create_groupings(NC, ngroups);
    std::vector<int> blocks = create_blocks(NC, 3);

    scran_markers::ScoreMarkersSummaryOptions opt;
    opt.compute_auc = do_auc; // false, if we want to check the running implementations.
    opt.block_weight_policy = policy;
    auto ref = scran_markers::score_markers_summary_blocked(*dense_row, groupings.data(), blocks.data(), opt);

    if (nthreads == 1) {
        scran_markers::ScoreMarkersPairwiseOptions popt;
        popt.compute_auc = do_auc;
        popt.block_weight_policy = policy;
        auto pairres = scran_markers::score_markers_pairwise_blocked(*dense_row, groupings.data(), blocks.data(), popt);
        compare_averages(ref.mean, pairres.mean);
        compare_averages(ref.detected, pairres.detected);

        scran_markers::SummarizeEffectsOptions sopt;
        size_t ngenes = dense_row->nrow();
        auto cohen_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.cohens_d.data(), sopt);
        compare_summaries_for_effect(ngroups, cohen_summ, ref.cohens_d);

        auto dm_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.delta_mean.data(), sopt);
        compare_summaries_for_effect(ngroups, dm_summ, ref.delta_mean);

        auto dd_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.delta_detected.data(), sopt);
        compare_summaries_for_effect(ngroups, dd_summ, ref.delta_detected);

        if (do_auc) {
            auto auc_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.auc.data(), sopt);
            compare_summaries_for_effect(ngroups, auc_summ, ref.auc);
        }

    } else {
        opt.num_threads = nthreads;
        auto dr = scran_markers::score_markers_summary_blocked(*dense_row, groupings.data(), blocks.data(), opt);
        compare_averages(ref.mean, dr.mean);
        compare_averages(ref.detected, dr.detected);
        compare_effects(ngroups, ref, dr, do_auc);
    }

    // Checking the other references.
    {
        auto dc = scran_markers::score_markers_summary_blocked(*dense_column, groupings.data(), blocks.data(), opt);
        compare_averages(ref.mean, dc.mean);
        compare_averages(ref.detected, dc.detected);
        compare_effects(ngroups, ref, dc, do_auc);

        auto sr = scran_markers::score_markers_summary_blocked(*sparse_row, groupings.data(), blocks.data(), opt);
        compare_averages(ref.mean, sr.mean);
        compare_averages(ref.detected, sr.detected);
        compare_effects(ngroups, ref, sr, do_auc);

        auto sc = scran_markers::score_markers_summary_blocked(*sparse_column, groupings.data(), blocks.data(), opt);
        compare_averages(ref.mean, sc.mean);
        compare_averages(ref.detected, sc.detected);
        compare_effects(ngroups, ref, sc, do_auc);
    }

    // Comparing to what happens if we disable minrank, which allows for more efficient AUC calculations.
    if (do_auc) {
        scran_markers::ScoreMarkersSummaryOptions opt_no_mr = opt;
        opt_no_mr.compute_min_rank = false;
        wipe_min_rank(ref);

        auto dr_nmr = scran_markers::score_markers_summary_blocked(*dense_row, groupings.data(), blocks.data(), opt_no_mr);
        compare_effects(ngroups, ref, dr_nmr, true);

        auto dc_nmr = scran_markers::score_markers_summary_blocked(*dense_column, groupings.data(), blocks.data(), opt_no_mr);
        compare_effects(ngroups, ref, dc_nmr, true);

        auto sr_nmr = scran_markers::score_markers_summary_blocked(*sparse_row, groupings.data(), blocks.data(), opt_no_mr);
        compare_effects(ngroups, ref, sr_nmr, true);

        auto sc_nmr = scran_markers::score_markers_summary_blocked(*sparse_column, groupings.data(), blocks.data(), opt_no_mr);
        compare_effects(ngroups, ref, sc_nmr, true);
    }
}

INSTANTIATE_TEST_SUITE_P(
    ScoreMarkersSummary,
    ScoreMarkersSummaryBlockedTest,
    ::testing::Combine(
        ::testing::Values(2, 3, 4, 5), // number of clusters
        ::testing::Values(false, true), // with or without the AUC?
        ::testing::Values(scran_blocks::WeightPolicy::NONE, scran_blocks::WeightPolicy::EQUAL), // block weighting method.
        ::testing::Values(1, 3) // number of threads
    )
);

/*********************************************/

// Checking that we get the same results for different cache sizes.

class ScoreMarkersSummaryCacheTest : public ScoreMarkersSummaryTestCore, public ::testing::TestWithParam<std::tuple<int, int> > {
protected:
    inline static std::shared_ptr<tatami::Matrix<double, int> > matrix;

    static void SetUpTestSuite() {
        size_t nr = 366, nc = 179;
        matrix.reset(
            new tatami::DenseRowMatrix<double, int>(
                nr,
                nc,
                scran_tests::simulate_vector(
                    nr * nc, 
                    []{
                        scran_tests::SimulationParameters sparam;
                        sparam.density = 0.1;
                        sparam.seed = 420;
                        return sparam;
                    }()
                )
            )
        );
    }

    inline static scran_markers::ScoreMarkersSummaryResults<double, int> ref;
    inline static int last_ngroups = 0;
};

TEST_P(ScoreMarkersSummaryCacheTest, Basic) {
    auto param = GetParam();

    auto ngroups = std::get<0>(param);
    std::vector<int> groupings = create_groupings(matrix->ncol(), ngroups);
    auto cache_size = std::get<1>(param);

    if (ngroups != last_ngroups) { // avoid unnecessary recompute across the test suite.
        scran_markers::ScoreMarkersSummaryOptions opt;
        opt.cache_size = 0;
        ref = scran_markers::score_markers_summary(*matrix, groupings.data(), opt);
        last_ngroups = ngroups;
    }

    scran_markers::ScoreMarkersSummaryOptions opt;
    opt.cache_size = cache_size;
    auto cached = scran_markers::score_markers_summary(*matrix, groupings.data(), opt);

    compare_averages(ref.mean, cached.mean);
    compare_averages(ref.detected, cached.detected);
    compare_effects(ngroups, ref, cached, true);
}

INSTANTIATE_TEST_SUITE_P(
    ScoreMarkersSummary,
    ScoreMarkersSummaryCacheTest,
    ::testing::Combine(
        ::testing::Values(2, 4, 8, 16), // number of clusters
        ::testing::Values(5, 10, 20, 100) // size of the cache
    )
);

/*********************************************/

TEST(ScoreMarkersSummaryScenarios, Thresholds) {
    int nrows = 291, ncols = 91;
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
    std::vector<int> groupings = create_groupings(mat.ncol(), ngroups);

    scran_markers::ScoreMarkersSummaryOptions sopt;
    auto ref = scran_markers::score_markers_summary(mat, groupings.data(), sopt);
    sopt.threshold = 1;
    auto out = scran_markers::score_markers_summary(mat, groupings.data(), sopt);

    int some_diff = 0;
    for (int l = 0; l < ngroups; ++l) {
        for (int g = 0; g < nrows; ++g) {
            EXPECT_GT(ref.cohens_d[l].min[g], out.cohens_d[l].min[g]);
            EXPECT_GT(ref.cohens_d[l].mean[g], out.cohens_d[l].mean[g]);
            EXPECT_GT(ref.cohens_d[l].median[g], out.cohens_d[l].median[g]);
            EXPECT_GT(ref.cohens_d[l].max[g], out.cohens_d[l].max[g]);

            // '>' is not guaranteed due to imprecision with ranks... but (see below).
            some_diff += (ref.auc[l].mean[g] > out.auc[l].mean[g]);
            EXPECT_GE(ref.auc[l].min[g], out.auc[l].min[g]);
            EXPECT_GE(ref.auc[l].mean[g], out.auc[l].mean[g]);
            EXPECT_GE(ref.auc[l].median[g], out.auc[l].median[g]);
            EXPECT_GE(ref.auc[l].max[g], out.auc[l].max[g]);

            // These are not affected.
            EXPECT_EQ(ref.delta_mean[l].mean[g], out.delta_mean[l].mean[g]);
            EXPECT_EQ(ref.delta_detected[l].mean[g], out.delta_detected[l].mean[g]);
        }

        // Also not affected.
        EXPECT_EQ(ref.mean[l], out.mean[l]);
        EXPECT_EQ(ref.detected[l], out.detected[l]);
    }

    EXPECT_TRUE(some_diff); // (from above)... at least one is '>', hopefully.
}

TEST(ScoreMarkersSummaryScenarios, BlockConfounded) {
    int nrows = 103, ncols = 290;
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

    // Block is fully confounded with the first group.
    std::vector<int> blocks(ncols);
    for (int c = 0; c < ncols; ++c) {
        blocks[c] = groupings[c] == 0;
    }

    scran_markers::ScoreMarkersSummaryOptions sopt;
    auto comres = scran_markers::score_markers_summary_blocked(mat, groupings.data(), blocks.data(), sopt);

    // Excluding the group and running on the remaining samples.
    std::vector<int> subgroups;
    std::vector<int> keep;
    for (int c = 0; c < ncols; ++c) {
        auto g = groupings[c];
        if (g != 0) {
            subgroups.push_back(g - 1);
            keep.push_back(c);
        }
    }

    auto sub = tatami::make_DelayedSubset(tatami::wrap_shared_ptr(&mat), std::move(keep), false);
    auto ref = scran_markers::score_markers_summary(*sub, subgroups.data(), sopt); 

    // First group is effectively all-NA.
    for (int g = 0; g < nrows; ++g) {
        EXPECT_TRUE(std::isnan(comres.cohens_d[0].min[g]));
        EXPECT_TRUE(std::isnan(comres.cohens_d[0].max[g]));
        EXPECT_TRUE(std::isnan(comres.cohens_d[0].mean[g]));
        EXPECT_TRUE(std::isnan(comres.cohens_d[0].median[g]));
        EXPECT_EQ(comres.cohens_d[0].min_rank[g], nrows);

        EXPECT_TRUE(std::isnan(comres.auc[0].min[g]));
        EXPECT_TRUE(std::isnan(comres.auc[0].max[g]));
        EXPECT_TRUE(std::isnan(comres.auc[0].mean[g]));
        EXPECT_TRUE(std::isnan(comres.auc[0].median[g]));
        EXPECT_EQ(comres.auc[0].min_rank[g], nrows);

        EXPECT_TRUE(std::isnan(comres.delta_mean[0].min[g]));
        EXPECT_TRUE(std::isnan(comres.delta_mean[0].max[g]));
        EXPECT_TRUE(std::isnan(comres.delta_mean[0].mean[g]));
        EXPECT_TRUE(std::isnan(comres.delta_mean[0].median[g]));
        EXPECT_EQ(comres.delta_mean[0].min_rank[g], nrows);

        EXPECT_TRUE(std::isnan(comres.delta_detected[0].min[g]));
        EXPECT_TRUE(std::isnan(comres.delta_detected[0].max[g]));
        EXPECT_TRUE(std::isnan(comres.delta_detected[0].mean[g]));
        EXPECT_TRUE(std::isnan(comres.delta_detected[0].median[g]));
        EXPECT_EQ(comres.delta_detected[0].min_rank[g], nrows);
    }

    // Expect all but the first group to give the same results.
    for (int l = 1; l < ngroups; ++l) {
        scran_tests::compare_almost_equal(ref.cohens_d[l-1].min, comres.cohens_d[l].min);
        scran_tests::compare_almost_equal(ref.cohens_d[l-1].mean, comres.cohens_d[l].mean);
        scran_tests::compare_almost_equal(ref.cohens_d[l-1].median, comres.cohens_d[l].median);
        scran_tests::compare_almost_equal(ref.cohens_d[l-1].max, comres.cohens_d[l].max);
        scran_tests::compare_almost_equal(ref.cohens_d[l-1].min_rank, comres.cohens_d[l].min_rank);

        scran_tests::compare_almost_equal(ref.auc[l-1].min, comres.auc[l].min);
        scran_tests::compare_almost_equal(ref.auc[l-1].mean, comres.auc[l].mean);
        scran_tests::compare_almost_equal(ref.auc[l-1].median, comres.auc[l].median);
        scran_tests::compare_almost_equal(ref.auc[l-1].max, comres.auc[l].max);
        scran_tests::compare_almost_equal(ref.auc[l-1].min_rank, comres.auc[l].min_rank);

        scran_tests::compare_almost_equal(ref.delta_mean[l-1].min, comres.delta_mean[l].min);
        scran_tests::compare_almost_equal(ref.delta_mean[l-1].mean, comres.delta_mean[l].mean);
        scran_tests::compare_almost_equal(ref.delta_mean[l-1].median, comres.delta_mean[l].median);
        scran_tests::compare_almost_equal(ref.delta_mean[l-1].max, comres.delta_mean[l].max);
        scran_tests::compare_almost_equal(ref.delta_mean[l-1].min_rank, comres.delta_mean[l].min_rank);

        scran_tests::compare_almost_equal(ref.delta_detected[l-1].min, comres.delta_detected[l].min);
        scran_tests::compare_almost_equal(ref.delta_detected[l-1].mean, comres.delta_detected[l].mean);
        scran_tests::compare_almost_equal(ref.delta_detected[l-1].median, comres.delta_detected[l].median);
        scran_tests::compare_almost_equal(ref.delta_detected[l-1].max, comres.delta_detected[l].max);
        scran_tests::compare_almost_equal(ref.delta_detected[l-1].min_rank, comres.delta_detected[l].min_rank);
    }
}

TEST(ScoreMarkersSummaryScenarios, MinRank) {
    // Checking that the minimum rank is somewhat sensible,
    // and we didn't feed in the wrong values somewhere.
    int ngenes = 10;
    int nsamples = 8;

    std::vector<double> buffer(ngenes * nsamples, 0);
    for (int i = 0; i < ngenes; ++i) {
        buffer[i * nsamples + 1] = i; // second group gets increasing rank with later genes
        buffer[i * nsamples + 3] = i+1;

        buffer[i * nsamples + 4] = -i; // third group gets decreasing rank with later genes
        buffer[i * nsamples + 6] = -i + 1;
    }

    std::vector<int> grouping{0, 1, 0, 1, 2, 3, 2, 3};
    tatami::DenseRowMatrix<double, int> mat(ngenes, nsamples, std::move(buffer));

    scran_markers::ScoreMarkersSummaryOptions opt;
    auto res = scran_markers::score_markers_summary(mat, grouping.data(), opt);

    for (int i = 0; i < ngenes; ++i) {
        EXPECT_EQ(res.cohens_d[1].min_rank[i], ngenes - i); // second group
        EXPECT_EQ(res.cohens_d[2].min_rank[i], i + 1);  // third group
    }
}

TEST(ScoreMarkersSummaryScenarios, DisabledSummaries) {
    int ngenes = 10, nsamples = 40;
    tatami::DenseRowMatrix<double, int> mat(ngenes, nsamples, std::vector<double>(ngenes * nsamples));
    int ngroups = 4;
    std::vector<int> groupings = create_groupings(nsamples, ngroups);

    scran_markers::ScoreMarkersSummaryOptions opt;
    opt.compute_min = false;
    opt.compute_mean = false;
    opt.compute_median = false;
    opt.compute_max = false;
    opt.compute_min_rank = false;

    auto empty = scran_markers::score_markers_summary(mat, groupings.data(), opt);
    EXPECT_EQ(empty.mean.size(), ngroups);
    EXPECT_EQ(empty.detected.size(), ngroups);
    EXPECT_EQ(empty.cohens_d.size(), ngroups);
    EXPECT_EQ(empty.auc.size(), ngroups);
    EXPECT_EQ(empty.delta_mean.size(), ngroups);
    EXPECT_EQ(empty.delta_detected.size(), ngroups);
}

TEST(ScoreMarkersSummary, CapCacheSize) {
    EXPECT_EQ(scran_markers::internal::cap_cache_size(100, 1), 0);
    EXPECT_EQ(scran_markers::internal::cap_cache_size(100, 0), 0);
    EXPECT_EQ(scran_markers::internal::cap_cache_size(100, 10), 45);
    EXPECT_EQ(scran_markers::internal::cap_cache_size(100, 11), 55);
    EXPECT_EQ(scran_markers::internal::cap_cache_size(100, 21), 100);
    EXPECT_EQ(scran_markers::internal::cap_cache_size(100, 20), 100);
}

/*********************************************/

class ScoreMarkersSummaryOneAtATimeTest :
    public ScoreMarkersSummaryTestCore,
    public ::testing::TestWithParam<int>
{
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
                        scran_tests::SimulationParameters sparam;
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

TEST_P(ScoreMarkersSummaryOneAtATimeTest, Basic) {
    auto NC = dense_row->ncol();
    int ngroups = 3;
    std::vector<int> groupings = create_groupings(NC, ngroups);

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

    scran_markers::ScoreMarkersSummaryOptions opt;
    auto ref = scran_markers::score_markers_summary<double>(*mat, groupings.data(), opt);

    // Only the group mean.
    {
        scran_markers::ScoreMarkersSummaryOptions opt;
        opt.compute_group_detected = false;
        opt.compute_cohens_d = false;
        opt.compute_auc = false;
        opt.compute_delta_mean = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_summary<double>(*mat, groupings.data(), opt);
        compare_averages(ref.mean, alt.mean);
        EXPECT_TRUE(alt.detected.empty());
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only the group detected proportions.
    {
        scran_markers::ScoreMarkersSummaryOptions opt;
        opt.compute_group_mean = false;
        opt.compute_cohens_d = false;
        opt.compute_auc = false;
        opt.compute_delta_mean = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_summary<double>(*mat, groupings.data(), opt);
        compare_averages(ref.detected, alt.detected);
        EXPECT_TRUE(alt.mean.empty());
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only Cohen's d.
    {
        scran_markers::ScoreMarkersSummaryOptions opt;
        opt.compute_group_mean = false;
        opt.compute_group_detected = false;
        opt.compute_auc = false;
        opt.compute_delta_mean = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_summary<double>(*mat, groupings.data(), opt);
        compare_summaries_for_effect(ngroups, alt.cohens_d, ref.cohens_d);
        EXPECT_TRUE(alt.mean.empty());
        EXPECT_TRUE(alt.detected.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only AUC.
    {
        scran_markers::ScoreMarkersSummaryOptions opt;
        opt.compute_group_mean = false;
        opt.compute_group_detected = false;
        opt.compute_cohens_d = false;
        opt.compute_delta_mean = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_summary<double>(*mat, groupings.data(), opt);
        compare_summaries_for_effect(ngroups, alt.auc, ref.auc);
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only delta-mean.
    {
        scran_markers::ScoreMarkersSummaryOptions opt;
        opt.compute_group_mean = false;
        opt.compute_group_detected = false;
        opt.compute_cohens_d = false;
        opt.compute_auc = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_summary<double>(*mat, groupings.data(), opt);
        compare_summaries_for_effect(ngroups, alt.delta_mean, ref.delta_mean);
        EXPECT_TRUE(alt.mean.empty());
        EXPECT_TRUE(alt.detected.empty());
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only delta-mean.
    {
        scran_markers::ScoreMarkersSummaryOptions opt;
        opt.compute_group_mean = false;
        opt.compute_group_detected = false;
        opt.compute_cohens_d = false;
        opt.compute_auc = false;
        opt.compute_delta_mean = false;

        auto alt = scran_markers::score_markers_summary<double>(*mat, groupings.data(), opt);
        compare_summaries_for_effect(ngroups, alt.delta_detected, ref.delta_detected);
        EXPECT_TRUE(alt.mean.empty());
        EXPECT_TRUE(alt.detected.empty());
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
    }
}

INSTANTIATE_TEST_SUITE_P(
    ScoreMarkersSummary,
    ScoreMarkersSummaryOneAtATimeTest,
    ::testing::Values(0, 1, 2, 3)
);

