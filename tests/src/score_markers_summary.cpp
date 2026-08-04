#include "scran_tests/scran_tests.hpp"
#include "tatami/tatami.hpp"

#include "scran_markers/score_markers_summary.hpp"
#include "scran_markers/score_markers_pairwise.hpp"
#include "scran_markers/summarize_effects.hpp"

#include "utils.h"

static void compare_averages(const std::vector<std::vector<double> >& res, const std::vector<std::vector<double> >& other) {
    const int ngroups = res.size();
    ASSERT_EQ(ngroups, other.size());
    for (int l = 0; l < ngroups; ++l) {
        scran_tests::compare_almost_equal(res[l], other[l]);
    }
}

template<class Summaries_>
static void compare_summaries_for_effect(int ngroups, const std::vector<Summaries_>& res, const std::vector<Summaries_>& other) {
    ASSERT_EQ(res.size(), ngroups);
    ASSERT_EQ(other.size(), ngroups);

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
static void compare_limited_minrank(int limit, const std::vector<Summaries_>& ref, const std::vector<Summaries_>& sum) {
    const int ngroups = ref.size();
    for (int l = 0; l < ngroups; ++l) {
        const auto& refmr = ref[l].min_rank;
        const auto& summr = sum[l].min_rank;

        const int ngenes = refmr.size();
        for (int ge = 0; ge < ngenes; ++ge) {
            auto smr = summr[ge];
            auto rmr = refmr[ge];
            if (smr <= limit) {
                EXPECT_EQ(smr, rmr);
            } else {
                EXPECT_GT(rmr, limit);
                EXPECT_EQ(smr, ngenes);
            }
        }
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

        // Relaxing by a tolerance to account for numerical imprecision when averaging multiple identical values. 
        EXPECT_LE(curmin, curmean + 1e-8);
        EXPECT_LE(curmin, curmed);
        EXPECT_GE(curmax, curmean - 1e-8);
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

/*********************************************/

class ScoreMarkersSummaryTest : public ::testing::TestWithParam<std::tuple<int, bool, int> > {
protected:
    inline static std::shared_ptr<tatami::Matrix<double, int> > dense_row, dense_column, sparse_row, sparse_column;

    static void SetUpTestSuite() {
        size_t nr = 548, nc = 192;
        dense_row.reset(
            new tatami::DenseRowMatrix<double, int>(
                nr,
                nc,
                scran_tests::simulate_vector(
                    nr * nc, 
                    []{
                        scran_tests::SimulateVectorParameters sparam;
                        sparam.density = 0.15;
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
    const auto ngroups = std::get<0>(param);
    const auto do_auc = std::get<1>(param);
    const auto nthreads = std::get<2>(param);

    // Using contiguous groups for some variety in the test suite.
    std::vector<int> groupings = create_contiguous_factor(dense_row->ncol(), ngroups);

    scran_markers::ScoreMarkersSummaryOptions opt;
    opt.compute_auc = do_auc; // false, if we want to check the running implementations.
    auto ref = scran_markers::score_markers_summary(*dense_row, groupings.data(), ngroups, opt);

    if (nthreads == 1) {
        const auto ngenes = dense_row->nrow();

        // Running some further checks on the effects.
        for (int l = 0; l < ngroups; ++l) {
            check_effects(ngenes, l, ref.cohens_d);
            check_effects(ngenes, l, ref.delta_mean);

            check_effects(ngenes, l, ref.delta_detected);
            for (int g = 0; g < ngenes; ++g) {
                double curmin = ref.delta_detected[l].min[g];
                double curmax = ref.delta_detected[l].max[g];
                EXPECT_GE(curmin, -1);
                EXPECT_LE(curmax, 1);
            }

            if (do_auc) {
                check_effects(ngenes, l, ref.auc);
                for (int g = 0; g < ngenes; ++g) {
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
        auto pairres = scran_markers::score_markers_pairwise(*dense_row, groupings.data(), ngroups, popt);
        compare_averages(ref.mean, pairres.mean);
        compare_averages(ref.detected, pairres.detected);

        scran_markers::SummarizeEffectsOptions sopt;
        auto cohen_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.cohens_d.data(), sopt);
        compare_summaries_for_effect(ngroups, cohen_summ, ref.cohens_d);
        compare_limited_minrank(opt.min_rank_limit, cohen_summ, ref.cohens_d);

        auto dm_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.delta_mean.data(), sopt);
        compare_summaries_for_effect(ngroups, dm_summ, ref.delta_mean);
        compare_limited_minrank(opt.min_rank_limit, dm_summ, ref.delta_mean);

        auto dd_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.delta_detected.data(), sopt);
        compare_summaries_for_effect(ngroups, dd_summ, ref.delta_detected);
        compare_limited_minrank(opt.min_rank_limit, dd_summ, ref.delta_detected);

        if (do_auc) {
            auto auc_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.auc.data(), sopt);
            compare_summaries_for_effect(ngroups, auc_summ, ref.auc);
            compare_limited_minrank(opt.min_rank_limit, auc_summ, ref.auc);
        }

    } else {
        opt.num_threads = nthreads;
        auto dr = scran_markers::score_markers_summary(*dense_row, groupings.data(), ngroups, opt);
        compare_averages(ref.mean, dr.mean);
        compare_averages(ref.detected, dr.detected);
        compare_effects(ngroups, ref, dr, do_auc);
    }

    // Comparing to all of the other matrix representations.
    {
        auto dc = scran_markers::score_markers_summary(*dense_column, groupings.data(), ngroups, opt);
        compare_averages(ref.mean, dc.mean);
        compare_averages(ref.detected, dc.detected);
        compare_effects(ngroups, ref, dc, do_auc);

        auto sr = scran_markers::score_markers_summary(*sparse_row, groupings.data(), ngroups, opt);
        compare_averages(ref.mean, sr.mean);
        compare_averages(ref.detected, sr.detected);
        compare_effects(ngroups, ref, sr, do_auc);

        auto sc = scran_markers::score_markers_summary(*sparse_column, groupings.data(), ngroups, opt);
        compare_averages(ref.mean, sc.mean);
        compare_averages(ref.detected, sc.detected);
        compare_effects(ngroups, ref, sc, do_auc);
    }

    // Trying with a more limited min-rank count.
    {
        auto mr_limited_opt = opt;
        mr_limited_opt.min_rank_limit = 10;
        auto limited = scran_markers::score_markers_summary(*dense_row, groupings.data(), ngroups, mr_limited_opt);

        compare_limited_minrank(mr_limited_opt.min_rank_limit, ref.cohens_d, limited.cohens_d);
        compare_limited_minrank(mr_limited_opt.min_rank_limit, ref.delta_mean, limited.delta_mean);
        compare_limited_minrank(mr_limited_opt.min_rank_limit, ref.delta_detected, limited.delta_detected);
        if (do_auc) {
            compare_limited_minrank(mr_limited_opt.min_rank_limit, ref.auc, limited.auc);
        }
    }

    // Comparing to quantile summaries; should be the same for 1 block.
    {
        auto qopt = opt;
        qopt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;

        auto qdr = scran_markers::score_markers_summary(*dense_row, groupings.data(), ngroups, qopt);
        compare_averages(ref.mean, qdr.mean);
        compare_averages(ref.detected, qdr.detected);
        compare_effects(ngroups, ref, qdr, do_auc);

        auto qdc = scran_markers::score_markers_summary(*dense_column, groupings.data(), ngroups, qopt);
        compare_averages(ref.mean, qdc.mean);
        compare_averages(ref.detected, qdc.detected);
        compare_effects(ngroups, ref, qdc, do_auc);

        auto qsr = scran_markers::score_markers_summary(*sparse_row, groupings.data(), ngroups, qopt);
        compare_averages(ref.mean, qsr.mean);
        compare_averages(ref.detected, qsr.detected);
        compare_effects(ngroups, ref, qsr, do_auc);

        auto qsc = scran_markers::score_markers_summary(*sparse_column, groupings.data(), ngroups, qopt);
        compare_averages(ref.mean, qsc.mean);
        compare_averages(ref.detected, qsc.detected);
        compare_effects(ngroups, ref, qsc, do_auc);
    }
}

TEST_P(ScoreMarkersSummaryTest, SummaryQuantile) {
    auto param = GetParam();
    const auto ngroups = std::get<0>(param);
    const auto do_auc = std::get<1>(param);
    const auto nthreads = std::get<2>(param);

    std::vector<int> groupings = create_interleaved_factor(dense_row->ncol(), ngroups);

    scran_markers::ScoreMarkersSummaryOptions opt;
    opt.compute_auc = do_auc; // false, if we want to check the running implementations.
    opt.num_threads = nthreads;
    opt.compute_summary_quantiles = std::vector<double>{ 0.0, 0.5, 1.0 };
    auto out = scran_markers::score_markers_summary(*dense_row, groupings.data(), ngroups, opt);

    const auto ngenes = dense_row->nrow();
    for (int l = 0; l < ngroups; ++l) {
        for (int g = 0; g < ngenes; ++g) {
            const auto& cq = *(out.cohens_d[l].quantiles);
            EXPECT_EQ(out.cohens_d[l].min[g], cq[0][g]);
            scran_tests::compare_almost_equal(out.cohens_d[l].median[g], cq[1][g], scran_tests::CompareAlmostEqualParameters{});
            EXPECT_EQ(out.cohens_d[l].max[g], cq[2][g]);

            if (do_auc) {
                const auto& aq = *(out.auc[l].quantiles);
                EXPECT_EQ(out.auc[l].min[g], aq[0][g]);
                scran_tests::compare_almost_equal(out.auc[l].median[g], aq[1][g], scran_tests::CompareAlmostEqualParameters{});
                EXPECT_EQ(out.auc[l].max[g], aq[2][g]);
            }

            const auto& mq = *(out.delta_mean[l].quantiles);
            EXPECT_EQ(out.delta_mean[l].min[g], mq[0][g]);
            scran_tests::compare_almost_equal(out.delta_mean[l].median[g], mq[1][g], scran_tests::CompareAlmostEqualParameters{});
            EXPECT_EQ(out.delta_mean[l].max[g], mq[2][g]);

            const auto& dq = *(out.delta_detected[l].quantiles);
            EXPECT_EQ(out.delta_detected[l].min[g], dq[0][g]);
            scran_tests::compare_almost_equal(out.delta_detected[l].median[g], dq[1][g], scran_tests::CompareAlmostEqualParameters{});
            EXPECT_EQ(out.delta_detected[l].max[g], dq[2][g]);
        }
    }
}

INSTANTIATE_TEST_SUITE_P(
    ScoreMarkersSummary,
    ScoreMarkersSummaryTest,
    ::testing::Combine(
        ::testing::Values(2, 5), // number of clusters
        ::testing::Values(false, true), // with or without the AUC?
        ::testing::Values(1, 3) // number of threads
    )
);

/*********************************************/

class ScoreMarkersSummaryBlockedTest : public ::testing::TestWithParam<std::tuple<int, bool, scran_blocks::WeightPolicy, int> > {
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
                        scran_tests::SimulateVectorParameters sparam;
                        sparam.density = 0.18;
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

TEST_P(ScoreMarkersSummaryBlockedTest, ReferenceMean) {
    auto param = GetParam();
    auto ngroups = std::get<0>(param);
    bool do_auc = std::get<1>(param);
    auto policy = std::get<2>(param);
    auto nthreads = std::get<3>(param);

    const auto ncols = dense_row->ncol();
    const int nblocks = 3;
    auto groupings = create_interleaved_factor(ncols, ngroups);
    auto blocks = create_contiguous_factor(ncols, nblocks);

    scran_markers::ScoreMarkersSummaryOptions opt;
    opt.compute_auc = do_auc; // false, if we want to check the running implementations.
    opt.block_weight_policy = policy;
    auto ref = scran_markers::score_markers_summary_blocked(*dense_row, groupings.data(), ngroups, blocks.data(), nblocks, opt);

    if (nthreads == 1) {
        scran_markers::ScoreMarkersPairwiseOptions popt;
        popt.compute_auc = do_auc;
        popt.block_weight_policy = policy;
        auto pairres = scran_markers::score_markers_pairwise_blocked(*dense_row, groupings.data(), ngroups, blocks.data(), nblocks, popt);
        compare_averages(ref.mean, pairres.mean);
        compare_averages(ref.detected, pairres.detected);

        scran_markers::SummarizeEffectsOptions sopt;
        const auto ngenes = dense_row->nrow();
        auto cohen_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.cohens_d.data(), sopt);
        compare_summaries_for_effect(ngroups, cohen_summ, ref.cohens_d);
        compare_limited_minrank(opt.min_rank_limit, cohen_summ, ref.cohens_d);

        auto dm_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.delta_mean.data(), sopt);
        compare_summaries_for_effect(ngroups, dm_summ, ref.delta_mean);
        compare_limited_minrank(opt.min_rank_limit, dm_summ, ref.delta_mean);

        auto dd_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.delta_detected.data(), sopt);
        compare_summaries_for_effect(ngroups, dd_summ, ref.delta_detected);
        compare_limited_minrank(opt.min_rank_limit, dd_summ, ref.delta_detected);

        if (do_auc) {
            auto auc_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.auc.data(), sopt);
            compare_summaries_for_effect(ngroups, auc_summ, ref.auc);
            compare_limited_minrank(opt.min_rank_limit, auc_summ, ref.auc);
        }

    } else {
        opt.num_threads = nthreads;
        auto dr = scran_markers::score_markers_summary_blocked(*dense_row, groupings.data(), ngroups, blocks.data(), nblocks, opt);
        compare_averages(ref.mean, dr.mean);
        compare_averages(ref.detected, dr.detected);
        compare_effects(ngroups, ref, dr, do_auc);
    }

    // Checking the other references.
    {
        auto dc = scran_markers::score_markers_summary_blocked(*dense_column, groupings.data(), ngroups, blocks.data(), nblocks, opt);
        compare_averages(ref.mean, dc.mean);
        compare_averages(ref.detected, dc.detected);
        compare_effects(ngroups, ref, dc, do_auc);

        auto sr = scran_markers::score_markers_summary_blocked(*sparse_row, groupings.data(), ngroups, blocks.data(), nblocks, opt);
        compare_averages(ref.mean, sr.mean);
        compare_averages(ref.detected, sr.detected);
        compare_effects(ngroups, ref, sr, do_auc);

        auto sc = scran_markers::score_markers_summary_blocked(*sparse_column, groupings.data(), ngroups, blocks.data(), nblocks, opt);
        compare_averages(ref.mean, sc.mean);
        compare_averages(ref.detected, sc.detected);
        compare_effects(ngroups, ref, sc, do_auc);
    }

    {
        // Trying with a more limited min-rank count.
        auto mr_limited_opt = opt;
        mr_limited_opt.min_rank_limit = 10;
        auto limited = scran_markers::score_markers_summary_blocked(*dense_row, groupings.data(), ngroups, blocks.data(), nblocks, mr_limited_opt);

        compare_limited_minrank(mr_limited_opt.min_rank_limit, ref.cohens_d, limited.cohens_d);
        compare_limited_minrank(mr_limited_opt.min_rank_limit, ref.delta_mean, limited.delta_mean);
        compare_limited_minrank(mr_limited_opt.min_rank_limit, ref.delta_detected, limited.delta_detected);
        if (do_auc) {
            compare_limited_minrank(mr_limited_opt.min_rank_limit, ref.auc, limited.auc);
        }
    }
}

TEST_P(ScoreMarkersSummaryBlockedTest, ReferenceQuantile) {
    auto param = GetParam();
    auto ngroups = std::get<0>(param);
    bool do_auc = std::get<1>(param);
    auto policy = std::get<2>(param);
    auto nthreads = std::get<3>(param);

    // Checking contiguous groupings with interleaved blocks, just for some variety.
    const auto ncols = dense_row->ncol();
    auto groupings = create_contiguous_factor(ncols, ngroups);
    const int nblocks = 3;
    auto blocks = create_interleaved_factor(ncols, nblocks);

    scran_markers::ScoreMarkersSummaryOptions opt;
    opt.compute_auc = do_auc; // false, if we want to check the running implementations.
    opt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
    opt.block_weight_policy = policy; // shouldn't really matter as quantile isn't weighted, but we might as well test it.
    auto ref = scran_markers::score_markers_summary_blocked(*dense_row, groupings.data(), ngroups, blocks.data(), nblocks, opt);

    if (nthreads == 1) {
        scran_markers::ScoreMarkersPairwiseOptions popt;
        popt.compute_auc = do_auc;
        popt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
        auto pairres = scran_markers::score_markers_pairwise_blocked(*dense_row, groupings.data(), ngroups, blocks.data(), nblocks, popt);
        compare_averages(ref.mean, pairres.mean);
        compare_averages(ref.detected, pairres.detected);

        scran_markers::SummarizeEffectsOptions sopt;
        const auto ngenes = dense_row->nrow();
        auto cohen_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.cohens_d.data(), sopt);
        compare_summaries_for_effect(ngroups, cohen_summ, ref.cohens_d);
        compare_limited_minrank(opt.min_rank_limit, cohen_summ, ref.cohens_d);

        auto dm_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.delta_mean.data(), sopt);
        compare_summaries_for_effect(ngroups, dm_summ, ref.delta_mean);
        compare_limited_minrank(opt.min_rank_limit, dm_summ, ref.delta_mean);

        auto dd_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.delta_detected.data(), sopt);
        compare_summaries_for_effect(ngroups, dd_summ, ref.delta_detected);
        compare_limited_minrank(opt.min_rank_limit, dd_summ, ref.delta_detected);

        if (do_auc) {
            auto auc_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.auc.data(), sopt);
            compare_summaries_for_effect(ngroups, auc_summ, ref.auc);
            compare_limited_minrank(opt.min_rank_limit, auc_summ, ref.auc);
        }

    } else {
        opt.num_threads = nthreads;
        auto qdr = scran_markers::score_markers_summary_blocked(*dense_row, groupings.data(), ngroups, blocks.data(), nblocks, opt);
        compare_averages(ref.mean, qdr.mean);
        compare_averages(ref.detected, qdr.detected);
        compare_effects(ngroups, ref, qdr, do_auc);
    }

    auto qdc = scran_markers::score_markers_summary_blocked(*dense_column, groupings.data(), ngroups, blocks.data(), nblocks, opt);
    compare_averages(ref.mean, qdc.mean);
    compare_averages(ref.detected, qdc.detected);
    compare_effects(ngroups, ref, qdc, do_auc);

    auto qsr = scran_markers::score_markers_summary_blocked(*sparse_row, groupings.data(), ngroups, blocks.data(), nblocks, opt);
    compare_averages(ref.mean, qsr.mean);
    compare_averages(ref.detected, qsr.detected);
    compare_effects(ngroups, ref, qsr, do_auc);

    auto qsc = scran_markers::score_markers_summary_blocked(*sparse_column, groupings.data(), ngroups, blocks.data(), nblocks, opt);
    compare_averages(ref.mean, qsc.mean);
    compare_averages(ref.detected, qsc.detected);
    compare_effects(ngroups, ref, qsc, do_auc);
}

INSTANTIATE_TEST_SUITE_P(
    ScoreMarkersSummary,
    ScoreMarkersSummaryBlockedTest,
    ::testing::Combine(
        ::testing::Values(2, 5), // number of clusters
        ::testing::Values(false, true), // with or without the AUC?
        ::testing::Values(scran_blocks::WeightPolicy::NONE, scran_blocks::WeightPolicy::EQUAL), // block weighting method.
        ::testing::Values(1, 3) // number of threads
    )
);

/*********************************************/

TEST(ScoreMarkersSummary, Thresholds) {
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

    const double threshold = 0.45;
    scran_markers::ScoreMarkersSummaryOptions opt;
    opt.threshold = threshold;
    auto out = scran_markers::score_markers_summary(mat, groupings.data(), ngroups, opt);

    // Comparing to pairwise + summarize_effects.
    {
        scran_markers::ScoreMarkersPairwiseOptions popt;
        popt.threshold = threshold;
        auto pairres = scran_markers::score_markers_pairwise(mat, groupings.data(), ngroups, popt);
        compare_averages(out.mean, pairres.mean);
        compare_averages(out.detected, pairres.detected);

        scran_markers::SummarizeEffectsOptions sopt;
        auto cohen_summ = scran_markers::summarize_effects(nrows, ngroups, pairres.cohens_d.data(), sopt);
        compare_summaries_for_effect(ngroups, cohen_summ, out.cohens_d);
        compare_limited_minrank(opt.min_rank_limit, cohen_summ, out.cohens_d);

        auto dm_summ = scran_markers::summarize_effects(nrows, ngroups, pairres.delta_mean.data(), sopt);
        compare_summaries_for_effect(ngroups, dm_summ, out.delta_mean);
        compare_limited_minrank(opt.min_rank_limit, dm_summ, out.delta_mean);

        auto dd_summ = scran_markers::summarize_effects(nrows, ngroups, pairres.delta_detected.data(), sopt);
        compare_summaries_for_effect(ngroups, dd_summ, out.delta_detected);
        compare_limited_minrank(opt.min_rank_limit, dd_summ, out.delta_detected);
    }

    // Check that the results are different from the no-threshold case.
    auto nothresh = scran_markers::score_markers_summary(mat, groupings.data(), ngroups, {});

    int some_diff = 0;
    for (int l = 0; l < ngroups; ++l) {
        for (int g = 0; g < nrows; ++g) {
            EXPECT_GT(nothresh.cohens_d[l].min[g], out.cohens_d[l].min[g]);
            EXPECT_GT(nothresh.cohens_d[l].mean[g], out.cohens_d[l].mean[g]);
            EXPECT_GT(nothresh.cohens_d[l].median[g], out.cohens_d[l].median[g]);
            EXPECT_GT(nothresh.cohens_d[l].max[g], out.cohens_d[l].max[g]);

            // '>' is not guaranteed due to imprecision with ranks... but (see below).
            some_diff += (nothresh.auc[l].mean[g] > out.auc[l].mean[g]);
            EXPECT_GE(nothresh.auc[l].min[g], out.auc[l].min[g]);
            EXPECT_GE(nothresh.auc[l].mean[g], out.auc[l].mean[g]);
            EXPECT_GE(nothresh.auc[l].median[g], out.auc[l].median[g]);
            EXPECT_GE(nothresh.auc[l].max[g], out.auc[l].max[g]);

            // These are not affected.
            EXPECT_EQ(nothresh.delta_mean[l].mean[g], out.delta_mean[l].mean[g]);
            EXPECT_EQ(nothresh.delta_detected[l].mean[g], out.delta_detected[l].mean[g]);
        }

        // Also not affected.
        EXPECT_EQ(nothresh.mean[l], out.mean[l]);
        EXPECT_EQ(nothresh.detected[l], out.detected[l]);
    }

    EXPECT_GT(some_diff, 0); // (from above)... at least one is '>', hopefully.

    // Quantile should give the same results for a single block.
    auto qopt = opt;
    qopt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
    auto qout = scran_markers::score_markers_summary(mat, groupings.data(), ngroups, qopt);
    compare_averages(out.mean, qout.mean);
    compare_averages(out.detected, qout.detected);
    compare_effects(ngroups, out, qout, true);
}

TEST(ScoreMarkersSummary, EmptyGroups) {
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

    scran_markers::ScoreMarkersSummaryOptions opt;
    auto ref = scran_markers::score_markers_summary(mat, groupings.data(), ngroups, opt);

    // Zero is effectively the missing group here.
    for (auto& g : groupings) {
        ++g;
    }
    const int ngroups_p1 = ngroups + 1;
    auto lost = scran_markers::score_markers_summary(mat, groupings.data(), ngroups_p1, opt);

    // First group is effectively all-NA.
    for (int g = 0; g < nrows; ++g) {
        EXPECT_TRUE(std::isnan(lost.cohens_d[0].min[g]));
        EXPECT_TRUE(std::isnan(lost.cohens_d[0].max[g]));
        EXPECT_TRUE(std::isnan(lost.cohens_d[0].mean[g]));
        EXPECT_TRUE(std::isnan(lost.cohens_d[0].median[g]));
        EXPECT_EQ(lost.cohens_d[0].min_rank[g], nrows);

        EXPECT_TRUE(std::isnan(lost.auc[0].min[g]));
        EXPECT_TRUE(std::isnan(lost.auc[0].max[g]));
        EXPECT_TRUE(std::isnan(lost.auc[0].mean[g]));
        EXPECT_TRUE(std::isnan(lost.auc[0].median[g]));
        EXPECT_EQ(lost.auc[0].min_rank[g], nrows);

        EXPECT_TRUE(std::isnan(lost.delta_mean[0].min[g]));
        EXPECT_TRUE(std::isnan(lost.delta_mean[0].max[g]));
        EXPECT_TRUE(std::isnan(lost.delta_mean[0].mean[g]));
        EXPECT_TRUE(std::isnan(lost.delta_mean[0].median[g]));
        EXPECT_EQ(lost.delta_mean[0].min_rank[g], nrows);

        EXPECT_TRUE(std::isnan(lost.delta_detected[0].min[g]));
        EXPECT_TRUE(std::isnan(lost.delta_detected[0].max[g]));
        EXPECT_TRUE(std::isnan(lost.delta_detected[0].mean[g]));
        EXPECT_TRUE(std::isnan(lost.delta_detected[0].median[g]));
        EXPECT_EQ(lost.delta_detected[0].min_rank[g], nrows);
    }

    // Expect all but the first group to give the same results.
    {
        auto copy = lost;
        copy.mean.erase(copy.mean.begin());
        compare_averages(ref.mean, copy.mean);

        copy.detected.erase(copy.detected.begin());
        compare_averages(ref.detected, copy.detected);

        copy.cohens_d.erase(copy.cohens_d.begin());
        copy.auc.erase(copy.auc.begin());
        copy.delta_mean.erase(copy.delta_mean.begin());
        copy.delta_detected.erase(copy.delta_detected.begin());
        compare_effects(ngroups, ref, copy, true);
    }

    // Quantile should give the same results for a single block.
    auto qopt = opt;
    qopt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
    auto qlost = scran_markers::score_markers_summary(mat, groupings.data(), ngroups_p1, qopt);
    compare_averages(lost.mean, qlost.mean);
    compare_averages(lost.detected, qlost.detected);
    compare_effects(ngroups + 1, lost, qlost, true);
}

TEST(ScoreMarkersSummary, BlockConfounded) {
    int nrows = 103, ncols = 290;
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

    // Block is fully confounded with the first group.
    std::vector<int> blocks(ncols);
    for (int c = 0; c < ncols; ++c) {
        blocks[c] = (groupings[c] == 0);
    }

    scran_markers::ScoreMarkersSummaryOptions sopt;
    auto comres = scran_markers::score_markers_summary_blocked(mat, groupings.data(), ngroups, blocks.data(), 2, sopt);

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
    auto ref = scran_markers::score_markers_summary(*sub, subgroups.data(), ngroups - 1, sopt); 

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
    {
        auto copy = comres;
        copy.mean.erase(copy.mean.begin());
        compare_averages(ref.mean, copy.mean);

        copy.detected.erase(copy.detected.begin());
        compare_averages(ref.detected, copy.detected);

        copy.cohens_d.erase(copy.cohens_d.begin());
        copy.auc.erase(copy.auc.begin());
        copy.delta_mean.erase(copy.delta_mean.begin());
        copy.delta_detected.erase(copy.delta_detected.begin());
        compare_effects(ngroups - 1, ref, copy, true);
    }

    // Quantile should give the same results, as there's basically only one block;
    // the second block is fully confounded.
    {
        auto qopt = sopt;
        qopt.block_average_policy = scran_markers::BlockAveragePolicy::QUANTILE;
        auto qcomres = scran_markers::score_markers_summary_blocked(mat, groupings.data(), ngroups, blocks.data(), 2, qopt);
        compare_averages(comres.mean, qcomres.mean);
        compare_averages(comres.detected, qcomres.detected);
        compare_effects(ngroups, comres, qcomres, true);
    }
}

/*********************************************/

static void check_one_at_a_time(
    const tatami::Matrix<double, int>& mat,
    const std::vector<int>& groupings,
    const std::size_t ngroups,
    const scran_markers::ScoreMarkersSummaryResults<double, int>& ref
) {
    // Only the group mean.
    {
        scran_markers::ScoreMarkersSummaryOptions opt;
        opt.compute_group_detected = false;
        opt.compute_cohens_d = false;
        opt.compute_auc = false;
        opt.compute_delta_mean = false;
        opt.compute_delta_detected = false;

        auto alt = scran_markers::score_markers_summary<double>(mat, groupings.data(), ngroups, opt);
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

        auto alt = scran_markers::score_markers_summary<double>(mat, groupings.data(), ngroups, opt);
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

        auto alt = scran_markers::score_markers_summary<double>(mat, groupings.data(), ngroups, opt);
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

        auto alt = scran_markers::score_markers_summary<double>(mat, groupings.data(), ngroups, opt);
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

        auto alt = scran_markers::score_markers_summary<double>(mat, groupings.data(), ngroups, opt);
        compare_summaries_for_effect(ngroups, alt.delta_mean, ref.delta_mean);
        EXPECT_TRUE(alt.mean.empty());
        EXPECT_TRUE(alt.detected.empty());
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }

    // Only delta-detected.
    {
        scran_markers::ScoreMarkersSummaryOptions opt;
        opt.compute_group_mean = false;
        opt.compute_group_detected = false;
        opt.compute_cohens_d = false;
        opt.compute_auc = false;
        opt.compute_delta_mean = false;

        auto alt = scran_markers::score_markers_summary<double>(mat, groupings.data(), ngroups, opt);
        compare_summaries_for_effect(ngroups, alt.delta_detected, ref.delta_detected);
        EXPECT_TRUE(alt.mean.empty());
        EXPECT_TRUE(alt.detected.empty());
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
    }

    // Nothing at all.
    {
        scran_markers::ScoreMarkersSummaryOptions opt;
        opt.compute_group_mean = false;
        opt.compute_group_detected = false;
        opt.compute_cohens_d = false;
        opt.compute_auc = false;
        opt.compute_delta_mean = false;
        opt.compute_delta_detected= false;

        auto alt = scran_markers::score_markers_summary<double>(mat, groupings.data(), ngroups, opt);
        EXPECT_TRUE(alt.mean.empty());
        EXPECT_TRUE(alt.detected.empty());
        EXPECT_TRUE(alt.cohens_d.empty());
        EXPECT_TRUE(alt.auc.empty());
        EXPECT_TRUE(alt.delta_mean.empty());
        EXPECT_TRUE(alt.delta_detected.empty());
    }
}

TEST_F(ScoreMarkersSummaryTest, OneAtATime) {
    const int NR = 228, NC = 318;
    auto dense_row = std::make_shared<tatami::DenseRowMatrix<double, int> >(
        NR,
        NC,
        scran_tests::simulate_vector(
            NR * NC, 
            []{
                scran_tests::SimulateVectorParameters sparam;
                sparam.density = 0.2;
                sparam.seed = 96;
                return sparam;
            }()
        )
    );

    const auto dense_column = tatami::convert_to_dense(dense_row.get(), false);
    const auto sparse_row = tatami::convert_to_compressed_sparse(dense_row.get(), true);
    const auto sparse_column = tatami::convert_to_compressed_sparse(dense_row.get(), false);

    int ngroups = 3;
    std::vector<int> groupings = create_groupings(NC, ngroups);
    scran_markers::ScoreMarkersSummaryOptions opt;
    auto ref = scran_markers::score_markers_summary<double>(*dense_row, groupings.data(), ngroups, opt);

    check_one_at_a_time(*dense_row, groupings, ngroups, ref);
    check_one_at_a_time(*sparse_row, groupings, ngroups, ref);
    check_one_at_a_time(*dense_column, groupings, ngroups, ref);
    check_one_at_a_time(*sparse_column, groupings, ngroups, ref);
}

/*********************************************/

TEST(ScoreMarkersSummary, MinRank) {
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
    auto res = scran_markers::score_markers_summary(mat, grouping.data(), 3, opt);

    for (int i = 0; i < ngenes; ++i) {
        EXPECT_EQ(res.cohens_d[1].min_rank[i], ngenes - i); // second group
        EXPECT_EQ(res.cohens_d[2].min_rank[i], i + 1);  // third group
    }
}

TEST(ScoreMarkersSummary, DisabledSummaries) {
    int ngenes = 10, nsamples = 40;
    tatami::DenseRowMatrix<double, int> mat(
        ngenes,
        nsamples, 
        scran_tests::simulate_vector(
            ngenes * nsamples,
            []{
                scran_tests::SimulateVectorParameters sparam;
                sparam.seed = 64;
                return sparam;
            }()
        )
    );

    int ngroups = 4;
    std::vector<int> groupings = create_groupings(nsamples, ngroups);

    // We won't test these one-by-one as this is already done in summarize_effects(). 
    scran_markers::ScoreMarkersSummaryOptions opt;
    auto ref = scran_markers::score_markers_summary(mat, groupings.data(), ngroups, opt);

    opt.compute_summary_min = false;
    opt.compute_summary_mean = false;
    opt.compute_summary_median = false;
    opt.compute_summary_max = false;
    opt.compute_summary_min_rank = false;
    auto empty = scran_markers::score_markers_summary(mat, groupings.data(), ngroups, opt);

    EXPECT_EQ(empty.mean.size(), ngroups);
    EXPECT_EQ(empty.detected.size(), ngroups);
    EXPECT_EQ(empty.cohens_d.size(), ngroups);
    EXPECT_EQ(empty.auc.size(), ngroups);
    EXPECT_EQ(empty.delta_mean.size(), ngroups);
    EXPECT_EQ(empty.delta_detected.size(), ngroups);

    for (int g = 0; g < ngroups; ++g) {
        EXPECT_EQ(empty.mean[g], ref.mean[g]);
        EXPECT_EQ(empty.detected[g], ref.detected[g]);

        EXPECT_TRUE(empty.cohens_d[g].min.empty());
        EXPECT_TRUE(empty.cohens_d[g].mean.empty());
        EXPECT_TRUE(empty.cohens_d[g].median.empty());
        EXPECT_TRUE(empty.cohens_d[g].max.empty());
        EXPECT_TRUE(empty.cohens_d[g].min_rank.empty());

        EXPECT_TRUE(empty.auc[g].min.empty());
        EXPECT_TRUE(empty.auc[g].mean.empty());
        EXPECT_TRUE(empty.auc[g].median.empty());
        EXPECT_TRUE(empty.auc[g].max.empty());
        EXPECT_TRUE(empty.auc[g].min_rank.empty());

        EXPECT_TRUE(empty.delta_mean[g].min.empty());
        EXPECT_TRUE(empty.delta_mean[g].mean.empty());
        EXPECT_TRUE(empty.delta_mean[g].median.empty());
        EXPECT_TRUE(empty.delta_mean[g].max.empty());
        EXPECT_TRUE(empty.delta_mean[g].min_rank.empty());

        EXPECT_TRUE(empty.delta_detected[g].min.empty());
        EXPECT_TRUE(empty.delta_detected[g].mean.empty());
        EXPECT_TRUE(empty.delta_detected[g].median.empty());
        EXPECT_TRUE(empty.delta_detected[g].max.empty());
        EXPECT_TRUE(empty.delta_detected[g].min_rank.empty());
    }
}

TEST(ScoreMarkersSummary, TiedMinRank) {
    // Simple set-up with two samples and unique delta-means for all genes.
    int nrows = 203, ncols = 2;
    std::vector<double> original(nrows * ncols);
    for (int r = 0; r < nrows; ++r) {
        original[r * 2] = r * 1.5;
    }
    tatami::DenseRowMatrix<double, int> omat(nrows, ncols, original);

    // The idea is to duplicate the matrix and check that the duplicated genes show up in the min_rank vector. 
    auto duplicated = original;
    duplicated.insert(duplicated.end(), original.begin(), original.end());
    const int nrows2 = 2 * nrows;
    tatami::DenseRowMatrix<double, int> dmat(nrows2, ncols, duplicated);

    scran_markers::ScoreMarkersSummaryOptions sopt;
    sopt.compute_cohens_d = false;
    sopt.compute_auc = false;
    sopt.compute_delta_detected = false;

    const int min_rank_limit = 10;
    sopt.min_rank_limit = min_rank_limit;
    sopt.min_rank_preserve_ties = true;

    std::vector<int> groupings{ 0, 1 };
    auto oout = scran_markers::score_markers_summary(omat, groupings.data(), 2, sopt);
    auto dout = scran_markers::score_markers_summary(dmat, groupings.data(), 2, sopt);

    for (int g = 0; g < 2; ++g) {
        std::vector<std::pair<int, int> > expected;
        for (int r = 0; r < nrows; ++r) {
            auto orank = oout.delta_mean[g].min_rank[r];
            if (orank <= min_rank_limit) {
                auto drank = (orank - 1) * 2 + 1;
                if (drank <= min_rank_limit) {
                    expected.emplace_back(drank, r);
                    expected.emplace_back(drank, r + nrows);
                }
            } else {
                EXPECT_EQ(orank, nrows);
            }
        }
        std::sort(expected.begin(), expected.end());

        std::vector<std::pair<int, int> > found_ties;
        for (int r = 0; r < nrows2; ++r) {
            auto drank = dout.delta_mean[g].min_rank[r];
            if (drank <= min_rank_limit) {
                found_ties.emplace_back(drank, r);
            } else {
                EXPECT_EQ(drank, nrows2);
            }
        }
        std::sort(found_ties.begin(), found_ties.end());

        EXPECT_EQ(expected, found_ties);
    }
}

TEST(ScoreMarkersSummary, UnsortedQuantiles) {
    tatami::DenseRowMatrix<double, int, std::vector<double> > empty(0, 0, std::vector<double>());
    scran_markers::ScoreMarkersSummaryOptions opts;
    opts.compute_summary_quantiles = std::vector<double>{ 0.5, 0.2 };
    scran_tests::expect_error([&]() -> void {
        scran_markers::score_markers_summary(empty, static_cast<int*>(NULL), 0, opts);
    }, "should be sorted");
}

TEST(ScoreMarkersSummary, NoGenes) {
    int nrows = 0, ncols = 91;
    tatami::DenseMatrix<double, int, std::vector<double> > mat(nrows, ncols, std::vector<double>(), true);

    int ngroups = 4;
    std::vector<int> groupings = create_groupings(ncols, ngroups);

    scran_markers::ScoreMarkersSummaryOptions opts;
    auto out = scran_markers::score_markers_summary(mat, groupings.data(), ngroups, opts);

    for (int g = 0; g < ngroups; ++g) {
        EXPECT_TRUE(out.mean[g].empty());
        EXPECT_TRUE(out.detected[g].empty());

        EXPECT_TRUE(out.cohens_d[g].mean.empty());
        EXPECT_TRUE(out.cohens_d[g].median.empty());
        EXPECT_TRUE(out.cohens_d[g].min.empty());
        EXPECT_TRUE(out.cohens_d[g].max.empty());
        EXPECT_TRUE(out.cohens_d[g].min_rank.empty());

        EXPECT_TRUE(out.auc[g].mean.empty());
        EXPECT_TRUE(out.auc[g].median.empty());
        EXPECT_TRUE(out.auc[g].min.empty());
        EXPECT_TRUE(out.auc[g].max.empty());
        EXPECT_TRUE(out.auc[g].min_rank.empty());

        EXPECT_TRUE(out.delta_mean[g].mean.empty());
        EXPECT_TRUE(out.delta_mean[g].median.empty());
        EXPECT_TRUE(out.delta_mean[g].min.empty());
        EXPECT_TRUE(out.delta_mean[g].max.empty());
        EXPECT_TRUE(out.delta_mean[g].min_rank.empty());

        EXPECT_TRUE(out.delta_detected[g].mean.empty());
        EXPECT_TRUE(out.delta_detected[g].median.empty());
        EXPECT_TRUE(out.delta_detected[g].min.empty());
        EXPECT_TRUE(out.delta_detected[g].max.empty());
        EXPECT_TRUE(out.delta_detected[g].min_rank.empty());
    }
}
