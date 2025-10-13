#include "scran_tests/scran_tests.hpp"
#include "tatami/tatami.hpp"
#include "topicks/topicks.hpp"

#include "scran_markers/score_markers_best.hpp"
#include "scran_markers/score_markers_pairwise.hpp"
#include "scran_markers/summarize_effects.hpp"

#include "utils.h"

class ScoreMarkersBestTestCore {
protected:
    template<class Left_, class Right_>
    static void compare_averages(const Left_& res, const Right_& other, int ngroups) {
        EXPECT_EQ(other.mean.size(), ngroups);
        EXPECT_EQ(other.detected.size(), ngroups);
        EXPECT_EQ(res.mean.size(), ngroups);
        EXPECT_EQ(res.detected.size(), ngroups);

        for (int l = 0; l < ngroups; ++l) {
            scran_tests::compare_almost_equal(res.mean[l], other.mean[l]);
            scran_tests::compare_almost_equal(res.detected[l], other.detected[l]);
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
                EXPECT_EQ(left[g1][g2], right[g1][g2]);
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

                topicks::PickTopGenesOptions<double> opt;
                opt.check_nan = true;
                opt.keep_ties = keep_ties;
                if (bound.has_value()) {
                    opt.bound = *bound;
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
                    } else {
                        return left.second > right.second;
                    }
                });

                EXPECT_EQ(expected, best[g1][g2]);
            }
        }
    }
};

/*********************************************/

class ScoreMarkersBestTest : public ScoreMarkersBestTestCore, public ::testing::TestWithParam<std::tuple<int, bool, int, bool, bool, bool, int> > {
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

    opt.compute_auc = do_auc; // false, if we want to check the running implementations.
    auto ref = scran_markers::score_markers_best<double>(*dense_row, groupings.data(), top, opt);

    if (nthreads == 1) {
        // Comparing the efficient ScoreMarkers implementation against the
        // PairwiseEffects + SummarizeEffects combination, which is less mind-bending
        // but requires holding a large 3D matrix in memory.
        scran_markers::ScoreMarkersPairwiseOptions popt;
        popt.compute_auc = do_auc;
        auto pairres = scran_markers::score_markers_pairwise(*dense_row, groupings.data(), popt);
        compare_averages(ref, pairres, ngroups);

        compare_best_to_pairwise(ref.cohens_d, pairres.cohens_d, ngenes, ngroups, top, opt.largest_cohens_d, keep_ties, opt.threshold_cohens_d);
        compare_best_to_pairwise(ref.delta_mean, pairres.delta_mean, ngenes, ngroups, top, opt.largest_delta_mean, keep_ties, opt.threshold_delta_mean);
        compare_best_to_pairwise(ref.delta_detected, pairres.delta_detected, ngenes, ngroups, top, opt.largest_delta_detected, keep_ties, opt.threshold_delta_detected);
        if (do_auc) {
            compare_best_to_pairwise(ref.auc, pairres.auc, ngenes, ngroups, top, opt.largest_auc, keep_ties, opt.threshold_auc);
        }

    } else {
        opt.num_threads = nthreads;
        auto dr = scran_markers::score_markers_best<double>(*dense_row, groupings.data(), top, opt);
        compare_averages(ref, dr, ngroups);
        compare_best(ref, dr);
    }

    // Comparing to all of the other matrix representations.
    {
        auto dc = scran_markers::score_markers_best<double>(*dense_column, groupings.data(), top, opt);
        compare_averages(ref, dc, ngroups);
        compare_best(ref, dc);

        auto sr = scran_markers::score_markers_best<double>(*sparse_row, groupings.data(), top, opt);
        compare_averages(ref, sr, ngroups);
        compare_best(ref, sr);

        auto sc = scran_markers::score_markers_best<double>(*sparse_column, groupings.data(), top, opt);
        compare_averages(ref, sc, ngroups);
        compare_best(ref, sc);
    }
}

INSTANTIATE_TEST_SUITE_P(
    ScoreMarkers,
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

///*********************************************/
//
//class ScoreMarkersBestBlockedTest : public ScoreMarkersBestTestCore, public ::testing::TestWithParam<std::tuple<int, bool, scran_blocks::WeightPolicy, int> > {
//protected:
//    inline static std::shared_ptr<tatami::Matrix<double, int> > dense_row, dense_column, sparse_row, sparse_column;
//
//    static void SetUpTestSuite() {
//        size_t nr = 248, nc = 287;
//        dense_row.reset(
//            new tatami::DenseRowMatrix<double, int>(
//                nr,
//                nc,
//                scran_tests::simulate_vector(
//                    nr * nc, 
//                    []{
//                        scran_tests::SimulationParameters sparam;
//                        sparam.density = 0.1;
//                        sparam.seed = 666;
//                        return sparam;
//                    }()
//                )
//            )
//        );
//
//        dense_column = tatami::convert_to_dense(dense_row.get(), false);
//        sparse_row = tatami::convert_to_compressed_sparse(dense_row.get(), true);
//        sparse_column = tatami::convert_to_compressed_sparse(dense_row.get(), false);
//    }
//};
//
//TEST_P(ScoreMarkersBestBlockedTest, AgainstPairwise) {
//    auto param = GetParam();
//    auto ngroups = std::get<0>(param);
//    bool do_auc = std::get<1>(param);
//    auto policy = std::get<2>(param);
//    auto nthreads = std::get<3>(param);
//
//    auto NC = dense_row->ncol();
//    std::vector<int> groupings = create_groupings(NC, ngroups);
//    std::vector<int> blocks = create_blocks(NC, 3);
//
//    scran_markers::ScoreMarkersBestOptions opt;
//    opt.compute_auc = do_auc; // false, if we want to check the running implementations.
//    opt.block_weight_policy = policy;
//    auto ref = scran_markers::score_markers_summary_blocked(*dense_row, groupings.data(), blocks.data(), opt);
//
//    if (nthreads == 1) {
//        scran_markers::ScoreMarkersPairwiseOptions popt;
//        popt.compute_auc = do_auc;
//        popt.block_weight_policy = policy;
//        auto pairres = scran_markers::score_markers_pairwise_blocked(*dense_row, groupings.data(), blocks.data(), popt);
//        compare_averages(ref, pairres, ngroups);
//
//        scran_markers::SummarizeEffectsOptions sopt;
//        size_t ngenes = dense_row->nrow();
//        auto cohen_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.cohens_d.data(), sopt);
//        compare_summaries_for_effect(ngroups, cohen_summ, ref.cohens_d);
//
//        auto dm_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.delta_mean.data(), sopt);
//        compare_summaries_for_effect(ngroups, dm_summ, ref.delta_mean);
//
//        auto dd_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.delta_detected.data(), sopt);
//        compare_summaries_for_effect(ngroups, dd_summ, ref.delta_detected);
//
//        if (do_auc) {
//            auto auc_summ = scran_markers::summarize_effects(ngenes, ngroups, pairres.auc.data(), sopt);
//            compare_summaries_for_effect(ngroups, auc_summ, ref.auc);
//        }
//
//    } else {
//        opt.num_threads = nthreads;
//        auto dr = scran_markers::score_markers_summary_blocked(*dense_row, groupings.data(), blocks.data(), opt);
//        compare_averages(ref, dr, ngroups);
//        compare_effects(ngroups, ref, dr, do_auc);
//    }
//
//    // Checking the other references.
//    {
//        auto dc = scran_markers::score_markers_summary_blocked(*dense_column, groupings.data(), blocks.data(), opt);
//        compare_averages(ref, dc, ngroups);
//        compare_effects(ngroups, ref, dc, do_auc);
//
//        auto sr = scran_markers::score_markers_summary_blocked(*sparse_row, groupings.data(), blocks.data(), opt);
//        compare_averages(ref, sr, ngroups);
//        compare_effects(ngroups, ref, sr, do_auc);
//
//        auto sc = scran_markers::score_markers_summary_blocked(*sparse_column, groupings.data(), blocks.data(), opt);
//        compare_averages(ref, sc, ngroups);
//        compare_effects(ngroups, ref, sc, do_auc);
//    }
//
//    // Comparing to what happens if we disable minrank, which allows for more efficient AUC calculations.
//    if (do_auc) {
//        scran_markers::ScoreMarkersBestOptions opt_no_mr = opt;
//        opt_no_mr.compute_min_rank = false;
//        wipe_min_rank(ref);
//
//        auto dr_nmr = scran_markers::score_markers_summary_blocked(*dense_row, groupings.data(), blocks.data(), opt_no_mr);
//        compare_effects(ngroups, ref, dr_nmr, true);
//
//        auto dc_nmr = scran_markers::score_markers_summary_blocked(*dense_column, groupings.data(), blocks.data(), opt_no_mr);
//        compare_effects(ngroups, ref, dc_nmr, true);
//
//        auto sr_nmr = scran_markers::score_markers_summary_blocked(*sparse_row, groupings.data(), blocks.data(), opt_no_mr);
//        compare_effects(ngroups, ref, sr_nmr, true);
//
//        auto sc_nmr = scran_markers::score_markers_summary_blocked(*sparse_column, groupings.data(), blocks.data(), opt_no_mr);
//        compare_effects(ngroups, ref, sc_nmr, true);
//    }
//}
//
//INSTANTIATE_TEST_SUITE_P(
//    ScoreMarkersBest,
//    ScoreMarkersBestBlockedTest,
//    ::testing::Combine(
//        ::testing::Values(2, 3, 4, 5), // number of clusters
//        ::testing::Values(false, true), // with or without the AUC?
//        ::testing::Values(scran_blocks::WeightPolicy::NONE, scran_blocks::WeightPolicy::EQUAL), // block weighting method.
//        ::testing::Values(1, 3) // number of threads
//    )
//);
