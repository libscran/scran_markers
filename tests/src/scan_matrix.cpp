#include "scran_tests/scran_tests.hpp"
#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"

#include "scran_markers/scan_matrix.hpp"
#include "scran_markers/create_combinations.hpp"

#include "utils.h"

class ScanMatrixSimpleTest : public ::testing::TestWithParam<std::tuple<int, int> > {
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
                        sparam.density = 0.2;
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

template<typename Group_>
static void simple_reference(
    const tatami::Matrix<double, int>& mat,
    const Group_* group,
    std::size_t num_groups,
    std::vector<double>& output_means,
    std::vector<double>& output_vars,
    std::vector<double>& output_detected
) {
    auto var_out = tatami_stats::group_variance(true, mat, group, num_groups, {});
    auto nonzero = tatami::DelayedUnaryIsometricOperation<double, double, int>(
        tatami::wrap_shared_ptr(&mat), 
        std::make_shared<tatami::DelayedUnaryIsometricCompareScalarHelper<tatami::CompareOperation::NOT_EQUAL, double, double, int, int> >(0)
    );
    auto det_out = tatami_stats::group_sum(true, nonzero, group, num_groups, {});

    const int num_genes = mat.nrow();
    auto group_sizes = scran_markers::tabulate_groups(mat.ncol(), group, num_groups);
    for (int r = 0; r < num_genes; ++r) {
        for (std::size_t g = 0; g < num_groups; ++g) {
            const auto pos = r * num_groups + g;
            output_means[pos] = var_out.mean[g][r];        
            output_vars[pos] = var_out.variance[g][r];        
            output_detected[pos] = det_out[g][r] / group_sizes[g];        
        }
    }
}

static void quick_compare_to_all_implementations(
    const std::shared_ptr<tatami::Matrix<double, int> >& dense_row,
    const std::shared_ptr<tatami::Matrix<double, int> >& dense_column, 
    const std::shared_ptr<tatami::Matrix<double, int> >& sparse_row,
    const std::shared_ptr<tatami::Matrix<double, int> >& sparse_column,
    const std::vector<int>& groupings,
    const int ngroups,
    const int nthreads,
    const std::vector<double>& ref_means,
    const std::vector<double>& ref_vars,
    const std::vector<double>& ref_detected
) {
    const auto full_size = ngroups * dense_row->nrow();
    auto gsizes = scran_markers::tabulate_groups(dense_row->ncol(), groupings.data(), ngroups);

    {
        std::vector<double> dr_means(full_size), dr_vars(full_size), dr_detected(full_size);
        scran_markers::internal::scan_matrix_by_row_full_auc<true>(
            *dense_row,
            groupings.data(),
            ngroups,
            static_cast<int*>(NULL),
            1,
            static_cast<int*>(NULL),
            ngroups,
            gsizes,
            scran_markers::internal::BlockAverageInfo<double>(),
            dr_means,
            dr_vars,
            dr_detected,
            static_cast<double*>(NULL),
            0.0,
            nthreads
        );

        scran_tests::compare_almost_equal_containers(ref_means, dr_means, {});
        scran_tests::compare_almost_equal_containers(ref_vars, dr_vars, {});
        scran_tests::compare_almost_equal_containers(ref_detected, dr_detected, {});
    }

    {
        std::vector<double> sr_means(full_size), sr_vars(full_size), sr_detected(full_size);
        scran_markers::internal::scan_matrix_by_row_full_auc<true>(
            *sparse_row,
            groupings.data(),
            ngroups,
            static_cast<int*>(NULL),
            1,
            static_cast<int*>(NULL),
            ngroups,
            gsizes,
            scran_markers::internal::BlockAverageInfo<double>(),
            sr_means,
            sr_vars,
            sr_detected,
            static_cast<double*>(NULL),
            0.0,
            nthreads
        );

        scran_tests::compare_almost_equal_containers(ref_means, sr_means, {});
        scran_tests::compare_almost_equal_containers(ref_vars, sr_vars, {});
        scran_tests::compare_almost_equal_containers(ref_detected, sr_detected, {});
    }

    {
        std::vector<double> dc_means(full_size), dc_vars(full_size), dc_detected(full_size);
        scran_markers::internal::scan_matrix_by_column(
            *dense_column,
            groupings.data(),
            ngroups,
            gsizes,
            dc_means,
            dc_vars,
            dc_detected,
            nthreads
        );

        scran_tests::compare_almost_equal_containers(ref_means, dc_means, {});
        scran_tests::compare_almost_equal_containers(ref_vars, dc_vars, {});
        scran_tests::compare_almost_equal_containers(ref_detected, dc_detected, {});
    }

    {
        std::vector<double> sc_means(full_size), sc_vars(full_size), sc_detected(full_size);
        scran_markers::internal::scan_matrix_by_column(
            *sparse_column,
            groupings.data(),
            ngroups,
            gsizes,
            sc_means,
            sc_vars,
            sc_detected,
            nthreads
        );

        scran_tests::compare_almost_equal_containers(ref_means, sc_means, {});
        scran_tests::compare_almost_equal_containers(ref_vars, sc_vars, {});
        scran_tests::compare_almost_equal_containers(ref_detected, sc_detected, {});
    }
}

TEST_P(ScanMatrixSimpleTest, InterleavedGroups) {
    const auto param = GetParam();
    const auto ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    // For interleaved groups, each group should be represented in each thread.
    auto groupings = create_interleaved_factor(dense_row->ncol(), ngroups);

    const auto full_size = ngroups * dense_row->nrow();
    std::vector<double> ref_means(full_size), ref_vars(full_size), ref_detected(full_size);
    simple_reference(*dense_row, groupings.data(), ngroups, ref_means, ref_vars, ref_detected);

    quick_compare_to_all_implementations(dense_row, dense_column, sparse_row, sparse_column, groupings, ngroups, nthreads, ref_means, ref_vars, ref_detected);
}

TEST_P(ScanMatrixSimpleTest, ContiguousGroups) {
    const auto param = GetParam();
    const auto ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    // Create contiguous groups where cells from the same group are in a single block.
    // This checks that multi-threaded column-wise scans correctly handle situations where a thread has a frequency of zero for a group.
    auto groupings = create_contiguous_factor(dense_row->ncol(), ngroups);
    std::reverse(groupings.begin(), groupings.end()); // reversing for some variety, so that the group with the highest index shows up first.  

    const auto full_size = ngroups * dense_row->nrow();
    std::vector<double> ref_means(full_size), ref_vars(full_size), ref_detected(full_size);
    simple_reference(*dense_row, groupings.data(), ngroups, ref_means, ref_vars, ref_detected);

    quick_compare_to_all_implementations(dense_row, dense_column, sparse_row, sparse_column, groupings, ngroups, nthreads, ref_means, ref_vars, ref_detected);
}

TEST_P(ScanMatrixSimpleTest, EmptyGroups) {
    const auto param = GetParam();
    const auto raw_ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);
    
    // Every odd-indexed group is empty.
    // This checks that all statistics are correctly set to NaN.
    auto groupings = create_interleaved_factor(dense_row->ncol(), raw_ngroups);
    for (auto& g : groupings) {
        g = 2 * g + 1;
    }
    const auto ngroups = raw_ngroups * 2 + 1;

    const int ngenes = dense_row->nrow();
    const auto full_size = ngroups * ngenes;
    std::vector<double> ref_means(full_size), ref_vars(full_size), ref_detected(full_size);
    simple_reference(*dense_row, groupings.data(), ngroups, ref_means, ref_vars, ref_detected);

    for (int r = 0; r < ngenes; ++r) {
        for (int g = 0; g < ngroups; ++g) {
            const auto pos = r * ngroups + g;
            bool is_missing = g % 2 == 0;
            EXPECT_EQ(std::isnan(ref_means[pos]), is_missing);
            EXPECT_EQ(std::isnan(ref_vars[pos]), is_missing);
            EXPECT_EQ(std::isnan(ref_detected[pos]), is_missing);
        }
    }

    quick_compare_to_all_implementations(dense_row, dense_column, sparse_row, sparse_column, groupings, ngroups, nthreads, ref_means, ref_vars, ref_detected);
}

TEST_P(ScanMatrixSimpleTest, OneCellGroups) {
    const auto param = GetParam();
    const auto raw_ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    // First and last groups only have one cell.
    // This checks for correct setting of variances to NaNs.
    auto groupings = create_contiguous_factor(dense_row->ncol(), raw_ngroups);
    for (auto& g : groupings) {
        ++g;
    }
    groupings.front() = raw_ngroups + 1;
    groupings.back() = 0;
    const auto ngroups = raw_ngroups + 2;

    const int ngenes = dense_row->nrow();
    const auto full_size = ngroups * ngenes;
    std::vector<double> ref_means(full_size), ref_vars(full_size), ref_detected(full_size);
    simple_reference(*dense_row, groupings.data(), ngroups, ref_means, ref_vars, ref_detected);

    for (int r = 0; r < ngenes; ++r) {
        for (int g = 0; g < ngroups; ++g) {
            const auto pos = r * ngroups + g;
            EXPECT_FALSE(std::isnan(ref_means[pos]));
            EXPECT_EQ(std::isnan(ref_vars[pos]), (g == 0 || g == ngroups - 1));
            EXPECT_FALSE(std::isnan(ref_detected[pos]));
        }
    }

    quick_compare_to_all_implementations(dense_row, dense_column, sparse_row, sparse_column, groupings, ngroups, nthreads, ref_means, ref_vars, ref_detected);
}

TEST_P(ScanMatrixSimpleTest, VarianceOnly) {
    const auto param = GetParam();
    const auto ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    auto groupings = create_interleaved_factor(dense_row->ncol(), ngroups);
    std::mt19937_64 rng(ngroups + nthreads * 100);
    std::shuffle(groupings.begin(), groupings.end(), rng); // shuffling for some variety.
    auto gsizes = scran_markers::tabulate_groups(dense_row->ncol(), groupings.data(), ngroups);

    const auto full_size = ngroups * dense_row->nrow();
    std::vector<double> ref_means(full_size), ref_vars(full_size), ref_detected(full_size);
    simple_reference(*dense_row, groupings.data(), ngroups, ref_means, ref_vars, ref_detected);

    // Here, we only compute the variance (and the mean), to check that no attempt is made to fill the other statistics (i.e., detected).
    auto quick_compare_row = [&](const auto& mat) -> void {
        std::vector<double> out_means(full_size), out_vars(full_size), out_detected;
        scran_markers::internal::scan_matrix_by_row_full_auc<true>(
            mat,
            groupings.data(),
            ngroups,
            static_cast<int*>(NULL),
            1,
            static_cast<int*>(NULL),
            ngroups,
            gsizes,
            scran_markers::internal::BlockAverageInfo<double>(),
            out_means,
            out_vars,
            out_detected,
            static_cast<double*>(NULL),
            0.0,
            nthreads
        );

        scran_tests::compare_almost_equal_containers(ref_means, out_means, {});
        scran_tests::compare_almost_equal_containers(ref_vars, out_vars, {});
        EXPECT_TRUE(out_detected.empty());
    };

    quick_compare_row(*dense_row);
    quick_compare_row(*sparse_row);

    auto quick_compare_col = [&](const auto& mat) -> void {
        std::vector<double> out_means(full_size), out_vars(full_size), out_detected;
        scran_markers::internal::scan_matrix_by_column(
            mat,
            groupings.data(),
            ngroups,
            gsizes,
            out_means,
            out_vars,
            out_detected,
            nthreads
        );
        scran_tests::compare_almost_equal_containers(ref_means, out_means, {});
        scran_tests::compare_almost_equal_containers(ref_vars, out_vars, {});
        EXPECT_TRUE(out_detected.empty());
    };

    quick_compare_col(*dense_column);
    quick_compare_col(*sparse_column);
}

TEST_P(ScanMatrixSimpleTest, MeanOnly) {
    const auto param = GetParam();
    const auto ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    auto groupings = create_contiguous_factor(dense_row->ncol(), ngroups);
    std::reverse(groupings.begin(), groupings.end()); // reversing for some variety.
    auto gsizes = scran_markers::tabulate_groups(dense_row->ncol(), groupings.data(), ngroups);

    const auto full_size = ngroups * dense_row->nrow();
    std::vector<double> ref_means(full_size), ref_vars(full_size), ref_detected(full_size);
    simple_reference(*dense_row, groupings.data(), ngroups, ref_means, ref_vars, ref_detected);

    // Here, we only compute the mean, to check that no attempt is made to fill the other statistics (i.e., variance, detected).
    auto quick_compare_row = [&](const auto& mat) -> void {
        std::vector<double> out_means(full_size), out_vars, out_detected;
        scran_markers::internal::scan_matrix_by_row_full_auc<true>(
            mat,
            groupings.data(),
            ngroups,
            static_cast<int*>(NULL),
            1,
            static_cast<int*>(NULL),
            ngroups,
            gsizes,
            scran_markers::internal::BlockAverageInfo<double>(),
            out_means,
            out_vars,
            out_detected,
            static_cast<double*>(NULL),
            0.0,
            nthreads
        );
        scran_tests::compare_almost_equal_containers(ref_means, out_means, {});
        EXPECT_TRUE(out_vars.empty());
        EXPECT_TRUE(out_detected.empty());
    };

    quick_compare_row(*dense_row);
    quick_compare_row(*sparse_row);

    auto quick_compare_col = [&](const auto& mat) -> void {
        std::vector<double> out_means(full_size), out_vars, out_detected;
        scran_markers::internal::scan_matrix_by_column(
            mat,
            groupings.data(),
            ngroups,
            gsizes,
            out_means,
            out_vars,
            out_detected,
            nthreads
        );
        scran_tests::compare_almost_equal_containers(ref_means, out_means, {});
        EXPECT_TRUE(out_vars.empty());
        EXPECT_TRUE(out_detected.empty());
    };

    quick_compare_col(*dense_column);
    quick_compare_col(*sparse_column);
}

TEST_P(ScanMatrixSimpleTest, MeanOnlyEmptyGroup) {
    const auto param = GetParam();
    const auto raw_ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    // First and last groups are empty.
    // The previous empty group tests only cover the empty-handling code in the variance calculation section.
    // They won't properly hit the mean-only code, hence our need for an explicit test suite.
    auto groupings = create_interleaved_factor(dense_row->ncol(), raw_ngroups);
    for (auto& g : groupings) {
        ++g;
    }
    const auto ngroups = raw_ngroups + 2;
    auto gsizes = scran_markers::tabulate_groups(dense_row->ncol(), groupings.data(), ngroups);

    const auto ngenes = dense_row->nrow();
    const auto full_size = ngroups * ngenes;
    std::vector<double> ref_means(full_size), ref_vars(full_size), ref_detected(full_size);
    simple_reference(*dense_row, groupings.data(), ngroups, ref_means, ref_vars, ref_detected);

    for (int r = 0; r < ngenes; ++r) {
        for (int g = 0; g < ngroups; ++g) {
            const auto pos = r * ngroups + g;
            EXPECT_EQ(std::isnan(ref_means[pos]), (g == 0 || g == ngroups - 1));
        }
    }

    auto quick_compare_row = [&](const auto& mat) -> void {
        std::vector<double> out_means(full_size), out_vars, out_detected;
        scran_markers::internal::scan_matrix_by_row_full_auc<true>(
            mat,
            groupings.data(),
            ngroups,
            static_cast<int*>(NULL),
            1,
            static_cast<int*>(NULL),
            ngroups,
            gsizes,
            scran_markers::internal::BlockAverageInfo<double>(),
            out_means,
            out_vars,
            out_detected,
            static_cast<double*>(NULL),
            0.0,
            nthreads
        );
        scran_tests::compare_almost_equal_containers(ref_means, out_means, {});
        EXPECT_TRUE(out_vars.empty());
        EXPECT_TRUE(out_detected.empty());
    };

    quick_compare_row(*dense_row);
    quick_compare_row(*sparse_row);

    auto quick_compare_col = [&](const auto& mat) -> void {
        std::vector<double> out_means(full_size), out_vars, out_detected;
        scran_markers::internal::scan_matrix_by_column(
            mat,
            groupings.data(),
            ngroups,
            gsizes,
            out_means,
            out_vars,
            out_detected,
            nthreads
        );
        scran_tests::compare_almost_equal_containers(ref_means, out_means, {});
        EXPECT_TRUE(out_vars.empty());
        EXPECT_TRUE(out_detected.empty());
    };

    quick_compare_col(*dense_column);
    quick_compare_col(*sparse_column);
}

TEST_P(ScanMatrixSimpleTest, DetectedOnly) {
    const auto param = GetParam();
    const auto ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    auto groupings = create_interleaved_factor(dense_row->ncol(), ngroups);
    std::mt19937_64 rng(ngroups + nthreads * 10);
    std::shuffle(groupings.begin(), groupings.end(), rng); // shuffling for some variety.
    auto gsizes = scran_markers::tabulate_groups(dense_row->ncol(), groupings.data(), ngroups);

    const auto full_size = ngroups * dense_row->nrow();
    std::vector<double> ref_means(full_size), ref_vars(full_size), ref_detected(full_size);
    simple_reference(*dense_row, groupings.data(), ngroups, ref_means, ref_vars, ref_detected);

    // Here, we only compute the detected, to check that no attempt is made to fill the other statistics (i.e., mean, variance).
    auto quick_compare_row = [&](const auto& mat) -> void {
        std::vector<double> out_means, out_vars, out_detected(full_size);
        scran_markers::internal::scan_matrix_by_row_full_auc<true>(
            mat,
            groupings.data(),
            ngroups,
            static_cast<int*>(NULL),
            1,
            static_cast<int*>(NULL),
            ngroups,
            gsizes,
            scran_markers::internal::BlockAverageInfo<double>(),
            out_means,
            out_vars,
            out_detected,
            static_cast<double*>(NULL),
            0.0,
            nthreads
        );
        EXPECT_TRUE(out_means.empty());
        EXPECT_TRUE(out_vars.empty());
        scran_tests::compare_almost_equal_containers(ref_detected, out_detected, {});
    };

    quick_compare_row(*dense_row);
    quick_compare_row(*sparse_row);

    auto quick_compare_col = [&](const auto& mat) -> void {
        std::vector<double> out_means, out_vars, out_detected(full_size);
        scran_markers::internal::scan_matrix_by_column(
            mat,
            groupings.data(),
            ngroups,
            gsizes,
            out_means,
            out_vars,
            out_detected,
            nthreads
        );
        EXPECT_TRUE(out_means.empty());
        EXPECT_TRUE(out_vars.empty());
        scran_tests::compare_almost_equal_containers(ref_detected, out_detected, {});
    };

    quick_compare_col(*dense_column);
    quick_compare_col(*sparse_column);
}

TEST_P(ScanMatrixSimpleTest, Blocked) {
    const auto param = GetParam();
    const auto ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    auto groupings = create_interleaved_factor(dense_row->ncol(), ngroups);
    const int nblocks = 3;
    auto blocks = create_contiguous_factor(dense_row->ncol(), nblocks);
    auto combo_out = scran_markers::create_combinations(dense_row->ncol(), groupings.data(), ngroups, blocks.data(), nblocks);

    const auto full_size = combo_out.num_combinations * dense_row->nrow();
    std::vector<double> ref_means(full_size), ref_vars(full_size), ref_detected(full_size);
    simple_reference(*dense_row, combo_out.combinations.data(), combo_out.num_combinations, ref_means, ref_vars, ref_detected);

    auto quick_compare_row = [&](const auto& mat) -> void {
        std::vector<double> out_means(full_size), out_vars(full_size), out_detected(full_size);
        scran_markers::internal::scan_matrix_by_row_full_auc<false>(
            mat,
            groupings.data(),
            ngroups,
            blocks.data(),
            nblocks,
            combo_out.combinations.data(),
            combo_out.num_combinations,
            combo_out.frequencies,
            scran_markers::internal::BlockAverageInfo<double>(),
            out_means,
            out_vars,
            out_detected,
            static_cast<double*>(NULL),
            0.0,
            nthreads
        );
        scran_tests::compare_almost_equal_containers(ref_means, out_means, {});
        scran_tests::compare_almost_equal_containers(ref_vars, out_vars, {});
        scran_tests::compare_almost_equal_containers(ref_detected, out_detected, {});
    };

    quick_compare_row(*dense_row);
    quick_compare_row(*sparse_row);

    auto quick_compare_col = [&](const auto& mat) -> void {
        std::vector<double> out_means(full_size), out_vars(full_size), out_detected(full_size);
        scran_markers::internal::scan_matrix_by_column(
            mat,
            combo_out.combinations.data(),
            combo_out.num_combinations,
            combo_out.frequencies,
            out_means,
            out_vars,
            out_detected,
            nthreads
        );
        scran_tests::compare_almost_equal_containers(ref_means, out_means, {});
        scran_tests::compare_almost_equal_containers(ref_vars, out_vars, {});
        scran_tests::compare_almost_equal_containers(ref_detected, out_detected, {});
    };

    quick_compare_col(*dense_column);
    quick_compare_col(*sparse_column);
}

TEST_P(ScanMatrixSimpleTest, BlockedConfounded) {
    const auto param = GetParam();
    const auto ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    // Each group is its own block.
    auto groupings = create_interleaved_factor(dense_row->ncol(), ngroups);
    auto combo_out = scran_markers::create_combinations(dense_row->ncol(), groupings.data(), ngroups, groupings.data(), ngroups);

    const auto ngenes = dense_row->nrow();
    const auto full_size = combo_out.num_combinations * ngenes;
    std::vector<double> ref_means(full_size), ref_vars(full_size), ref_detected(full_size);
    simple_reference(*dense_row, combo_out.combinations.data(), combo_out.num_combinations, ref_means, ref_vars, ref_detected);

    for (int r = 0; r < ngenes; ++r) {
        for (int g = 0; g < ngroups; ++g) {
            for (int b = 0; b < ngroups; ++b) {
                const auto pos = (r * ngroups + b) * ngroups + g;
                const bool is_missing = g != b;
                EXPECT_EQ(std::isnan(ref_means[pos]), is_missing);
                EXPECT_EQ(std::isnan(ref_vars[pos]), is_missing);
                EXPECT_EQ(std::isnan(ref_detected[pos]), is_missing);
            }
        }
    }

    auto quick_compare_row = [&](const auto& mat) -> void {
        std::vector<double> out_means(full_size), out_vars(full_size), out_detected(full_size);
        scran_markers::internal::scan_matrix_by_row_full_auc<false>(
            mat,
            groupings.data(),
            ngroups,
            groupings.data(),
            ngroups,
            combo_out.combinations.data(),
            combo_out.num_combinations,
            combo_out.frequencies,
            scran_markers::internal::BlockAverageInfo<double>(),
            out_means,
            out_vars,
            out_detected,
            static_cast<double*>(NULL),
            0.0,
            nthreads
        );
        scran_tests::compare_almost_equal_containers(ref_means, out_means, {});
        scran_tests::compare_almost_equal_containers(ref_vars, out_vars, {});
        scran_tests::compare_almost_equal_containers(ref_detected, out_detected, {});
    };

    quick_compare_row(*dense_row);
    quick_compare_row(*sparse_row);

    auto quick_compare_col = [&](const auto& mat) -> void {
        std::vector<double> out_means(full_size), out_vars(full_size), out_detected(full_size);
        scran_markers::internal::scan_matrix_by_column(
            mat,
            combo_out.combinations.data(),
            combo_out.num_combinations,
            combo_out.frequencies,
            out_means,
            out_vars,
            out_detected,
            nthreads
        );
        scran_tests::compare_almost_equal_containers(ref_means, out_means, {});
        scran_tests::compare_almost_equal_containers(ref_vars, out_vars, {});
        scran_tests::compare_almost_equal_containers(ref_detected, out_detected, {});
    };

    quick_compare_col(*dense_column);
    quick_compare_col(*sparse_column);
}

INSTANTIATE_TEST_SUITE_P(
    ScanMatrix,
    ScanMatrixSimpleTest,
    ::testing::Combine(
        ::testing::Values(2, 3, 4, 5), // number of groups
        ::testing::Values(1, 3) // number of threads
    )
);

/*********************************************/

class ScanMatrixAucTest : public ::testing::TestWithParam<std::tuple<int, int> > {
protected:
    inline static std::shared_ptr<tatami::Matrix<double, int> > dense_row, sparse_row;

    static void SetUpTestSuite() {
        size_t nr = 598, nc = 376;
        dense_row.reset(
            new tatami::DenseRowMatrix<double, int>(
                nr,
                nc,
                scran_tests::simulate_vector(
                    nr * nc, 
                    []{
                        scran_tests::SimulateVectorParameters sparam;
                        sparam.density = 0.25;
                        sparam.seed = 453298;
                        return sparam;
                    }()
                )
            )
        );

        sparse_row = tatami::convert_to_compressed_sparse(dense_row.get(), true);
    }
};

TEST_P(ScanMatrixAucTest, BasicConsistency) {
    const auto param = GetParam();
    const auto ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    auto groupings = create_interleaved_factor(dense_row->ncol(), ngroups);
    auto gsizes = scran_markers::tabulate_groups(dense_row->ncol(), groupings.data(), ngroups);

    const auto ngenes = dense_row->nrow();
    const auto full_size = ngroups * ngroups * ngenes;
    std::vector<double> mock_means, mock_vars, mock_detected;

    auto quick_compute = [&](const auto& mat, scran_markers::internal::BlockAverageInfo<double>& ave_info) -> std::vector<double>  {
        std::vector<double> out_auc(full_size);
        scran_markers::internal::scan_matrix_by_row_full_auc<true>(
            mat,
            groupings.data(),
            ngroups,
            static_cast<int*>(NULL),
            1,
            static_cast<int*>(NULL),
            ngroups,
            gsizes,
            ave_info,
            mock_means,
            mock_vars,
            mock_detected,
            out_auc.data(),
            0.0,
            nthreads
        );
        return out_auc;
    };

    scran_markers::internal::BlockAverageInfo<double> mean_info(
        scran_blocks::compute_weights<double>(
            gsizes,
            scran_blocks::WeightPolicy::EQUAL,
            {}
        )
    );
    auto ref_auc = quick_compute(*dense_row, mean_info);

    for (int r = 0; r < ngenes; ++r) {
        for (int g1 = 0; g1 < ngroups; ++g1) {
            for (int g2 = 0; g2 < g1; ++g2) {
                const auto val = ref_auc[(r * ngroups + g1) * ngroups + g2];
                EXPECT_FALSE(std::isnan(val));
                const auto revval = ref_auc[(r * ngroups + g2) * ngroups + g1];
                scran_tests::compare_almost_equal(1 - revval, val, scran_tests::CompareAlmostEqualParameters());
            }
            EXPECT_EQ(ref_auc[(r * ngroups + g1) * ngroups + g1], 0);
        }
    }

    // Checking that dense and sparse results are the same.
    auto sparse_auc = quick_compute(*sparse_row, mean_info);
    EXPECT_EQ(sparse_auc, ref_auc);

    // Checking that quantile and mean results are the same when there's only one block.
    scran_markers::internal::BlockAverageInfo<double> quantile_info(0.5);
    auto quantile_auc = quick_compute(*dense_row, quantile_info);
    EXPECT_EQ(quantile_auc, ref_auc);

    auto sparse_quantile_auc = quick_compute(*sparse_row, quantile_info);
    EXPECT_EQ(sparse_quantile_auc, ref_auc);
}

TEST_P(ScanMatrixAucTest, EmptyGroups) {
    const auto param = GetParam();
    const auto raw_ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    auto raw_groupings = create_contiguous_factor(dense_row->ncol(), raw_ngroups);
    auto raw_gsizes = scran_markers::tabulate_groups(dense_row->ncol(), raw_groupings.data(), raw_ngroups);

    const auto ngenes = dense_row->nrow();
    std::vector<double> mock_means, mock_vars, mock_detected;
    std::vector<double> raw_ref_auc(raw_ngroups * raw_ngroups * ngenes);

    scran_markers::internal::BlockAverageInfo<double> raw_mean_info(
        scran_blocks::compute_weights<double>(
            raw_gsizes,
            scran_blocks::WeightPolicy::EQUAL,
            {}
        )
    );
    
    scran_markers::internal::scan_matrix_by_row_full_auc<true>(
        *dense_row,
        raw_groupings.data(),
        raw_ngroups,
        static_cast<int*>(NULL),
        1,
        static_cast<int*>(NULL),
        raw_ngroups,
        raw_gsizes,
        raw_mean_info,
        mock_means,
        mock_vars,
        mock_detected,
        raw_ref_auc.data(),
        0.0,
        nthreads
    );

    // Inject empty groups, namely at 0 and and 'ngroups - 1'.
    auto groupings = raw_groupings;
    for (auto& g : groupings) {
        ++g;
    }
    auto gsizes = raw_gsizes;
    gsizes.insert(gsizes.begin(), 0);
    gsizes.push_back(0);
    const auto ngroups = raw_ngroups + 2;

    const auto full_size = ngroups * ngroups * ngenes;
    std::vector<double> ref_auc(full_size, std::numeric_limits<double>::quiet_NaN());
    for (int r = 0; r < ngenes; ++r) {
        for (int g1 = 0; g1 < raw_ngroups; ++g1) {
            for (int g2 = 0; g2 < raw_ngroups; ++g2) {
                ref_auc[(r * ngroups + g1 + 1) * ngroups + g2 + 1] = raw_ref_auc[(r * raw_ngroups + g1) * raw_ngroups + g2];
            }
        }
        ref_auc[r * ngroups * ngroups] = 0;
        ref_auc[(r + 1) * ngroups * ngroups - 1] = 0;
    }

    // Checking what happens when we leave the empty groups in there.
    auto quick_compare = [&](const auto& mat, scran_markers::internal::BlockAverageInfo<double>& ave_info) -> void {
        std::vector<double> output_auc(full_size);
        scran_markers::internal::scan_matrix_by_row_full_auc<true>(
            mat,
            groupings.data(),
            ngroups,
            static_cast<int*>(NULL),
            1,
            static_cast<int*>(NULL),
            ngroups,
            gsizes,
            ave_info,
            mock_means,
            mock_vars,
            mock_detected,
            output_auc.data(),
            0.0,
            nthreads
        );
        scran_tests::compare_almost_equal_containers(ref_auc, output_auc, {});
    };

    scran_markers::internal::BlockAverageInfo<double> mean_info(
        scran_blocks::compute_weights<double>(
            gsizes,
            scran_blocks::WeightPolicy::EQUAL,
            {}
        )
    );
    quick_compare(*dense_row, mean_info);
    quick_compare(*sparse_row, mean_info);

    // Also testing the quantile code. The result should be the same as we only have one block.
    scran_markers::internal::BlockAverageInfo<double> quantile_info(0.5);
    quick_compare(*dense_row, quantile_info);
    quick_compare(*sparse_row, quantile_info);
}

TEST_P(ScanMatrixAucTest, Threshold) {
    const auto param = GetParam();
    const auto ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    auto groupings = create_contiguous_factor(dense_row->ncol(), ngroups);
    auto gsizes = scran_markers::tabulate_groups(dense_row->ncol(), groupings.data(), ngroups);

    const auto ngenes = dense_row->nrow();
    const auto full_size = ngroups * ngroups * ngenes;
    std::vector<double> mock_means, mock_vars, mock_detected;
    std::vector<double> ref_auc(full_size);

    scran_markers::internal::BlockAverageInfo<double> mean_info(
        scran_blocks::compute_weights<double>(
            gsizes,
            scran_blocks::WeightPolicy::EQUAL,
            {}
        )
    );

    scran_markers::internal::scan_matrix_by_row_full_auc<true>(
        *dense_row,
        groupings.data(),
        ngroups,
        static_cast<int*>(NULL),
        1,
        static_cast<int*>(NULL),
        ngroups,
        gsizes,
        mean_info,
        mock_means,
        mock_vars,
        mock_detected,
        ref_auc.data(),
        0.,
        nthreads
    );

    auto quick_compute = [&](const auto& mat, scran_markers::internal::BlockAverageInfo<double>& ave_info) -> std::vector<double> {
        std::vector<double> output_auc(full_size);
        scran_markers::internal::scan_matrix_by_row_full_auc<true>(
            mat,
            groupings.data(),
            ngroups,
            static_cast<int*>(NULL),
            1,
            static_cast<int*>(NULL),
            ngroups,
            gsizes,
            ave_info,
            mock_means,
            mock_vars,
            mock_detected,
            output_auc.data(),
            0.5,
            nthreads
        );
        return output_auc;
    };

    // Check that thresholding has some effect and that the AUCs are pushed towards zero.
    auto threshold_auc = quick_compute(*dense_row, mean_info);

    EXPECT_NE(threshold_auc, ref_auc);
    for (int r = 0; r < ngenes; ++r) {
        for (int g1 = 0; g1 < ngroups; ++g1) {
            for (int g2 = 0; g2 < ngroups; ++g2) {
                auto thresholded = threshold_auc[(r * ngroups + g1) * ngroups + g2];
                if (g1 == g2) {
                    EXPECT_EQ(thresholded, 0);
                } else {
                    auto original = ref_auc[(r * ngroups + g1) * ngroups + g2];
                    EXPECT_GE(original, thresholded);
                }
            }
        }
    }

    auto sparse_threshold_auc = quick_compute(*sparse_row, mean_info);
    EXPECT_EQ(threshold_auc, sparse_threshold_auc);

    // Checking that quantile and mean results are the same when there's only one block.
    scran_markers::internal::BlockAverageInfo<double> quantile_info(0.5);

    auto quantile_threshold_auc = quick_compute(*dense_row, quantile_info);
    EXPECT_EQ(threshold_auc, quantile_threshold_auc);

    auto sparse_quantile_threshold_auc = quick_compute(*dense_row, quantile_info);
    EXPECT_EQ(threshold_auc, sparse_quantile_threshold_auc);
}

static std::vector<std::vector<double> > compute_individual_aucs(
    const tatami::Matrix<double, int>& mat,
    const int* grouping,
    const int num_groups,
    const int* blocks,
    const int num_blocks
) {
    std::vector<std::vector<double> > output(num_blocks);
    const auto ngenes = mat.nrow();
    const auto ncells = mat.ncol();
    const auto full_size = num_groups * num_groups * ngenes;

    std::vector<double> mock_means, mock_vars, mock_detected;
    for (int b = 0; b < num_blocks; ++b) {
        std::vector<int> subset, subgroup;
        for (int c = 0; c < ncells; ++c) {
            if (blocks[c] == b) {
                subset.push_back(c);
                subgroup.push_back(grouping[c]);
            }
        }

        auto subgsizes = scran_markers::tabulate_groups<int, int>(subset.size(), subgroup.data(), num_groups);
        scran_markers::internal::BlockAverageInfo<double> mean_info(
            scran_blocks::compute_weights<double>(
                subgsizes,
                scran_blocks::WeightPolicy::EQUAL,
                {}
            )
        );

        output[b].resize(full_size);
        auto sub = tatami::make_DelayedSubset<double, int>(tatami::wrap_shared_ptr(&mat), std::move(subset), false);
        scran_markers::internal::scan_matrix_by_row_full_auc<true>(
            *sub,
            subgroup.data(),
            num_groups,
            static_cast<int*>(NULL),
            1,
            static_cast<int*>(NULL),
            num_groups,
            subgsizes,
            mean_info,
            mock_means,
            mock_vars,
            mock_detected,
            output[b].data(),
            0.,
            1 
        );
    }

    return output;
}

TEST_P(ScanMatrixAucTest, BlockedMean) {
    const auto param = GetParam();
    const auto ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    auto groupings = create_interleaved_factor(dense_row->ncol(), ngroups);
    const int nblocks = 3;
    auto blocks = create_contiguous_factor(dense_row->ncol(), nblocks);
    auto combo_out = scran_markers::create_combinations(dense_row->ncol(), groupings.data(), ngroups, blocks.data(), nblocks);

    auto per_block_aucs = compute_individual_aucs(*dense_row, groupings.data(), ngroups, blocks.data(), nblocks);

    const auto ngenes = dense_row->nrow();
    const int ngroups2 = ngroups * ngroups;
    const std::size_t full_size = ngroups2 * ngenes;
    std::vector<double> mock_means, mock_vars, mock_detected;

    auto quick_compute = [&](const auto& mat, scran_markers::internal::BlockAverageInfo<double>& ave_info) -> std::vector<double> {
        std::vector<double> output_auc(full_size);
        scran_markers::internal::scan_matrix_by_row_full_auc<false>(
            mat,
            groupings.data(),
            ngroups,
            blocks.data(),
            nblocks,
            combo_out.combinations.data(),
            combo_out.num_combinations,
            combo_out.frequencies,
            ave_info,
            mock_means,
            mock_vars,
            mock_detected,
            output_auc.data(),
            0.,
            nthreads
        );
        return output_auc;
    };

    // Using equally weighted blocks.
    {
        scran_markers::internal::BlockAverageInfo<double> mean_info(
            scran_blocks::compute_weights<double>(
                combo_out.frequencies,
                scran_blocks::WeightPolicy::EQUAL,
                {}
            )
        );
        auto block_auc = quick_compute(*dense_row, mean_info);

        std::vector<double> expected_auc(per_block_aucs.front());
        for (int b = 1; b < nblocks; ++b) {
            for (std::size_t f = 0; f < full_size; ++f) {
                expected_auc[f] += per_block_aucs[b][f];
            }
        }
        for (std::size_t f = 0; f < full_size; ++f) {
            expected_auc[f] /= nblocks;
        }
        scran_tests::compare_almost_equal_containers(expected_auc, block_auc, {});
    }

    // Using blocks weighted by the product of the group sizes.
    {
        scran_markers::internal::BlockAverageInfo<double> mean_info(
            scran_blocks::compute_weights<double>(
                combo_out.frequencies,
                scran_blocks::WeightPolicy::SIZE,
                {}
            )
        );
        auto block_auc = quick_compute(*sparse_row, mean_info); // using sparse matrix for some variety.

        // Calculation of the expected AUC is a bit tricky, as we need to convert the AUC back into its pre-scaled value.
        // Doing so requires multiplication of the denominators, then addition, then division by the summed denominators.
        std::vector<double> expected_auc(full_size);
        std::vector<double> actual_total(ngroups2); 
        for (int b = 0; b < nblocks; ++b) {
            for (int r = 0; r < ngenes; ++r) {
                for (int g1 = 0; g1 < ngroups; ++g1) {
                    const double w1 = combo_out.frequencies[b * ngroups + g1];
                    for (int g2 = 0; g2 < ngroups; ++g2) {
                        const double w2 = combo_out.frequencies[b * ngroups + g2];
                        const auto pos = (r * ngroups + g1) * ngroups + g2;
                        expected_auc[pos] += per_block_aucs[b][pos] * w1 * w2;
                    }
                }
            }
            for (int g1 = 0; g1 < ngroups; ++g1) {
                const double w1 = combo_out.frequencies[b * ngroups + g1];
                for (int g2 = 0; g2 < ngroups; ++g2) {
                    const double w2 = combo_out.frequencies[b * ngroups + g2];
                    actual_total[g1 * ngroups + g2] += w1 * w2;
                }
            }
        }
        for (int r = 0; r < ngenes; ++r) {
            for (int h = 0; h < ngroups2; ++h) {
                expected_auc[r * ngroups2 + h] /= actual_total[h];
            }
        }
        scran_tests::compare_almost_equal_containers(expected_auc, block_auc, {});
    }
}

TEST_P(ScanMatrixAucTest, BlockedQuantile) {
    const auto param = GetParam();
    const auto ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    auto groupings = create_interleaved_factor(dense_row->ncol(), ngroups);
    std::reverse(groupings.begin(), groupings.end()); // reversing for some variety.
    const int nblocks = 3;
    auto blocks = create_contiguous_factor(dense_row->ncol(), nblocks);
    std::reverse(blocks.begin(), blocks.end()); // reversing for some variety.
    auto combo_out = scran_markers::create_combinations(dense_row->ncol(), groupings.data(), ngroups, blocks.data(), nblocks);

    auto per_block_aucs = compute_individual_aucs(*dense_row, groupings.data(), ngroups, blocks.data(), nblocks);

    const auto ngenes = dense_row->nrow();
    const int ngroups2 = ngroups * ngroups;
    const std::size_t full_size = ngroups2 * ngenes;
    std::vector<double> mock_means, mock_vars, mock_detected;

    auto quick_compute = [&](const auto& mat, scran_markers::internal::BlockAverageInfo<double>& ave_info) -> std::vector<double> {
        std::vector<double> output_auc(full_size);
        scran_markers::internal::scan_matrix_by_row_full_auc<false>(
            mat,
            groupings.data(),
            ngroups,
            blocks.data(),
            nblocks,
            combo_out.combinations.data(),
            combo_out.num_combinations,
            combo_out.frequencies,
            ave_info,
            mock_means,
            mock_vars,
            mock_detected,
            output_auc.data(),
            0.,
            nthreads
        );
        return output_auc;
    };

    // Using the median.
    {
        scran_markers::internal::BlockAverageInfo<double> quantile_info(0.5);
        auto block_auc = quick_compute(*dense_row, quantile_info);

        std::vector<double> expected_auc(per_block_aucs.front());
        std::vector<double> buffer;
        for (std::size_t f = 0; f < full_size; ++f) {
            buffer.clear();
            for (int b = 0; b < nblocks; ++b) {
                buffer.push_back(per_block_aucs[b][f]);
            }
            expected_auc[f] = quickstats::median(buffer.size(), buffer.data());
        }
        scran_tests::compare_almost_equal_containers(expected_auc, block_auc, {});
    }

    // Using the minimum.
    {
        scran_markers::internal::BlockAverageInfo<double> quantile_info(0);
        auto block_auc = quick_compute(*sparse_row, quantile_info); // using sparse matrix for some variety.

        std::vector<double> expected_auc(per_block_aucs.front());
        std::vector<double> buffer;
        for (std::size_t f = 0; f < full_size; ++f) {
            buffer.clear();
            for (int b = 0; b < nblocks; ++b) {
                buffer.push_back(per_block_aucs[b][f]);
            }
            expected_auc[f] = *std::min_element(buffer.begin(), buffer.end());
        }
        scran_tests::compare_almost_equal_containers(expected_auc, block_auc, {});
    }
}

TEST_P(ScanMatrixAucTest, BlockedOverlap) {
    const auto param = GetParam();
    const auto raw_ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    // Here, we're interested in a situation where some of the groups overlap between the blocks,
    // while some of the groups are unique to each block.
    const auto ncells = dense_row->ncol();
    const int midpt = ncells / 2;
    auto raw_groupings1 = create_interleaved_factor(midpt, raw_ngroups);
    auto raw_groupings2 = create_contiguous_factor(ncells - midpt, raw_ngroups);

    const int shift = raw_ngroups / 2;
    auto groupings = raw_groupings1;
    for (auto& g : groupings) {
        g += shift;
    }
    groupings.insert(groupings.end(), raw_groupings2.begin(), raw_groupings2.end());
    const int ngroups = raw_ngroups + shift;

    // Specifically, the last 'shift' groups are unique to the first block while the first 'shift' groups are unique to the second block.
    const int nblocks = 2;
    std::vector<int> blocks(ncells);
    std::fill(blocks.begin() + midpt, blocks.end(), 1);

    auto per_block_aucs = compute_individual_aucs(*dense_row, groupings.data(), ngroups, blocks.data(), nblocks);
    const auto ngenes = dense_row->nrow(); 
    for (int r = 0; r < ngenes; ++r) {
        for (int g1 = 0; g1 < ngroups; ++g1) {
            for (int g2 = 0; g2 < ngroups; ++g2) {
                const auto pos = (r * ngroups + g1) * ngroups + g2;
                if (g1 == g2) {
                    EXPECT_EQ(per_block_aucs[0][pos], 0);
                    EXPECT_EQ(per_block_aucs[1][pos], 0);
                } else {
                    EXPECT_EQ(std::isnan(per_block_aucs[0][pos]), g1 < shift || g2 < shift);
                    EXPECT_EQ(std::isnan(per_block_aucs[1][pos]), g1 >= raw_ngroups || g2 >= raw_ngroups);
                }
            }
        }
    }

    const std::size_t full_size = ngroups * ngroups * ngenes;
    std::vector<double> mock_means, mock_vars, mock_detected;
    auto combo_out = scran_markers::create_combinations(dense_row->ncol(), groupings.data(), ngroups, blocks.data(), nblocks);

    auto quick_compute = [&](const auto& mat, scran_markers::internal::BlockAverageInfo<double>& ave_info) -> std::vector<double> {
        std::vector<double> output_auc(full_size);
        scran_markers::internal::scan_matrix_by_row_full_auc<false>(
            mat,
            groupings.data(),
            ngroups,
            blocks.data(),
            nblocks,
            combo_out.combinations.data(),
            combo_out.num_combinations,
            combo_out.frequencies,
            ave_info,
            mock_means,
            mock_vars,
            mock_detected,
            output_auc.data(),
            0.,
            nthreads
        );
        return output_auc;
    };

    // Comparing with block means.
    {
        scran_markers::internal::BlockAverageInfo<double> mean_info(
            scran_blocks::compute_weights<double>(
                combo_out.frequencies,
                scran_blocks::WeightPolicy::EQUAL,
                {}
            )
        );
        auto mean_auc = quick_compute(*sparse_row, mean_info); // using a sparse matrix for some variety.

        std::vector<double> expected_auc(full_size);
        for (std::size_t f = 0; f < full_size; ++f) {
            const auto val1 = per_block_aucs[0][f];
            const auto val2 = per_block_aucs[1][f];
            expected_auc[f] = ((std::isnan(val1) ? 0 : val1) + (std::isnan(val2) ? 0 : val2)) / (!std::isnan(val1) + !std::isnan(val2));
        }
        scran_tests::compare_almost_equal_containers(expected_auc, mean_auc, {});
    }

    // Comparing with block quantiles.
    {
        scran_markers::internal::BlockAverageInfo<double> quantile_info(1);
        auto quantile_auc = quick_compute(*dense_row, quantile_info);

        std::vector<double> expected_auc(full_size);
        for (std::size_t f = 0; f < full_size; ++f) {
            const auto val1 = per_block_aucs[0][f];
            const auto val2 = per_block_aucs[1][f];
            if (std::isnan(val1)) {
                expected_auc[f] = val2;
            } else if (std::isnan(val2)) {
                expected_auc[f] = val1;
            } else {
                expected_auc[f] = std::max(val1, val2);
            }
        }
        scran_tests::compare_almost_equal_containers(expected_auc, quantile_auc, {});
    }
}

TEST_P(ScanMatrixAucTest, BlockedConfounded) {
    const auto param = GetParam();
    const auto ngroups = std::get<0>(param);
    const auto nthreads = std::get<1>(param);

    // Here, we're interested in a situation where the groups are completely confounded with the blocks.
    // while some of the groups are unique to each block.
    const auto ncells = dense_row->ncol();
    auto groupings = create_interleaved_factor(ncells, ngroups);
    std::mt19937_64 rng(342 + nthreads * ngroups); // shuffling for some variety.
    std::shuffle(groupings.begin(), groupings.end(), rng);

    const int nblocks = 2;
    std::vector<int> blocks;
    blocks.reserve(ncells);
    for (auto g : groupings) {
        blocks.push_back(g % 2 == 0); // every even group goes in the second block.
    }

    auto per_block_aucs = compute_individual_aucs(*dense_row, groupings.data(), ngroups, blocks.data(), nblocks);
    const auto ngenes = dense_row->nrow(); 
    for (int r = 0; r < ngenes; ++r) {
        for (int g1 = 0; g1 < ngroups; ++g1) {
            for (int g2 = 0; g2 < ngroups; ++g2) {
                const auto pos = (r * ngroups + g1) * ngroups + g2;
                if (g1 == g2) {
                    EXPECT_EQ(per_block_aucs[0][pos], 0);
                    EXPECT_EQ(per_block_aucs[1][pos], 0);
                } else {
                    EXPECT_EQ(std::isnan(per_block_aucs[0][pos]), g1 % 2 == 0 || g2 % 2 == 0);
                    EXPECT_EQ(std::isnan(per_block_aucs[1][pos]), g1 % 2 == 1 || g2 % 2 == 1);
                }
            }
        }
    }

    const std::size_t full_size = ngroups * ngroups * ngenes;
    std::vector<double> mock_means, mock_vars, mock_detected;
    auto combo_out = scran_markers::create_combinations(dense_row->ncol(), groupings.data(), ngroups, blocks.data(), nblocks);

    auto quick_compute = [&](const auto& mat, scran_markers::internal::BlockAverageInfo<double>& ave_info) -> std::vector<double> {
        std::vector<double> output_auc(full_size);
        scran_markers::internal::scan_matrix_by_row_full_auc<false>(
            mat,
            groupings.data(),
            ngroups,
            blocks.data(),
            nblocks,
            combo_out.combinations.data(),
            combo_out.num_combinations,
            combo_out.frequencies,
            ave_info,
            mock_means,
            mock_vars,
            mock_detected,
            output_auc.data(),
            0.,
            nthreads
        );
        return output_auc;
    };

    // Comparing with block means.
    {
        scran_markers::internal::BlockAverageInfo<double> mean_info(
            scran_blocks::compute_weights<double>(
                combo_out.frequencies,
                scran_blocks::WeightPolicy::EQUAL,
                {}
            )
        );
        auto mean_auc = quick_compute(*sparse_row, mean_info); // using a sparse matrix for some variety.

        std::vector<double> expected_auc(full_size);
        for (std::size_t f = 0; f < full_size; ++f) {
            const auto val1 = per_block_aucs[0][f];
            const auto val2 = per_block_aucs[1][f];
            expected_auc[f] = (std::isnan(val1) ? val2 : val1);
        }
        scran_tests::compare_almost_equal_containers(expected_auc, mean_auc, {});
    }

    // Comparing with block quantiles.
    {
        scran_markers::internal::BlockAverageInfo<double> quantile_info(1);
        auto quantile_auc = quick_compute(*dense_row, quantile_info);

        std::vector<double> expected_auc(full_size);
        for (std::size_t f = 0; f < full_size; ++f) {
            const auto val1 = per_block_aucs[0][f];
            const auto val2 = per_block_aucs[1][f];
            expected_auc[f] = (std::isnan(val1) ? val2 : val1);
        }
        scran_tests::compare_almost_equal_containers(expected_auc, quantile_auc, {});
    }
}

INSTANTIATE_TEST_SUITE_P(
    ScanMatrix,
    ScanMatrixAucTest,
    ::testing::Combine(
        ::testing::Values(2, 3, 4, 5), // number of groups
        ::testing::Values(1, 3) // number of threads
    )
);
