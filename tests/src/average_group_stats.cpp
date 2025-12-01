#include "scran_tests/scran_tests.hpp"
#include "scran_blocks/scran_blocks.hpp"

#include "scran_markers/average_group_stats.hpp"

class AverageGroupStatsTest : public ::testing::Test {
protected:
    inline static std::size_t ngroups = 5, nblocks = 3, ngenes = 100; 
    inline static std::vector<double> means, weights;

    static void SetUpTestSuite() {
        means = scran_tests::simulate_vector(ngenes * ngroups * nblocks, []{ 
            scran_tests::SimulationParameters sparam;
            sparam.seed = 43210;
            return sparam;
        }());

        weights = scran_tests::simulate_vector(ngroups * nblocks, []{ 
            scran_tests::SimulationParameters sparam;
            sparam.lower = 0.1;
            sparam.upper = 1;
            sparam.seed = 69;
            return sparam;
        }());
    }
};

TEST_F(AverageGroupStatsTest, Mean) {
    std::vector<std::vector<double> > average_mean(ngroups);
    std::vector<double*> mean_ptrs(ngroups); 
    for (std::size_t g = 0; g < ngroups; ++g) {
        average_mean[g].resize(ngenes);
        mean_ptrs[g] = average_mean[g].data();
    }

    auto total_weights = scran_markers::internal::compute_total_weight_per_group(ngroups, nblocks, weights.data());
    for (std::size_t r = 0; r < ngenes; ++r) {
        size_t offset = r * ngroups * nblocks;
        scran_markers::internal::average_group_stats_blockmean(
            r, 
            ngroups,
            nblocks,
            means.data() + offset,
            weights.data(),
            total_weights.data(),
            mean_ptrs
        );
    }

    // Doing a reference calculation.
    std::vector<std::vector<double> > buffers(nblocks);
    std::vector<double*> ptrs;
    std::vector<double> current_weights(nblocks);
    for (auto& bk : buffers) {
        bk.resize(ngenes);
        ptrs.push_back(bk.data());
    }

    for (std::size_t g = 0; g < ngroups; ++g) {
        for (std::size_t b = 0; b < nblocks; ++b) {
            current_weights[b] = weights[b * ngroups + g];
        }

        for (std::size_t r = 0; r < ngenes; ++r) {
            for (std::size_t b = 0; b < nblocks; ++b) {
                buffers[b][r] = means[(r * nblocks + b) * ngroups + g];
            }
        }
        auto expected_m = scran_blocks::parallel_weighted_means(ngenes, ptrs, current_weights.data(), /* skip_nan = */ false);
        scran_tests::compare_almost_equal(expected_m, average_mean[g]);
    }
}

TEST_F(AverageGroupStatsTest, Quantile) {
    std::vector<std::vector<double> > average_quantile(ngroups);
    std::vector<double*> quantile_ptrs(ngroups); 
    for (size_t g = 0; g < ngroups; ++g) {
        average_quantile[g].resize(ngenes);
        quantile_ptrs[g] = average_quantile[g].data();
    }

    std::vector<double> buffer;
    scran_blocks::SingleQuantileVariable<double, typename std::vector<double>::iterator> qcalc(nblocks, 0.5);
    for (std::size_t r = 0; r < ngenes; ++r) {
        std::size_t offset = r * ngroups * nblocks;
        scran_markers::internal::average_group_stats_blockquantile(
            r, 
            ngroups,
            nblocks,
            means.data() + offset,
            buffer,
            qcalc,
            quantile_ptrs
        );
    }

    // Doing a reference calculation.
    std::vector<std::vector<double> > buffers(nblocks);
    std::vector<double*> ptrs;
    std::vector<double> current_weights(nblocks);
    for (auto& bk : buffers) {
        bk.resize(ngenes);
        ptrs.push_back(bk.data());
    }

    for (std::size_t g = 0; g < ngroups; ++g) {
        for (std::size_t b = 0; b < nblocks; ++b) {
            current_weights[b] = weights[b * ngroups + g];
        }

        for (std::size_t r = 0; r < ngenes; ++r) {
            for (std::size_t b = 0; b < nblocks; ++b) {
                buffers[b][r] = means[(r * nblocks + b) * ngroups + g];
            }
        }
        auto expected_q = scran_blocks::parallel_quantiles(ngenes, ptrs, 0.5, /* skip_nan = */ false);
        scran_tests::compare_almost_equal(expected_q, average_quantile[g]);
    }
}

TEST_F(AverageGroupStatsTest, ZeroedBlockMean) {
    // Zeroing all weights in the second block.
    auto copy_weights = weights;
    std::fill_n(copy_weights.begin() + ngroups, ngroups, 0);

    std::vector<std::vector<double> > average_mean(ngroups);
    std::vector<double*> mean_ptrs(ngroups); 
    for (size_t g = 0; g < ngroups; ++g) {
        average_mean[g].resize(ngenes);
        mean_ptrs[g] = average_mean[g].data();
    }

    auto means_copy = means;
    auto total_weights = scran_markers::internal::compute_total_weight_per_group(ngroups, nblocks, copy_weights.data());
    for (std::size_t r = 0; r < ngenes; ++r) {
        std::size_t offset = r * ngroups * nblocks;

        // Also setting all associated means to NaN, to ensure that
        // we explicitly don't use the zero weights for anything.
        auto mcptr = means_copy.data() + offset;
        std::fill_n(mcptr + ngroups, ngroups, std::numeric_limits<double>::quiet_NaN());

        scran_markers::internal::average_group_stats_blockmean(
            r, 
            ngroups,
            nblocks,
            mcptr,
            copy_weights.data(),
            total_weights.data(),
            mean_ptrs
        );
    }

    // Doing a reference calculation, but skipping the middle block.
    std::vector<std::vector<double> > buffers(nblocks - 1);
    std::vector<double*> ptrs;
    std::vector<double> current_weights(nblocks - 1);
    for (auto& bk : buffers) {
        bk.resize(ngenes);
        ptrs.push_back(bk.data());
    }
    std::vector<std::pair<int, int> > block_remap{ { 0, 0 }, { 1, 2 }};

    for (std::size_t g = 0; g < ngroups; ++g) {
        for (const auto& br : block_remap) {
            current_weights[br.first] = copy_weights[br.second * ngroups + g];
        }

        for (std::size_t r = 0; r < ngenes; ++r) {
            for (const auto& br : block_remap) {
                buffers[br.first][r] = means[(r * nblocks + br.second) * ngroups + g];
            }
        }
        auto expected_m = scran_blocks::parallel_weighted_means(ngenes, ptrs, current_weights.data(), /* skip_nan = */ false);
        scran_tests::compare_almost_equal(expected_m, average_mean[g]);
    }
}

TEST_F(AverageGroupStatsTest, ZeroedBlockQuantile) {
    std::vector<std::vector<double> > average_quantile(ngroups);
    std::vector<double*> quantile_ptrs(ngroups); 
    for (size_t g = 0; g < ngroups; ++g) {
        average_quantile[g].resize(ngenes);
        quantile_ptrs[g] = average_quantile[g].data();
    }

    auto means_copy = means;
    std::vector<double> buffer;
    scran_blocks::SingleQuantileVariable<double, typename std::vector<double>::iterator> qcalc(nblocks, 0.5);
    for (std::size_t r = 0; r < ngenes; ++r) {
        // Setting all means in the second block to NaN.
        size_t offset = r * ngroups * nblocks;
        auto mcptr = means_copy.data() + offset;
        std::fill_n(mcptr + ngroups, ngroups, std::numeric_limits<double>::quiet_NaN());

        scran_markers::internal::average_group_stats_blockquantile(
            r, 
            ngroups,
            nblocks,
            mcptr,
            buffer,
            qcalc,
            quantile_ptrs
        );
    }

    // Doing a reference calculation, but skipping the middle block.
    std::vector<std::vector<double> > buffers(nblocks - 1);
    std::vector<double*> ptrs;
    for (auto& bk : buffers) {
        bk.resize(ngenes);
        ptrs.push_back(bk.data());
    }
    std::vector<std::pair<int, int> > block_remap{ { 0, 0 }, { 1, 2 }};

    for (std::size_t g = 0; g < ngroups; ++g) {
        for (std::size_t r = 0; r < ngenes; ++r) {
            for (const auto& br : block_remap) {
                buffers[br.first][r] = means[(r * nblocks + br.second) * ngroups + g];
            }
        }

        auto expected_q = scran_blocks::parallel_quantiles(ngenes, ptrs, 0.5, /* skip_nan = */ false);
        scran_tests::compare_almost_equal(expected_q, average_quantile[g]);
    }
}

TEST_F(AverageGroupStatsTest, ZeroedGroupMean) {
    // Zeroing all weights in the first group.
    auto copy_weights = weights;
    for (size_t b = 0; b < nblocks; ++b) {
        copy_weights[b * ngroups] = 0;
    }

    std::vector<std::vector<double> > average_mean(ngroups);
    std::vector<double*> mean_ptrs(ngroups);
    for (size_t g = 0; g < ngroups; ++g) {
        average_mean[g].resize(ngenes);
        mean_ptrs[g] = average_mean[g].data();
    }

    auto means_copy = means;
    auto total_weights = scran_markers::internal::compute_total_weight_per_group(ngroups, nblocks, copy_weights.data());
    for (size_t r = 0; r < ngenes; ++r) {
        // Setting all means in the first group to NaN as well, to make sure they're unused.
        std::size_t offset = r * ngroups * nblocks;
        auto mcptr = means_copy.data() + offset;
        for (std::size_t b = 0; b < nblocks; ++b) {
            mcptr[b * ngroups] = std::numeric_limits<double>::quiet_NaN();
        }

        scran_markers::internal::average_group_stats_blockmean(
            r, 
            ngroups,
            nblocks,
            mcptr,
            copy_weights.data(),
            total_weights.data(),
            mean_ptrs
        );
    }

    int num_na = 0;
    for (auto x : average_mean[0]) {
        num_na += std::isnan(x);
    }
    EXPECT_EQ(num_na, ngenes);

    // Other groups are not affected.
    num_na = 0;
    for (auto x : average_mean[1]) {
        num_na += std::isnan(x);
    }
    EXPECT_EQ(num_na, 0);
}

TEST_F(AverageGroupStatsTest, ZeroedGroupQuantile) {
    std::vector<std::vector<double> > average_quantile(ngroups);
    std::vector<double*> quantile_ptrs(ngroups);
    for (size_t g = 0; g < ngroups; ++g) {
        average_quantile[g].resize(ngenes);
        quantile_ptrs[g] = average_quantile[g].data();
    }

    auto means_copy = means;
    std::vector<double> buffer;
    scran_blocks::SingleQuantileVariable<double, typename std::vector<double>::iterator> qcalc(nblocks, 0.5);
    for (size_t r = 0; r < ngenes; ++r) {
        // Setting all quantiles in the first group to NaN. 
        std::size_t offset = r * ngroups * nblocks;
        auto mcptr = means_copy.data() + offset;
        for (std::size_t b = 0; b < nblocks; ++b) {
            mcptr[b * ngroups] = std::numeric_limits<double>::quiet_NaN();
        }

        scran_markers::internal::average_group_stats_blockquantile(
            r, 
            ngroups,
            nblocks,
            mcptr,
            buffer,
            qcalc,
            quantile_ptrs
        );
    }

    int num_na = 0;
    for (auto x : average_quantile[0]) {
        num_na += std::isnan(x);
    }
    EXPECT_EQ(num_na, ngenes);

    // Other groups are not affected.
    num_na = 0;
    for (auto x : average_quantile[1]) {
        num_na += std::isnan(x);
    }
    EXPECT_EQ(num_na, 0);
}
