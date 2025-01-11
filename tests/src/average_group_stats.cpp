#include "scran_tests/scran_tests.hpp"
#include "scran_blocks/scran_blocks.hpp"

#include "scran_markers/average_group_stats.hpp"

class AverageGroupStatsTest : public ::testing::Test {
protected:
    inline static size_t ngroups = 5, nblocks = 3, ngenes = 100; 
    inline static std::vector<double> means, detected, weights;

    static void SetUpTestSuite() {
        means = scran_tests::simulate_vector(ngenes * ngroups * nblocks, []{ 
            scran_tests::SimulationParameters sparam;
            sparam.seed = 43210;
            return sparam;
        }());

        detected = scran_tests::simulate_vector(ngenes * ngroups * nblocks, []{ 
            scran_tests::SimulationParameters sparam;
            sparam.seed = 98765;
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

TEST_F(AverageGroupStatsTest, Basic) {
    std::vector<std::vector<double> > average_mean(ngroups), average_detected(ngroups);
    std::vector<double*> mean_ptrs(ngroups), detected_ptrs(ngroups);
    for (size_t g = 0; g < ngroups; ++g) {
        average_mean[g].resize(ngenes);
        mean_ptrs[g] = average_mean[g].data();
        average_detected[g].resize(ngenes);
        detected_ptrs[g] = average_detected[g].data();
    }

    auto total_weights = scran_markers::internal::compute_total_weight_per_group(ngroups, nblocks, weights.data());
    for (size_t r = 0; r < ngenes; ++r) {
        size_t offset = r * ngroups * nblocks;
        scran_markers::internal::average_group_stats(
            r, 
            ngroups,
            nblocks,
            means.data() + offset,
            detected.data() + offset,
            weights.data(),
            total_weights.data(),
            mean_ptrs,
            detected_ptrs
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

    for (size_t g = 0; g < ngroups; ++g) {
        for (size_t b = 0; b < nblocks; ++b) {
            current_weights[b] = weights[b * ngroups + g];
        }

        for (size_t r = 0; r < ngenes; ++r) {
            for (size_t b = 0; b < nblocks; ++b) {
                buffers[b][r] = means[(r * nblocks + b) * ngroups + g];
            }
        }
        auto expected_m = scran_blocks::average_vectors_weighted(ngenes, ptrs, current_weights.data(), /* skip_nan = */ false);
        scran_tests::compare_almost_equal(expected_m, average_mean[g]);

        for (size_t r = 0; r < ngenes; ++r) {
            for (size_t b = 0; b < nblocks; ++b) {
                buffers[b][r] = detected[(r * nblocks + b) * ngroups + g];
            }
        }
        auto expected_d = scran_blocks::average_vectors_weighted(ngenes, ptrs, current_weights.data(), /* skip_nan = */ false);
        scran_tests::compare_almost_equal(expected_d, average_detected[g]);
    }
}

TEST_F(AverageGroupStatsTest, Zeroed) {
    // Zeroing all weights in the second block.
    auto copy_weights = weights;
    std::fill_n(copy_weights.begin() + ngroups, ngroups, 0);

    std::vector<std::vector<double> > average_mean(ngroups), average_detected(ngroups);
    std::vector<double*> mean_ptrs(ngroups), detected_ptrs(ngroups);
    for (size_t g = 0; g < ngroups; ++g) {
        average_mean[g].resize(ngenes);
        mean_ptrs[g] = average_mean[g].data();
        average_detected[g].resize(ngenes);
        detected_ptrs[g] = average_detected[g].data();
    }

    auto total_weights = scran_markers::internal::compute_total_weight_per_group(ngroups, nblocks, copy_weights.data());
    for (size_t r = 0; r < ngenes; ++r) {
        size_t offset = r * ngroups * nblocks;
        scran_markers::internal::average_group_stats(
            r, 
            ngroups,
            nblocks,
            means.data() + offset,
            detected.data() + offset,
            copy_weights.data(),
            total_weights.data(),
            mean_ptrs,
            detected_ptrs
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

    for (size_t g = 0; g < ngroups; ++g) {
        for (const auto& br : block_remap) {
            current_weights[br.first] = copy_weights[br.second * ngroups + g];
        }

        for (size_t r = 0; r < ngenes; ++r) {
            for (const auto& br : block_remap) {
                buffers[br.first][r] = means[(r * nblocks + br.second) * ngroups + g];
            }
        }
        auto expected_m = scran_blocks::average_vectors_weighted(ngenes, ptrs, current_weights.data(), /* skip_nan = */ false);
        scran_tests::compare_almost_equal(expected_m, average_mean[g]);

        for (size_t r = 0; r < ngenes; ++r) {
            for (const auto& br : block_remap) {
                buffers[br.first][r] = detected[(r * nblocks + br.second) * ngroups + g];
            }
        }
        auto expected_d = scran_blocks::average_vectors_weighted(ngenes, ptrs, current_weights.data(), /* skip_nan = */ false);
        scran_tests::compare_almost_equal(expected_d, average_detected[g]);
    }
}

TEST_F(AverageGroupStatsTest, AllZeroed) {
    // Zeroing all weights in the first group.
    auto copy_weights = weights;
    for (size_t b = 0; b < nblocks; ++b) {
        copy_weights[b * ngroups] = 0;
    }

    std::vector<std::vector<double> > average_mean(ngroups), average_detected(ngroups);
    std::vector<double*> mean_ptrs(ngroups), detected_ptrs(ngroups);
    for (size_t g = 0; g < ngroups; ++g) {
        average_mean[g].resize(ngenes);
        mean_ptrs[g] = average_mean[g].data();
        average_detected[g].resize(ngenes);
        detected_ptrs[g] = average_detected[g].data();
    }

    auto total_weights = scran_markers::internal::compute_total_weight_per_group(ngroups, nblocks, copy_weights.data());
    for (size_t r = 0; r < ngenes; ++r) {
        size_t offset = r * ngroups * nblocks;
        scran_markers::internal::average_group_stats(
            r, 
            ngroups,
            nblocks,
            means.data() + offset,
            detected.data() + offset,
            copy_weights.data(),
            total_weights.data(),
            mean_ptrs,
            detected_ptrs
        );
    }

    int num_na = 0;
    for (auto x : average_mean[0]) {
        num_na += std::isnan(x);
    }
    EXPECT_EQ(num_na, ngenes);

    num_na = 0;
    for (auto x : average_detected[0]) {
        num_na += std::isnan(x);
    }
    EXPECT_EQ(num_na, ngenes);

    // Other groups are not affected.
    num_na = 0;
    for (auto x : average_detected[1]) {
        num_na += std::isnan(x);
    }
    EXPECT_EQ(num_na, 0);

}
