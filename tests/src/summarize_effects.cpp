#include "scran_tests/scran_tests.hpp"

#include <numeric>
#include <random>

#include "scran_markers/summarize_effects.hpp"

TEST(SummarizeEffects, Basic) {
    size_t ngenes = 111, ngroups = 3;
    std::vector<double> stuff = scran_tests::simulate_vector(ngroups * ngroups * ngenes, scran_tests::SimulateVectorParameters());

    scran_markers::SummarizeEffectsOptions opts;
    auto res = scran_markers::summarize_effects(ngenes, ngroups, stuff.data(), opts);
    EXPECT_EQ(res.size(), ngroups);

    // Summaries are actually different for each group.
    EXPECT_NE(res[0].mean, res[1].mean);
    EXPECT_NE(res[1].mean, res[2].mean);

    // Summaries are actually filled.
    for (size_t g = 0; g < ngroups; ++g) {
        const auto& gres = res[g];

        auto check = [&](const auto& x) {
            EXPECT_EQ(x.size(), ngenes);
            bool nonzero = false;
            for (auto v : x) {
                if (v != 0) {
                    nonzero = true;
                }
            }
            EXPECT_TRUE(nonzero);
        };
        check(gres.min);
        check(gres.mean);
        check(gres.median);
        check(gres.max);

        // Min, max and mean make sense.
        for (size_t i = 0; i < ngenes; ++i) {
            EXPECT_LE(gres.min[i], gres.mean[i]);
            EXPECT_GE(gres.max[i], gres.mean[i]);
            EXPECT_LE(gres.min[i], gres.median[i]);
            EXPECT_GE(gres.max[i], gres.median[i]);
        }

        // Minimum rank makes sense.
        const auto& ranks = gres.min_rank;
        EXPECT_EQ(ranks.size(), ngenes);
        for (auto x : ranks) {
            EXPECT_TRUE(x >= 1);
            EXPECT_TRUE(x <= static_cast<int>(ngenes));
        }
    }

    // Same results with multiple threads.
    opts.num_threads = 3;
    auto res2 = scran_markers::summarize_effects(ngenes, ngroups, stuff.data(), opts);
    EXPECT_EQ(res2.size(), ngroups);

    for (size_t g = 0; g < ngroups; ++g) {
        EXPECT_EQ(res[g].min, res2[g].min);
        EXPECT_EQ(res[g].mean, res2[g].mean);
        EXPECT_EQ(res[g].median, res2[g].median);
        EXPECT_EQ(res[g].max, res2[g].max);
        EXPECT_EQ(res[g].min_rank, res2[g].min_rank);
    }
}

TEST(SummarizeEffects, Quantile) {
    size_t ngenes = 111, ngroups = 3;
    std::vector<double> stuff = scran_tests::simulate_vector(ngroups * ngroups * ngenes, scran_tests::SimulateVectorParameters());

    scran_markers::SummarizeEffectsOptions opts;
    opts.compute_quantiles.emplace(std::vector<double>{ 0.0, 0.5, 1.0 });
    auto res = scran_markers::summarize_effects(ngenes, ngroups, stuff.data(), opts);
    EXPECT_EQ(res.size(), ngroups);

    // Check that min, median and max make sense.
    for (size_t g = 0; g < ngroups; ++g) {
        const auto& gres = res[g];
        const auto& gq = *(gres.quantiles);
        for (size_t i = 0; i < ngenes; ++i) {
            EXPECT_EQ(gres.min[i], gq[0][i]);
            scran_tests::compare_almost_equal(gres.median[i], gq[1][i], scran_tests::CompareAlmostEqualParameters{});
            EXPECT_EQ(gres.max[i], gq[2][i]);
        }
    }

    // Same results with multiple threads.
    opts.num_threads = 3;
    auto res2 = scran_markers::summarize_effects(ngenes, ngroups, stuff.data(), opts);
    EXPECT_EQ(res2.size(), ngroups);

    for (size_t g = 0; g < ngroups; ++g) {
        const auto& gq = *(res[g].quantiles);
        const auto& gq2 = *(res2[g].quantiles);
        for (size_t i = 0; i < ngenes; ++i) {
            EXPECT_EQ(gq[0][i], gq2[0][i]);
            EXPECT_EQ(gq[1][i], gq2[1][i]);
            EXPECT_EQ(gq[2][i], gq2[2][i]);
        }
    }
}

TEST(SummarizeEffects, None) {
    size_t ngenes = 10, ngroups = 3;
    std::vector<double> stuff(ngroups * ngroups * ngenes);

    scran_markers::SummarizeEffectsOptions opts;
    opts.compute_min = false;
    opts.compute_mean = false;
    opts.compute_median = false;
    opts.compute_max = false;
    opts.compute_min_rank = false;

    auto res = scran_markers::summarize_effects(ngenes, ngroups, stuff.data(), opts);

    EXPECT_EQ(res.size(), ngroups);
    for (const auto& r : res) {
        EXPECT_EQ(r.min.size(), 0);
        EXPECT_EQ(r.mean.size(), 0);
        EXPECT_EQ(r.median.size(), 0);
        EXPECT_EQ(r.max.size(), 0);
        EXPECT_EQ(r.min_rank.size(), 0);
    }
}
