#include <gtest/gtest.h>

#include "scran_markers/summarize_comparisons.hpp"

#include <numeric>
#include <random>

class SummarizeComparisonsTest : public ::testing::TestWithParam<std::tuple<int, int> > {
protected:
    static void create_outputs(int ngenes, int ngroups, std::vector<double>& output, std::vector<std::vector<double*> >& ptrs) {
        output.resize(ngroups * static_cast<size_t>(scran_markers::Summary::n_summaries) * ngenes);

        auto ptr = output.data();
        ptrs.resize(static_cast<size_t>(scran_markers::Summary::n_summaries));
        for (size_t s = 0; s < ptrs.size(); ++s) {
            ptrs[s].clear();
            for (int g = 0; g < ngroups; ++g, ptr += ngenes) {
                ptrs[s].push_back(ptr);
            }
        }
    }

    /* For each group, the effect sizes start from `group_id * gene` for the
     * comparison to group 0 and increase consecutively until the comparison to
     * group `ngroups - 1`. This should make it relatively straightforward to
     * compute the various metrics exactly for testing purposes. We add a `gene`
     * multiplier to ensure that `summarize_comparisons` responds to it properly.
     */
    static std::vector<double> spawn_simple_values(int ngenes, int ngroups) {
        std::vector<double> output(ngroups * ngroups * ngenes);
        auto start = output.begin();
        for (int gene = 0; gene < ngenes; ++gene) {
            for (int g = 0; g < ngroups; ++g, start += ngroups) {
                std::iota(start, start + ngroups, g * gene);
            }
        }
        return output;
    }

    static std::vector<double> spawn_missing_values(int ngenes, int ngroups, int lost) {
        std::vector<double> output(ngroups * ngroups * ngenes);
        auto start = output.begin();
        for (int gene = 0; gene < ngenes; ++gene) {
            for (int g = 0; g < ngroups; ++g, start += ngroups) {
                if (g == lost) {
                    std::fill(start, start + ngroups, std::numeric_limits<double>::quiet_NaN());
                } else {
                    std::iota(start, start + ngroups, g * gene);
                    *(start + lost) = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
        return output;
    }
};

TEST_P(SummarizeComparisonsTest, Basic) {
    auto params = GetParam();
    int ngenes = 100;
    int ngroups = std::get<0>(params);

    std::vector<double> output;
    std::vector<std::vector<double*> > ptrs;
    create_outputs(ngenes, ngroups, output, ptrs);

    auto values = spawn_simple_values(ngenes, ngroups);
    auto threads = std::get<1>(params);
    scran_markers::internal::summarize_comparisons(ngenes, ngroups, values.data(), ptrs, threads);

    for (int gene = 0; gene < ngenes; ++gene) {
        for (int g = 0; g < ngroups; ++g) {
            // Checking that the minimum is correct.
            EXPECT_FLOAT_EQ(ptrs[0][g][gene], g * gene + (g == 0));

            // Checking that the mean is correct.
            EXPECT_FLOAT_EQ(ptrs[1][g][gene], g * gene + ((ngroups - 1)*ngroups/2.0 - g)/(ngroups-1));

            // Checking that the median is between the min and max.
            EXPECT_TRUE(ptrs[2][g][gene] >= ptrs[0][g][gene]);
            EXPECT_TRUE(ptrs[2][g][gene] <= ptrs[3][g][gene]);

            // Checking that the maximum is correct.
            EXPECT_FLOAT_EQ(ptrs[3][g][gene], g * gene + ngroups - 1 - (g == ngroups - 1));
        }
    }

    // Checking the serial version for consistency.
    if (threads > 1) {
        auto parallelized = output;
        std::fill(output.begin(), output.end(), 0);
        scran_markers::internal::summarize_comparisons(ngenes, ngroups, values.data(), ptrs, 1); 
        EXPECT_EQ(parallelized, output);
    }
}

TEST_P(SummarizeComparisonsTest, Missing) {
    auto params = GetParam();
    int ngenes = 100;
    int ngroups = std::get<0>(params);

    std::vector<double> output;
    std::vector<std::vector<double*> > ptrs;
    create_outputs(ngenes, ngroups, output, ptrs);

    auto threads = std::get<1>(GetParam());

    for (int lost = 0; lost < ngroups; ++lost) {
        auto values = spawn_missing_values(ngenes, ngroups, lost);
        scran_markers::internal::summarize_comparisons(ngenes, ngroups, values.data(), ptrs, threads);

        for (int gene = 0; gene < ngenes; ++ gene) {
            for (int g = 0; g < ngroups; ++g) {
                if (g == lost) {
                    for (int i = 0; i < 3; ++i) {
                        EXPECT_TRUE(std::isnan(ptrs[i][g][gene]));
                    }
                    continue;
                }

                // Checking that the minimum is correct.
                double baseline = g * gene;
                if ((g==0 && lost==1) || (g==1 && lost==0)) {
                    baseline += 2;
                } else if (g==0 || lost==0) {
                    baseline += 1;
                }
                EXPECT_FLOAT_EQ(ptrs[0][g][gene], baseline);

                // Checking that the mean is correct.
                baseline = g * gene;
                if (lost == g) {
                    baseline += ((ngroups - 1)*ngroups/2.0 - g)/(ngroups-1);
                } else {
                    baseline += ((ngroups - 1)*ngroups/2.0 - g - lost)/(ngroups-2);
                }
                EXPECT_FLOAT_EQ(ptrs[1][g][gene], baseline);

                // Checking that the maximum is correct.
                baseline = g * gene + ngroups - 1;
                if ((g==ngroups - 1 && lost==ngroups - 2) || (g==ngroups - 2  && lost==ngroups - 1)) {
                    baseline -= 2;
                } else if (g==ngroups - 1 || lost==ngroups - 1) {
                    baseline -= 1;
                }
                EXPECT_FLOAT_EQ(ptrs[3][g][gene], baseline);
            }
        }
    }
}

INSTANTIATE_TEST_SUITE_P(
    SummarizeComparisons,
    SummarizeComparisonsTest,
    ::testing::Combine(
        ::testing::Values(3, 5, 7), // number of groups
        ::testing::Values(1, 3) // number of threads
    )
);

/*********************************************/

class ComputeMinRankTest : public ::testing::TestWithParam<int> {
protected:
    static void configure(int ngenes, int ngroups, std::vector<double>& output, std::vector<double*>& ptrs) {
        output.resize(ngroups * ngenes);
        auto ptr = output.data();
        for (int g = 0; g < ngroups; ++g, ptr += ngenes) {
            ptrs.push_back(ptr);
        }
    }

    void compare_vectors(const std::vector<double>& expected, const double* ptr) {
        std::vector<double> observed(ptr, ptr + expected.size());
        EXPECT_EQ(expected, observed);
    }
};

TEST_P(ComputeMinRankTest, Basic) {
    size_t ngenes = 4, ngroups = 3;
    std::vector<double> effects { // rows are genes, columns are group * group of pairwise effects.
        0, 1, 1, 1, 0, 1, 1, 1, 0,
        0, 2, 2, 2, 0, 2, 2, 2, 0,
        0, 3, 3, 3, 0, 3, 3, 3, 0,
        0, 4, 4, 4, 0, 4, 4, 4, 0
    };

    std::vector<double> output;
    std::vector<double*> ptrs;
    configure(ngenes, ngroups, output, ptrs);

    auto threads = GetParam();
    scran_markers::internal::compute_min_rank_pairwise(ngenes, ngroups, effects.data(), ptrs, threads);
    for (size_t i = 0; i < ngroups; ++i) {
        compare_vectors({4, 3, 2, 1}, ptrs[i]); // reversed, for maximum effect sizes.
    }
}

TEST_P(ComputeMinRankTest, LessBasic) {
    size_t ngenes = 4, ngroups = 3;
    std::vector<double> effects { // rows are genes, columns are group * group of pairwise effects.
        0, 1, 2, 2, 0, 4, 1, 3, 0,
        0, 2, 4, 3, 0, 3, 2, 1, 0,
        0, 3, 1, 1, 0, 2, 3, 2, 0,
        0, 4, 3, 4, 0, 1, 4, 4, 0
     /* 1  1  1  2  2  2  3  3  3 => comparisons for each group */
    };
    for (auto& e : effects) { e *= -1; }

    std::vector<double> output;
    std::vector<double*> ptrs;
    configure(ngenes, ngroups, output, ptrs);

    auto threads = GetParam();
    scran_markers::internal::compute_min_rank_pairwise(ngenes, ngroups, effects.data(), ptrs, threads);
    compare_vectors({1, 2, 1, 3}, output.data());
    compare_vectors({2, 3, 1, 1}, output.data() + ngenes);
    compare_vectors({1, 1, 2, 4}, output.data() + ngenes * 2);
}

TEST_P(ComputeMinRankTest, Missing) {
    size_t ngenes = 4, ngroups = 3;
    auto n = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> effects { // rows are genes, columns are group * group of pairwise effects.
        0, 1, 2, 3, 0, 3, 1, n, 0,
        0, n, 4, 2, 0, 2, 2, n, 0,
        0, 3, 3, n, 0, 4, 3, n, 0,
        0, 4, 1, 1, 0, n, 4, n, 0
     /* 1  1  1  2  2  2  3  3  3 => comparisons for each group */
    };
    for (auto& e : effects) { e *= -1; }

    /* Implicitly, the ranks become:
        0, 1, 2, 3, 0, 2, 1, n, 0,
        0, X, 4, 2, 0, 1, 2, n, 0,
        0, 2, 3, n, 0, 3, 3, n, 0,
        0, 3, 1, 1, 0, X, 4, n, 0
     * after we remove the NA and promote all subsequent entries.
     */

    std::vector<double> output;
    std::vector<double*> ptrs;
    configure(ngenes, ngroups, output, ptrs);

    auto threads = GetParam();
    scran_markers::internal::compute_min_rank_pairwise(ngenes, ngroups, effects.data(), ptrs, threads);
    compare_vectors({1, 4, 2, 1}, output.data());
    compare_vectors({2, 1, 3, 1}, output.data() + ngenes);
    compare_vectors({1, 2, 3, 4}, output.data() + ngenes * 2);
}

INSTANTIATE_TEST_SUITE_P(
    ComputeMinRank,
    ComputeMinRankTest,
    ::testing::Values(1, 3) // number of threads
);

// Also checking that our multiple min_rank variants are consistent
// with each other, especially when threading gets involved.
class ComputeMinRankTestThreaded : public ::testing::TestWithParam<std::tuple<int, int, bool> > {
protected:
    std::vector<double*> configure(int ngenes, int ngroups, double* ptr) {
        std::vector<double*> ptrs;
        for (int g = 0; g < ngroups; ++g, ptr += ngenes) {
            ptrs.push_back(ptr);
        }
        return ptrs;
    }
};

TEST_P(ComputeMinRankTestThreaded, Consistency) {
    auto param = GetParam();
    size_t ngenes = 20;
    size_t ngroups = std::get<0>(param);
    size_t nthreads = std::get<1>(param);
    bool add_nans = std::get<2>(param);

    std::mt19937_64 rng(/* seed */ ngroups * nthreads + add_nans);
    std::vector<double> effects(ngenes * ngroups * ngroups);
    std::normal_distribution<double> dist;
    std::uniform_real_distribution<double> unif;
    for (auto& e : effects) {
        e = dist(rng);
        if (add_nans && unif(rng) <= 0.05) {
            e = std::numeric_limits<double>::quiet_NaN();
        }
    }

    std::vector<double> ref_output(ngroups * ngenes);
    std::vector<double*> ref_ptrs = configure(ngenes, ngroups, ref_output.data());
    scran_markers::internal::compute_min_rank_pairwise(ngenes, ngroups, effects.data(), ref_ptrs, 1);

    std::vector<double> threaded_output(ngroups * ngenes);
    std::vector<double*> threaded_ptrs = configure(ngenes, ngroups, threaded_output.data());
    scran_markers::internal::compute_min_rank_pairwise(ngenes, ngroups, effects.data(), threaded_ptrs, nthreads);
    EXPECT_EQ(threaded_output, ref_output);

    std::vector<double> pergroup_output(ngroups * ngenes);
    for (size_t g = 0; g < ngroups; ++g) {
        std::vector<double> copy(ngroups * ngenes);
        for (size_t i = 0; i < ngenes; ++i) {
            auto base = effects.data() + i * ngroups * ngroups + g * ngroups;
            std::copy(base, base + ngroups, copy.data() + i * ngroups);
        }
        scran_markers::internal::compute_min_rank_for_group(ngenes, ngroups, g, copy.data(), pergroup_output.data() + g * ngenes, nthreads);
    }
    EXPECT_EQ(pergroup_output, ref_output);
}

INSTANTIATE_TEST_SUITE_P(
    ComputeMinRankThreaded,
    ComputeMinRankTestThreaded,
    ::testing::Combine(
        ::testing::Values(1, 2, 3), // number of threads
        ::testing::Values(2, 3, 4, 5), // number of groups
        ::testing::Values(false, true) // whether to spike in NaN's.
    )
);
