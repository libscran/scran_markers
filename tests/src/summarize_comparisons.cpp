#include "scran_tests/scran_tests.hpp"

#include <numeric>
#include <random>

#include "scran_markers/summarize_comparisons.hpp"

class SummarizeComparisonsTest : public ::testing::TestWithParam<std::tuple<int, int> > {
protected:
    static void create_outputs(int ngenes, int ngroups, std::vector<double>& output, std::vector<scran_markers::SummaryBuffers<double, int> >& ptrs) {
        output.resize(ngroups * 4 * ngenes
#ifdef SCRAN_MARKERS_TEST_INIT
            , SCRAN_MARKERS_TEST_INIT
#endif
        );

        auto optr = output.data();
        ptrs.resize(ngroups);
        for (int g = 0; g < ngroups; ++g) {
            auto& current = ptrs[g];
            current.min = optr;
            optr += ngenes;
            current.mean = optr;
            optr += ngenes;
            current.median = optr;
            optr += ngenes;
            current.max = optr;
            optr += ngenes;
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
    std::vector<scran_markers::SummaryBuffers<double, int> > ptrs;
    create_outputs(ngenes, ngroups, output, ptrs);

    auto values = spawn_simple_values(ngenes, ngroups);
    auto threads = std::get<1>(params);
    scran_markers::internal::summarize_comparisons(ngenes, ngroups, values.data(), ptrs, threads);

    for (int gene = 0; gene < ngenes; ++gene) {
        for (int g = 0; g < ngroups; ++g) {
            // Checking that the minimum is correct.
            EXPECT_FLOAT_EQ(ptrs[g].min[gene], g * gene + (g == 0));

            // Checking that the mean is correct.
            EXPECT_FLOAT_EQ(ptrs[g].mean[gene], g * gene + ((ngroups - 1)*ngroups/2.0 - g)/(ngroups-1));

            // Checking that the median is between the min and max.
            EXPECT_TRUE(ptrs[g].median[gene] >= ptrs[g].min[gene]);
            EXPECT_TRUE(ptrs[g].median[gene] <= ptrs[g].max[gene]);

            // Checking that the maximum is correct.
            EXPECT_FLOAT_EQ(ptrs[g].max[gene], g * gene + ngroups - 1 - (g == ngroups - 1));
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
    std::vector<scran_markers::SummaryBuffers<double, int> > ptrs;
    create_outputs(ngenes, ngroups, output, ptrs);

    auto threads = std::get<1>(GetParam());

    for (int lost = 0; lost < ngroups; ++lost) {
        auto values = spawn_missing_values(ngenes, ngroups, lost);
        scran_markers::internal::summarize_comparisons(ngenes, ngroups, values.data(), ptrs, threads);

        for (int gene = 0; gene < ngenes; ++ gene) {
            for (int g = 0; g < ngroups; ++g) {
                if (g == lost) {
                    EXPECT_TRUE(std::isnan(ptrs[g].min[gene]));
                    EXPECT_TRUE(std::isnan(ptrs[g].mean[gene]));
                    EXPECT_TRUE(std::isnan(ptrs[g].median[gene]));
                    EXPECT_TRUE(std::isnan(ptrs[g].max[gene]));
                    continue;
                }

                // Checking that the minimum is correct.
                double baseline = g * gene;
                if ((g==0 && lost==1) || (g==1 && lost==0)) {
                    baseline += 2; // first non-self, non-NaN effect is 2.
                } else if (g==0 || lost==0) {
                    baseline += 1; // self is NaN, so first non-self, non-NaN effect is 1.
                }
                EXPECT_FLOAT_EQ(ptrs[g].min[gene], baseline);

                // Checking that the mean is correct.
                baseline = g * gene;
                if (lost == g) {
                    baseline += ((ngroups - 1)*ngroups/2.0 - g)/(ngroups-1);
                } else {
                    baseline += ((ngroups - 1)*ngroups/2.0 - g - lost)/(ngroups-2);
                }
                EXPECT_FLOAT_EQ(ptrs[g].mean[gene], baseline);

                // Checking that the maximum is correct.
                baseline = g * gene;
                if ((g==ngroups - 1 && lost==ngroups - 2) || (g==ngroups - 2  && lost==ngroups - 1)) {
                    baseline += ngroups - 3; // i.e., the last non-self, non-NaN effect. 
                } else if (g==ngroups - 1 || lost==ngroups - 1) {
                    baseline += ngroups - 2; // self is NaN, so we get a different maximum.
                } else {
                    baseline += ngroups - 1;
                }
                EXPECT_FLOAT_EQ(ptrs[g].max[gene], baseline);
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

class ComputeMinRankTest : public ::testing::TestWithParam<std::tuple<int, int> > {
protected:
    static void configure(int ngenes, int ngroups, std::vector<int>& output, std::vector<scran_markers::SummaryBuffers<double, int> >& ptrs) {
        output.resize(ngroups * ngenes);
        auto optr = output.data();
        ptrs.resize(ngroups);
        for (int g = 0; g < ngroups; ++g, optr += ngenes) {
            ptrs[g].min_rank = optr;
        }
    }

    void compare_vectors(const std::vector<int>& expected, const int* ptr) {
        EXPECT_EQ(expected, scran_tests::vector_n(ptr, expected.size()));
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

    std::vector<int> output;
    std::vector<scran_markers::SummaryBuffers<double, int> > ptrs;
    configure(ngenes, ngroups, output, ptrs);

    const auto params = GetParam();
    const auto preserve_ties = std::get<0>(params);
    const auto threads = std::get<1>(params);

    scran_markers::internal::compute_min_rank_pairwise(ngenes, ngroups, effects.data(), ptrs, preserve_ties, threads);
    for (size_t i = 0; i < ngroups; ++i) {
        compare_vectors({4, 3, 2, 1}, ptrs[i].min_rank); // reversed, for maximum effect sizes.
    }
}

TEST_P(ComputeMinRankTest, Tied) {
    size_t ngenes = 4, ngroups = 3;
    std::vector<double> effects { // rows are genes, columns are group * group of pairwise effects.
        0, 1, 1, 1, 0, 1, 1, 1, 0,
        0, 4, 4, 4, 0, 4, 4, 4, 0,
        0, 1, 1, 1, 0, 1, 1, 1, 0,
        0, 4, 4, 4, 0, 4, 4, 4, 0
    };

    std::vector<int> output;
    std::vector<scran_markers::SummaryBuffers<double, int> > ptrs;
    configure(ngenes, ngroups, output, ptrs);

    const auto params = GetParam();
    const auto preserve_ties = std::get<0>(params);
    const auto threads = std::get<1>(params);
    scran_markers::internal::compute_min_rank_pairwise(ngenes, ngroups, effects.data(), ptrs, preserve_ties, threads);

    if (preserve_ties) {
        for (size_t i = 0; i < ngroups; ++i) {
            compare_vectors({3, 1, 3, 1}, ptrs[i].min_rank);
        }
    } else {
        for (size_t i = 0; i < ngroups; ++i) {
            compare_vectors({3, 1, 4, 2}, ptrs[i].min_rank);
        }
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

    std::vector<int> output;
    std::vector<scran_markers::SummaryBuffers<double, int> > ptrs;
    configure(ngenes, ngroups, output, ptrs);

    const auto params = GetParam();
    const auto preserve_ties = std::get<0>(params);
    const auto threads = std::get<1>(params);

    scran_markers::internal::compute_min_rank_pairwise(ngenes, ngroups, effects.data(), ptrs, preserve_ties, threads);
    compare_vectors({1, 2, 1, 3}, ptrs[0].min_rank);
    compare_vectors({2, 3, 1, 1}, ptrs[1].min_rank);
    compare_vectors({1, 1, 2, 4}, ptrs[2].min_rank);
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
        0, n, 4, 2, 0, 1, 2, n, 0,
        0, 2, 3, n, 0, 3, 3, n, 0,
        0, 3, 1, 1, 0, n, 4, n, 0
     * after we remove the NA and promote all subsequent entries.
     */

    std::vector<int> output;
    std::vector<scran_markers::SummaryBuffers<double, int> > ptrs;
    configure(ngenes, ngroups, output, ptrs);

    const auto params = GetParam();
    const auto preserve_ties = std::get<0>(params);
    const auto threads = std::get<1>(params);

    scran_markers::internal::compute_min_rank_pairwise(ngenes, ngroups, effects.data(), ptrs, preserve_ties, threads);
    compare_vectors({1, 4, 2, 1}, ptrs[0].min_rank);
    compare_vectors({2, 1, 3, 1}, ptrs[1].min_rank);
    compare_vectors({1, 2, 3, 4}, ptrs[2].min_rank);
}

INSTANTIATE_TEST_SUITE_P(
    ComputeMinRank,
    ComputeMinRankTest,
    ::testing::Combine(
        ::testing::Values(false, true), // consider ties
        ::testing::Values(1, 3) // number of threads
    )
);

// Also checking that our multiple min_rank variants are consistent
// with each other, especially when threading gets involved.
class ComputeMinRankTestThreaded : public ::testing::TestWithParam<std::tuple<int, bool> > {
protected:
    static void configure(int ngenes, int ngroups, std::vector<int>& output, std::vector<scran_markers::SummaryBuffers<double, int> >& ptrs) {
        output.resize(ngroups * ngenes);
        auto optr = output.data();
        ptrs.resize(ngroups);
        for (int g = 0; g < ngroups; ++g, optr += ngenes) {
            ptrs[g].min_rank = optr;
        }
    }
};

TEST_P(ComputeMinRankTestThreaded, Consistency) {
    auto param = GetParam();
    size_t ngenes = 200;
    size_t ngroups = std::get<0>(param);
    bool add_nans = std::get<1>(param);

    scran_tests::SimulateVectorParameters sparams;
    sparams.seed = ngroups * 10 + add_nans * 100;
    auto effects = scran_tests::simulate_vector(ngenes * ngroups * ngroups, sparams);
    if (add_nans) {
        std::mt19937_64 rng(sparams.seed + 1);
        std::uniform_real_distribution<double> unif;
        for (auto& e : effects) {
            if (unif(rng) <= 0.05) {
                e = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }

    {
        std::vector<int> ref_output;
        std::vector<scran_markers::SummaryBuffers<double, int> > ref_ptrs;
        configure(ngenes, ngroups, ref_output, ref_ptrs);
        scran_markers::internal::compute_min_rank_pairwise(ngenes, ngroups, effects.data(), ref_ptrs, false, 1);

        std::vector<int> threaded_output;
        std::vector<scran_markers::SummaryBuffers<double, int> > threaded_ptrs;
        configure(ngenes, ngroups, threaded_output, threaded_ptrs);
        scran_markers::internal::compute_min_rank_pairwise(ngenes, ngroups, effects.data(), threaded_ptrs, false, 3);
        EXPECT_EQ(threaded_output, ref_output);
    }

    // Testing the tied logic.
    {
        std::vector<int> ref_output;
        std::vector<scran_markers::SummaryBuffers<double, int> > ref_ptrs;
        configure(ngenes, ngroups, ref_output, ref_ptrs);
        scran_markers::internal::compute_min_rank_pairwise(ngenes, ngroups, effects.data(), ref_ptrs, true, 1);

        std::vector<int> threaded_output;
        std::vector<scran_markers::SummaryBuffers<double, int> > threaded_ptrs;
        configure(ngenes, ngroups, threaded_output, threaded_ptrs);
        scran_markers::internal::compute_min_rank_pairwise(ngenes, ngroups, effects.data(), threaded_ptrs, true, 3);
        EXPECT_EQ(threaded_output, ref_output);
    }
}

INSTANTIATE_TEST_SUITE_P(
    ComputeMinRankThreaded,
    ComputeMinRankTestThreaded,
    ::testing::Combine(
        ::testing::Values(2, 3, 4, 5), // number of groups
        ::testing::Values(false, true) // whether to spike in NaN's.
    )
);
