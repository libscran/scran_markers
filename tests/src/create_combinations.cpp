#include "scran_tests/scran_tests.hpp"

#include "scran_markers/create_combinations.hpp"

TEST(CreateCombinations, Basic) {
    std::vector<int> group{ 0, 1, 2, 3 };
    std::vector<int> block{ 3, 2, 1, 0 };
    auto out = scran_markers::internal::create_combinations(10, group.data(), 4, block.data(), group.size());

    std::vector<size_t> expected{ 30, 21, 12, 3 };
    EXPECT_EQ(out, expected);

    auto counts = scran_markers::internal::tabulate_combinations<int>(10, 5, out);
    EXPECT_EQ(counts.size(), 50);
}
