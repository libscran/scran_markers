#include "scran_tests/scran_tests.hpp"

#include "scran_markers/create_combinations.hpp"

TEST(CreateCombinations, Basic) {
    std::vector<int> group{ 0, 1, 2, 3 };
    std::vector<int> block{ 3, 2, 1, 0 };
    auto out = scran_markers::create_combinations(4, group.data(), 10, block.data(), 5);

    std::vector<std::size_t> expected{ 30, 21, 12, 3 };
    EXPECT_EQ(out.combinations, expected);
    EXPECT_EQ(out.num_combinations, 50);
}
