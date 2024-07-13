#include <gtest/gtest.h>

#include <vector>
#include <algorithm>

#include "scran_markers/auc.hpp"

class AucTest : public ::testing::Test {
protected:
    static double slow_reference(
        const std::vector<double>& left,
        const std::vector<double>& right,
        double threshold = 0) 
    {
        double auc = 0;
        for (auto l : left) {
            size_t rx = 0;
            double collected = 0;
            while (rx < right.size() && right[rx] + threshold < l) {
                ++collected;
                ++rx;
            }

            double ties = 0;
            while (rx < right.size() && right[rx] + threshold == l) {
                ++ties;
                ++rx;
            }

            auc += collected + ties * 0.5;
        }

        return auc / (left.size() * right.size());
    }

    static void add_to_store(
        std::vector<double>& contents, 
        std::vector<std::pair<double, int> >& paired, 
        std::vector<int>& num_zeros, 
        std::vector<int>& totals)
    {
        std::sort(contents.begin(), contents.end());
        size_t group = num_zeros.size();
        num_zeros.resize(group + 1);
        totals.push_back(contents.size());
        for (auto c : contents) {
            if (c) {
                paired.push_back(std::make_pair(c, group));
            } else {
                ++num_zeros[group];
            }
        }
    }
};

TEST_F(AucTest, Self) {
    std::vector<double> group1 { 0, -0.1, 1, 2.2, 3.5, 5 }; 

    {
        std::vector<double> output(4);
        scran_markers::internal::AucWorkspace<double, int, double> input(2, output.data());

        std::vector<int> num_zeros, totals;
        add_to_store(group1, input.paired, num_zeros, totals); 
        add_to_store(group1, input.paired, num_zeros, totals);
        EXPECT_FLOAT_EQ(slow_reference(group1, group1), 0.5); // checking that the default calculation is correct.

        scran_markers::internal::compute_pairwise_auc(input, num_zeros, totals, true);

        EXPECT_FLOAT_EQ(output[1], 0.5);
        EXPECT_FLOAT_EQ(output[2], 0.5);
    }

    // Trying again with 3 groups.
    {
        std::vector<double> output(9);
        scran_markers::internal::AucWorkspace<double, int, double> input(3, output.data());

        std::vector<int> num_zeros, totals;
        add_to_store(group1, input.paired, num_zeros, totals); 
        add_to_store(group1, input.paired, num_zeros, totals);
        add_to_store(group1, input.paired, num_zeros, totals);

        scran_markers::internal::compute_pairwise_auc(input, num_zeros, totals, true);

        EXPECT_FLOAT_EQ(output[0 + 1], 0.5); 
        EXPECT_FLOAT_EQ(output[0 + 2], 0.5); 
        EXPECT_FLOAT_EQ(output[3 + 0], 0.5); 
        EXPECT_FLOAT_EQ(output[3 + 2], 0.5); 
        EXPECT_FLOAT_EQ(output[6 + 0], 0.5); 
        EXPECT_FLOAT_EQ(output[6 + 1], 0.5); 
    }
}

TEST_F(AucTest, NoZeros) {
    std::vector<double> group1 { -0.1, 1, 1, 2.3, 4 };
    std::vector<double> group2 { 1, 5, 3.4, 5, -0.1, 5, -0.2, -5 };
    std::vector<double> group3 { -0.12, 4, -0.1, 5, 2, -0.1, 3, 5, 6.2, 1.2, 1.11 };

    std::vector<double> output(9);
    scran_markers::internal::AucWorkspace<double, int, double> input(3, output.data());

    std::vector<int> num_zeros, totals;
    add_to_store(group1, input.paired, num_zeros, totals); 
    add_to_store(group2, input.paired, num_zeros, totals); 
    add_to_store(group3, input.paired, num_zeros, totals); 

    scran_markers::internal::compute_pairwise_auc(input, num_zeros, totals, true);

    EXPECT_FLOAT_EQ(output[3 + 0], slow_reference(group2, group1)); 
    EXPECT_FLOAT_EQ(output[6 + 0], slow_reference(group3, group1)); 
    EXPECT_FLOAT_EQ(output[6 + 1], slow_reference(group3, group2)); 

    EXPECT_FLOAT_EQ(output[0 + 1], 1 - output[3 + 0]);
    EXPECT_FLOAT_EQ(output[0 + 2], 1 - output[6 + 0]);
    EXPECT_FLOAT_EQ(output[3 + 2], 1 - output[6 + 1]); 
}

TEST_F(AucTest, Zeros) {
    std::vector<double> group1 { 0, -0.1, 0, 0, 2.3, 4, 0 };
    std::vector<double> group2 { 0, 1, 5, 0, 5, 0, 5, -0.2, -5 };
    std::vector<double> group3 { -0.12, 4, 0, 5, 2, 0, 3, 5, 0, 1.2, 1.11 };

    std::vector<double> output(9);
    scran_markers::internal::AucWorkspace<double, int, double> input(3, output.data());

    std::vector<int> num_zeros, totals;
    add_to_store(group1, input.paired, num_zeros, totals); 
    add_to_store(group2, input.paired, num_zeros, totals); 
    add_to_store(group3, input.paired, num_zeros, totals); 

    scran_markers::internal::compute_pairwise_auc(input, num_zeros, totals, true);

    EXPECT_FLOAT_EQ(output[3 + 0], slow_reference(group2, group1)); 
    EXPECT_FLOAT_EQ(output[6 + 0], slow_reference(group3, group1)); 
    EXPECT_FLOAT_EQ(output[6 + 1], slow_reference(group3, group2)); 
}

TEST_F(AucTest, ThresholdSelf) {
    std::vector<double> group { -1, 0, 1, 4, 3, 2, 5, 6, 7, 9 };

    std::vector<double> output(4);
    scran_markers::internal::AucWorkspace<double, int, double> input(2, output.data());

    std::vector<int> num_zeros, totals;
    add_to_store(group, input.paired, num_zeros, totals); 
    add_to_store(group, input.paired, num_zeros, totals); 

    for (double threshold = 0.5; threshold <= 2; ++threshold) { 
        scran_markers::internal::compute_pairwise_auc(input, num_zeros, totals, threshold, true);
        EXPECT_FLOAT_EQ(output[2], slow_reference(group, group, threshold));
        EXPECT_FLOAT_EQ(output[1], output[2]);
    }

    // Consistent results with a threshold of zero.
    scran_markers::internal::compute_pairwise_auc(input, num_zeros, totals, 0, true);

    EXPECT_FLOAT_EQ(output[2], 0.5);
    EXPECT_FLOAT_EQ(output[1], output[2]);
}

TEST_F(AucTest, ThresholdNoZero) {
    // Use 0.5 increments so that we get some juicy ties after adding the threshold.
    std::vector<double> group1 { 0.5, -0.5, 3, 2, -1.5 };
    std::vector<double> group2 { -0.5, 1.5, 1.5, 1.5, 2.5, -0.5, -0.5 };
    std::vector<double> group3 { -0.5, 6, 2, -1.5, 0.5, 0.15, 1, 2, 5 };

    std::vector<double> output(9);
    scran_markers::internal::AucWorkspace<double, int, double> input(3, output.data());

    std::vector<int> num_zeros, totals;
    add_to_store(group1, input.paired, num_zeros, totals); 
    add_to_store(group2, input.paired, num_zeros, totals); 
    add_to_store(group3, input.paired, num_zeros, totals); 

    for (double threshold = 0; threshold <= 2; threshold += 0.5) {
        scran_markers::internal::compute_pairwise_auc(input, num_zeros, totals, threshold, true);

        EXPECT_FLOAT_EQ(output[0 + 1], slow_reference(group1, group2, threshold)); 
        EXPECT_FLOAT_EQ(output[0 + 2], slow_reference(group1, group3, threshold)); 
        EXPECT_FLOAT_EQ(output[3 + 2], slow_reference(group2, group3, threshold)); 
        EXPECT_FLOAT_EQ(output[3 + 0], slow_reference(group2, group1, threshold)); 
        EXPECT_FLOAT_EQ(output[6 + 0], slow_reference(group3, group1, threshold)); 
        EXPECT_FLOAT_EQ(output[6 + 1], slow_reference(group3, group2, threshold)); 
    }
}

TEST_F(AucTest, ThresholdZeros) {
    std::vector<double> group1 { 0, 0.5, -0.5, 3, 2, -1.5, 0 };
    std::vector<double> group2 { -0.5, 0, 1.5, 1.5, 0, 1.5, 2.5, 0, -0.5, -0.5 };
    std::vector<double> group3 { -0.5, 6, 2, 0, -1.5, 0.5, 0.15, 1, 2, 0, 5 };

    std::vector<double> output(9);
    scran_markers::internal::AucWorkspace<double, int, double> input(3, output.data());

    std::vector<int> num_zeros, totals;
    add_to_store(group1, input.paired, num_zeros, totals); 
    add_to_store(group2, input.paired, num_zeros, totals); 
    add_to_store(group3, input.paired, num_zeros, totals); 

    for (double threshold = 0; threshold <= 0; threshold += 0.5) {
        scran_markers::internal::compute_pairwise_auc(input, num_zeros, totals, threshold, true);

        EXPECT_FLOAT_EQ(output[0 + 1], slow_reference(group1, group2, threshold)); 
        EXPECT_FLOAT_EQ(output[0 + 2], slow_reference(group1, group3, threshold)); 
        EXPECT_FLOAT_EQ(output[3 + 2], slow_reference(group2, group3, threshold)); 
        EXPECT_FLOAT_EQ(output[3 + 0], slow_reference(group2, group1, threshold)); 
        EXPECT_FLOAT_EQ(output[6 + 0], slow_reference(group3, group1, threshold)); 
        EXPECT_FLOAT_EQ(output[6 + 1], slow_reference(group3, group2, threshold)); 
    }
}
