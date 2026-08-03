#ifndef SCRAN_MARKERS_TEST_UTILS_H
#define SCRAN_MARKERS_TEST_UTILS_H

#include <vector>
#include <algorithm>
#include <cstddef>

inline std::vector<int> create_interleaved_factor(std::size_t n, int num_factors) {
    std::vector<int> factor(n);
    for (std::size_t x = 0; x < n; ++x) {
        factor[x] = x % num_factors;
    }
    return factor;
}

inline std::vector<int> create_contiguous_factor(std::size_t n, int num_factors) {
    const std::size_t per_factor = n / num_factors;
    const int remainder = n % num_factors;
    std::vector<int> factor;
    factor.reserve(n);
    for (int b = 0; b < num_factors; ++b) {
        factor.insert(factor.end(), per_factor + (b < remainder), b);
    }
    assert(factor.size() == n);
    return factor;
}

inline std::vector<int> create_groupings(size_t n, int ngroups) {
    return create_interleaved_factor(n, ngroups);
}

inline std::vector<int> create_blocks(size_t n, int nblocks) {
    return create_contiguous_factor(n, nblocks);
}

#endif
