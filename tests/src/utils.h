#ifndef SCRAN_MARKERS_TEST_UTILS_H
#define SCRAN_MARKERS_TEST_UTILS_H

#include <vector>
#include <algorithm>

inline std::vector<int> create_groupings(size_t n, int ngroups) {
    std::vector<int> groupings(n);
    for (size_t g = 0; g < groupings.size(); ++g) {
        groupings[g] = g % ngroups;
    }
    return groupings;
}

inline std::vector<int> create_blocks(size_t n, int nblocks) {
    size_t per_block = (n / nblocks) + (n % nblocks > 0);
    std::vector<int> blocks;
    blocks.reserve(n);
    for (int b = 0; b < nblocks; ++b) {
        size_t extend_to = std::min(per_block + blocks.size(), n);
        blocks.resize(extend_to, b);
    }
    return blocks;
}

#endif
