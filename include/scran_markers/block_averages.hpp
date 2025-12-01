#ifndef SCRAN_MARKERS_BLOCK_AVERAGES_HPP
#define SCRAN_MARKERS_BLOCK_AVERAGES_HPP

#include <vector>
#include <cstddef>
#include <optional>

#include "sanisizer/sanisizer.hpp"

#include "utils.hpp"

namespace scran_markers {

enum class AveragePolicy : unsigned char { MEAN, QUANTILE };

namespace internal {

template<typename Stat_>
class PrecomputedPairwiseWeights {
public:
    // 'combo_weights' are expected to be 'ngroups * nblocks' arrays where
    // groups are the faster-changing dimension and the blocks are slower.
    PrecomputedPairwiseWeights(const std::size_t ngroups, const std::size_t nblocks, const Stat_* const combo_weights) :
        my_total(sanisizer::product<decltype(I(my_total.size()))>(ngroups, ngroups)),
        my_by_block(sanisizer::product<decltype(I(my_by_block.size()))>(my_total.size(), nblocks)),
        my_ngroups(ngroups),
        my_nblocks(nblocks)
    {
        for (decltype(I(nblocks)) b = 0; b < nblocks; ++b) {
            for (decltype(I(ngroups)) g1 = 1; g1 < ngroups; ++g1) {
                const auto w1 = combo_weights[sanisizer::nd_offset<std::size_t>(g1, ngroups, b)];
                for (decltype(I(g1)) g2 = 0; g2 < g1; ++g2) {
                    const Stat_ combined = w1 * combo_weights[sanisizer::nd_offset<std::size_t>(g2, ngroups, b)];

                    // Storing it as a 3D array where the blocks are the fastest changing, 
                    // and then the two groups are the next fastest changing.
                    const auto out_offset1 = sanisizer::nd_offset<std::size_t>(g2, ngroups, g1);
                    my_by_block[sanisizer::nd_offset<std::size_t>(b, nblocks, out_offset1)] = combined;
                    my_by_block[sanisizer::nd_offset<std::size_t>(b, nblocks, g1, ngroups, g2)] = combined;

                    my_total[out_offset1] += combined;
                }
            }
        }

        // Filling the other side, for completeness.
        for (decltype(I(ngroups)) g1 = 1; g1 < ngroups; ++g1) {
            for (decltype(I(g1)) g2 = 0; g2 < g1; ++g2) {
                my_total[sanisizer::nd_offset<std::size_t>(g1, ngroups, g2)] = my_total[sanisizer::nd_offset<std::size_t>(g2, ngroups, g1)];
            }
        }
    }

public:
    std::pair<const Stat_*, Stat_> get(const std::size_t g1, const std::size_t g2) const {
        const auto offset = sanisizer::nd_offset<std::size_t>(g2, my_ngroups, g1);
        return std::make_pair(
            my_by_block.data() + offset * my_nblocks,
            my_total[offset]
        );
    }

private:
    std::vector<Stat_> my_total;
    std::vector<Stat_> my_by_block;
    std::size_t my_ngroups, my_nblocks;
};

template<typename Stat_>
class BlockAverageInfo {
public:
    BlockAverageInfo() = default;
    BlockAverageInfo(std::vector<Stat_> combo_weights) : my_combo_weights(std::move(combo_weights)) {}
    BlockAverageInfo(const double quantile) : my_quantile(quantile) {}

private:
    std::optional<std::vector<Stat_> > my_combo_weights;
    double my_quantile = 0;

public:
    bool use_mean() const {
        return my_combo_weights.has_value();
    }

    const std::vector<Stat_>& combo_weights() const {
        return *my_combo_weights;
    }

    double quantile() const {
        return my_quantile;
    }
};

}

}

#endif
