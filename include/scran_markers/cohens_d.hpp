#ifndef SCRAN_MARKERS_COHEN_D_HPP
#define SCRAN_MARKERS_COHEN_D_HPP

#include <vector>
#include <limits>
#include <cmath>
#include <type_traits>

namespace scran_markers {

namespace internal {

template<typename Stat_>
Stat_ compute_cohens_d(Stat_ m1, Stat_ m2, Stat_ sd, Stat_ threshold) {
    if (std::isnan(sd)) {
        return std::numeric_limits<Stat_>::quiet_NaN();
    } 
    
    Stat_ delta = m1 - m2 - threshold;
    if (sd == 0 && delta == 0) {
        return 0;
    } else if (sd == 0) {
        if (delta > 0) {
            return std::numeric_limits<Stat_>::infinity();
        } else {
            return -std::numeric_limits<Stat_>::infinity();
        }
    } else {
        return delta / sd;
    }
}

template<typename Stat_>
Stat_ cohen_denominator(Stat_ left_var, Stat_ right_var) {
    if (std::isnan(left_var) && std::isnan(right_var)) {
        return std::numeric_limits<Stat_>::quiet_NaN();
    } else if (std::isnan(left_var)) {
        return std::sqrt(right_var);
    } else if (std::isnan(right_var)) {
        return std::sqrt(left_var);
    } else {
        return std::sqrt((left_var + right_var)/2);
    }
}

// 'values' and 'combo_weights' are expected to be 'ngroups * nblocks' arrays
// where groups are the faster-changing dimension and the blocks are slower.
template<typename Stat_, typename Weight_, class Output_>
void compute_pairwise_cohens_d_internal(size_t g1, size_t g2, const Stat_* means, const Stat_* vars, size_t ngroups, size_t nblocks, const Weight_* combo_weights, Stat_ threshold, Output_& output) {
    Stat_ total_weight = 0;
    constexpr bool do_both_sides = !std::is_same<Stat_, Output_>::value;

    size_t offset1 = g1, offset2 = g2; // no need to cast, everything's already a size_t.
    for (size_t b = 0; b < nblocks; ++b, offset1 += ngroups, offset2 += ngroups) {
        auto left_weight = combo_weights[offset1];
        if (!left_weight) {
            continue;
        }

        auto right_weight = combo_weights[offset2];
        if (!right_weight) {
            continue;
        }

        auto left_var = vars[offset1];
        auto right_var = vars[offset2];
        Stat_ denom = cohen_denominator(left_var, right_var);
        if (std::isnan(denom)) {
            continue;
        }

        Stat_ weight = static_cast<Stat_>(left_weight) * static_cast<Stat_>(right_weight);
        total_weight += weight;

        auto left_mean = means[offset1];
        auto right_mean = means[offset2]; 
        Stat_ extra = compute_cohens_d(left_mean, right_mean, denom, threshold) * weight;

        if constexpr(do_both_sides) {
            output.first += extra;
            if (threshold) {
                output.second += compute_cohens_d(right_mean, left_mean, denom, threshold) * weight;
            }
        } else {
            output += extra;
        }
    }

    if constexpr(do_both_sides) {
        if (total_weight) {
            output.first /= total_weight;
            if (threshold) {
                output.second /= total_weight;
            } else {
                output.second = -output.first;
            }
        } else {
            output.first = std::numeric_limits<Stat_>::quiet_NaN();
            output.second = std::numeric_limits<Stat_>::quiet_NaN();
        }
    } else {
        if (total_weight) {
            output /= total_weight;
        } else {
            output = std::numeric_limits<Stat_>::quiet_NaN();
        }
    }
}

template<typename Stat_, typename Weight_>
Stat_ compute_pairwise_cohens_d_one_sided(size_t g1, size_t g2, const Stat_* means, const Stat_* vars, size_t ngroups, size_t nblocks, const Weight_* combo_weights, Stat_ threshold) {
    Stat_ output = 0;
    compute_pairwise_cohens_d_internal(g1, g2, means, vars, ngroups, nblocks, combo_weights, threshold, output);
    return output;
}

template<typename Stat_, typename Weight_>
std::pair<Stat_, Stat_> compute_pairwise_cohens_d_two_sided(size_t g1, size_t g2, const Stat_* means, const Stat_* vars, size_t ngroups, size_t nblocks, const Weight_* combo_weights, Stat_ threshold) {
    std::pair<Stat_, Stat_> output(0, 0);
    compute_pairwise_cohens_d_internal(g1, g2, means, vars, ngroups, nblocks, combo_weights, threshold, output);
    return output;
}

template<typename Stat_, typename Weight_>
void compute_pairwise_cohens_d(const Stat_* means, const Stat_* vars, size_t ngroups, size_t nblocks, const Weight_* combo_weights, Stat_ threshold, Stat_* output) {
    for (size_t g1 = 1; g1 < ngroups; ++g1) { // might as well skip 0, the inner loop doesn't do anything.
        for (size_t g2 = 0; g2 < g1; ++g2) {
            auto tmp = compute_pairwise_cohens_d_two_sided(g1, g2, means, vars, ngroups, nblocks, combo_weights, threshold);
            output[g1 * ngroups + g2] = tmp.first;
            output[g2 * ngroups + g1] = tmp.second;
        }
    }
}

}

}

#endif
