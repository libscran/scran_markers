#ifndef SCRAN_MARKERS_AUC_HPP
#define SCRAN_MARKERS_AUC_HPP

#include <vector>
#include <algorithm>

namespace scran_markers {

namespace internal {

template<typename Value_, typename Group_, typename Output_>
struct AucWorkspace {
    AucWorkspace(size_t ngroups, Output_* buffer) : less_than(ngroups), equal(ngroups), outputs(ngroups) {
        output_start = buffer;
        for (size_t i = 0; i < ngroups; ++i) {
            outputs[i] = buffer;
            buffer += ngroups;
        }
        output_end = buffer;
    }

    std::vector<std::pair<Value_, Group_> > paired;

    std::vector<Output_> less_than, equal;

    std::vector<Output_*> outputs;

    Output_* output_start;

    Output_* output_end;
};

template<typename Value_, typename Group_, typename Output_>
void prepare_auc_workspace(AucWorkspace<Value_, Group_, Output_>& work) {
    std::sort(work.paired.begin(), work.paired.end());
    std::fill(work.less_than.begin(), work.less_than.end(), 0);
    std::fill(work.equal.begin(), work.equal.end(), 0);
    std::fill(work.output_start, work.output_end, 0);
}

template<typename Value_, typename Group_, typename Count_, typename Output_>
void compute_pairwise_auc(AucWorkspace<Value_, Group_, Output_>& work, const std::vector<Count_>& num_zeros, const std::vector<Count_>& totals, bool normalize) {
    prepare_auc_workspace(work);

    auto& input = work.paired;
    auto& less_than = work.less_than;
    auto& equal = work.equal;
    auto& outputs = work.outputs;
    size_t ngroups = num_zeros.size();

    auto inner_loop = [&](size_t& pos) -> void {
        const auto& current = input[pos];

        ++pos;
        bool tied = false;
        while (pos != input.size() && input[pos].first == current.first) {
            tied = true;
            ++equal[input[pos].second];
            ++pos;
        }

        if (tied) {
            ++equal[current.second]; // contribution of the current (tied) observation.

            for (size_t l = 1; l < ngroups; ++l) { // starting from 1 as zero has no work for the g < l condition anyway.
                auto num_eq = equal[l];
                if (num_eq) {
                    auto outptr = outputs[l];
                    for (size_t g = 0; g < l; ++g) {
                        outptr[g] += num_eq * (less_than[g] + 0.5 * equal[g]);
                    }
                }
            }

            for (size_t l = 0; l < ngroups; ++l) {
                less_than[l] += equal[l];
                equal[l] = 0;
            }
        } else {
            auto outptr = outputs[current.second];
            for (Group_ g = 0; g < current.second; ++g) { // only computing the top half of the pairwise AUC matrix for now.
                outptr[g] += less_than[g];
            }
            ++less_than[current.second];
        }
    };

    size_t pos = 0;

    // Values < 0.
    while (pos != input.size() && input[pos].first < 0) {
        inner_loop(pos);
    }

    // Values == 0. 
    for (size_t l = 1; l < ngroups; ++l) { // starting from 1, see above.
        auto num_z = num_zeros[l];
        if (num_z) {
            auto outptr = outputs[l];
            for (size_t g = 0; g < l; ++g) {
                outptr[g] += num_z * (less_than[g] + 0.5 * num_zeros[g]);
            }
        }
    }

    for (size_t l = 0; l < ngroups; ++l) {
        less_than[l] += num_zeros[l];
    }

    // Values > 0. We expect that all values in 'input' are already non-zero;
    // any zeros should have been moved into the 'num_zeros' already.
    while (pos != input.size()) {
        inner_loop(pos);
    }

    // Filling in the other side.
    for (size_t l = 1; l < ngroups; ++l) { // starting from 1, see above.
        for (size_t g = 0; g < l; ++g) {
            Output_ prod = static_cast<Output_>(totals[l]) * static_cast<Output_>(totals[g]);
            outputs[g][l] = prod - outputs[l][g];

            if (normalize) {
                if (prod) {
                    outputs[l][g] /= prod;
                    outputs[g][l] /= prod;
                } else {
                    outputs[l][g] = std::numeric_limits<Output_>::quiet_NaN();
                    outputs[g][l] = std::numeric_limits<Output_>::quiet_NaN();
                }
            }
        }
    }
}

template<typename Value_, typename Group_, typename Count_, typename Output_, typename Threshold_>
void compute_pairwise_auc(AucWorkspace<Value_, Group_, Output_>& work, const std::vector<Count_>& num_zeros, const std::vector<Count_>& totals, Threshold_ threshold, bool normalize) {
    prepare_auc_workspace(work);

    auto& input = work.paired;
    auto& less_than = work.less_than;
    auto& equal = work.equal;
    auto& outputs = work.outputs;
    size_t ngroups = num_zeros.size();

    auto inner_loop = [&](size_t& pos, size_t& comp) -> void {
        const auto& current = input[pos];
        Threshold_ limit = current.first - threshold;

        // Hunting all entities less than the limit.
        while (comp != input.size() && static_cast<Threshold_>(input[comp].first) < limit) {
            ++less_than[input[comp].second];
            ++comp;
        }

        // Checking for ties with the limit.
        bool tied = false;
        while (comp != input.size() && static_cast<Threshold_>(input[comp].first) == limit) {
            tied = true;
            ++equal[input[comp].second];
            ++comp;
        }

        if (tied) {
            do {
                auto outptr = outputs[input[pos].second];
                for (size_t g = 0; g < ngroups; ++g) {
                    outptr[g] += less_than[g] + 0.5 * equal[g];
                }
                ++pos;
            } while (pos != input.size() && input[pos].first == current.first);

            for (size_t l = 0; l < ngroups; ++l) {
                less_than[l] += equal[l];
                equal[l] = 0;
            }
        } else {
            do {
                auto outptr = outputs[input[pos].second];
                for (size_t g = 0; g < ngroups; ++g) {
                    outptr[g] += less_than[g];
                }
                ++pos;
            } while (pos != input.size() && input[pos].first == current.first);
        }
    };

//    auto print_output = [&]() -> void {
//        std::cout << "Less than is: ";
//        for (size_t m = 0; m < ngroups; ++m) {
//            std::cout << " " << less_than[m];
//        }
//        std::cout << std::endl;
//        std::cout << "Equal is: ";
//        for (size_t m = 0; m < ngroups; ++m) {
//            std::cout << " " << equal[m];
//        }
//        std::cout << std::endl;
//        std::cout << "AUCs are:\n";
//        for (size_t n = 0; n < ngroups; ++n) {
//            std::cout << "  Group " << n << ": ";
//            for (size_t m = 0; m < ngroups; ++m) {
//                std::cout << " " << outputs[n][m];
//            }
//            std::cout << std::endl;
//        }
//    };

    size_t pos = 0, comp = 0;

    while (pos != input.size() && input[pos].first < 0) {
        inner_loop(pos, comp);
    }

    // Adding the contribution of zeros (in terms of the things they're greater than).
    // This effectively replicates the inner_loop but accounts for lots of zeros.
    // Again, we expect that all values in 'input' are already non-zero;
    // any zeros should have been moved into the 'num_zeros' already.
    {
        while (comp != input.size() && static_cast<Threshold_>(input[comp].first) < -threshold) {
            ++less_than[input[comp].second];
            ++comp;
        }

        for (size_t l = 0; l < ngroups; ++l) {
            auto num_z = num_zeros[l];
            if (num_z) {
                auto outptr = outputs[l];
                for (size_t g = 0; g < ngroups; ++g) {
                    outptr[g] += less_than[g] * num_z;
                }
            }
        }

        // Handling ties at the threshold boundary...
        bool tied = false;
        while (comp != input.size() && static_cast<Threshold_>(input[comp].first) == -threshold) {
            tied = true;
            ++equal[input[comp].second];
            ++comp;
        }

        if (tied) {
            for (size_t l = 0; l < ngroups; ++l) {
                auto num_z = num_zeros[l];
                if (num_z) {
                    auto outptr = outputs[l];
                    for (size_t g = 0; g < ngroups; ++g) {
                        outptr[g] += 0.5 * equal[g] * num_z;
                    }
                }
            }

            for (size_t l = 0; l < ngroups; ++l) {
                less_than[l] += equal[l];
                equal[l] = 0;
            }
        }

        // Or to each other, if the threshold is zero.
        if (threshold == 0) {
            for (size_t l = 0; l < ngroups; ++l) {
                auto num_z = num_zeros[l];
                if (num_z) {
                    auto outptr = outputs[l];
                    for (size_t g = 0; g < ngroups; ++g) {
                        outptr[g] += num_z * 0.5 * num_zeros[g];
                    }
                }
            }
        }
    }

    while (pos != input.size() && static_cast<Threshold_>(input[pos].first) < threshold) {
        inner_loop(pos, comp);
    }

    // Adding the contribution of zeros (in terms of the limit _being_ at zero + threshold)
    while (pos != input.size() && static_cast<Threshold_>(input[pos].first) == threshold) {
        auto outptr = outputs[input[pos].second];
        for (size_t g = 0; g < ngroups; ++g) {
            outptr[g] += less_than[g] + 0.5 * num_zeros[g];
        }
        ++pos;
    }

    for (size_t l = 0; l < ngroups; ++l) {
        less_than[l] += num_zeros[l];
    }

    while (pos != input.size()) {
        inner_loop(pos, comp);
    }

    // Dividing by the product of sizes.
    if (normalize) {
        for (size_t l = 0; l < ngroups; ++l) {
            for (size_t g = 0; g < ngroups; ++g) {
                Output_ prod = static_cast<Output_>(totals[l]) * static_cast<Output_>(totals[g]);
                if (prod) {
                    outputs[l][g] /= prod;
                } else {
                    outputs[l][g] = std::numeric_limits<Output_>::quiet_NaN();
                }
            }
        }
    }
}

}

}

#endif
