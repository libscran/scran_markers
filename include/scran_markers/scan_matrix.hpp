#ifndef SCRAN_MARKERS_SCAN_MATRIX_HPP 
#define SCRAN_MARKERS_SCAN_MATRIX_HPP

#include "cohens_d.hpp"
#include "auc.hpp"
#include "simple_diff.hpp"

#include <vector>
#include <cassert>
#include <algorithm>

#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"

namespace scran_markers {

namespace internal {

template<typename Value_, typename Group_, typename Index_, typename Stat_>
struct AucWorkspace {
    std::vector<PairedStore<Value_, Group_> > paired;
    std::vector<std::vector<Index_> > num_zeros;
    std::vector<std::vector<Index_> > totals;
    std::vector<Stat_> auc_buffer;
    std::vector<std::vector<Stat_> > block_scale;
    std::vector<Stat_> full_weight;
};

template<typename Value_, typename Group_, typename Index_, typename Stat_, typename Weight_>
void initialize_auc_workspace(
    AucWorkspace<Value_, Group_, Index_, Stat_>& work,
    size_t ngroups,
    size_t nblocks,
    const std::vector<Index_>& combo_size,
    const std::vector<Weight_>& combo_weight) 
{
    work.paired.resize(nblocks);

    work.num_zeros.reserve(nblocks);
    work.totals.reserve(nblocks);
    for (size_t g = 0; g < ngroups; ++g) {
        work.num_zeros.emplace_back(ngroups);
        work.totals.emplace_back(ngroups);
    }

    auto lsIt = combo_size.begin();
    for (size_t g = 0; g < ngroups; ++g) {
        for (size_t b = 0; b < nblocks; ++b, ++lsIt) {
            work.totals[b][g] = *lsIt;
        }
    }

    work.auc_buffer.resize(ngroups * ngroups);
    work.block_scale.reserve(nblocks);
    work.full_weight.resize(ngroups * ngroups)
    for (size_t b = 0; b < nblocks; ++b) {
        work.block_scale.emplace_back(ngroups * ngroups);

        for (size_t g1 = 0; g1 < ngroups; ++g1) {
            auto w1 = combo_weight[g1 * nblocks + b]; // all size_t's already, so no need to cast.

            for (size_t g2 = 0; g2 < g1; ++g2) {
                Stat_ block_denom = static_cast<Stat_>(work.totals[b][g1]) * static_cast<Stat_>(work.totals[b][g2]);
                if (block_denom == 0) {
                    continue;
                }

                Stat_ block_weight = w1 * combo_weight[g2 * nblocks + b];
                Stat_ block_scaling = block_denom / block_weight;

                auto offset1 = g1 * ngroups + g2;
                work.block_scale[b][offset1] = block_scaling;
                work.full_weight[offset1] += block_weight;

                auto offset2 = g2 * ngroups + g1;
                work.block_scale[b][offset2] = block_scaling;
                work.full_weight[offset2] += block_weight;
            }
        }
    }
}

template<typename Value_, typename Group_, typename Index_, typename Stat_, typename Threshold_>
void process_auc_for_rows(
    AucWorkspace<Value_, Group_, Index_, Stat_>& work,
    size_t ngroups,
    size_t nblocks,
    Threshold_ threshold,
    Stat_* output) 
{
    size_t ngroups2 = ngroups * ngroups;
    std::fill_n(output, ngroups2, 0);
    auto& auc_buffer = work.auc_buffer;

    for (size_t b = 0; b < nblocks; ++b) {
        auto& pr = work.paired[b];
        auto& nz = work.num_zeros[b];
        const auto& tt = work.totals[b];

        std::fill(auc_buffer.begin(), auc_buffer.end(), 0);
        if (threshold) {
            compute_pairwise_auc(pr, nz, tt, auc_buffer.data(), threshold, false);
        } else {
            compute_pairwise_auc(pr, nz, tt, auc_buffer.data(), false);
        }

        // Adding to the blocks.
        const auto& block_scale = work.block_scale[b];
        for (size_t g = 0; g < ngroups2; ++g) {
            if (block_scale[g]) {
                output[g] += auc_buffer[g] / block_scale[g];
            }
        }
    }

    size_t g = 0;
    for (size_t g1 = 0; g1 < ngroups; ++g1) {
        for (size_t g2 = 0; g2 < ngroups; ++g2, ++g) {
            auto& current = output[g];
            if (work.full_weight[g]) {
                current /= work.full_weight[g];
            } else if (g1 != g2) {
                current = std::numeric_limits<Stat_>::quiet_NaN();
            } // g1 == g2 gets a current = 0, which is technically wrong, but no one should be using the self-comparison effect size anyway.
        }
    }
}

template<bool single_block_, typename Value_, typename Index_, typename Group_, typename Block_, typename Stat_, typename Weight_>
void scan_matrix_by_row(
    const tatami::Matrix<Value_, Index_>& matrix, 
    size_t ngroups,
    const Group_* group,
    size_t nblocks, // should be equal to 1 if single_block_ = 1.
    const Block_* block, // ignored if single_block_ = true.
    size_t ncombos, // should be equal to ngroups if single_block_ = true. 
    const size_t* combinations, // ignored if single_block_ = true.
    std::vector<Stat_>& combo_means,
    std::vector<Stat_>& combo_vars,
    std::vector<Stat_>& combo_detected,
    Stat_* auc,
    const std::vector<Index_>& combo_size,
    const std::vector<Weight_>& combo_weights,
    int num_threads)
{
    Index_ NC = matrix.ncol();
    const auto* grouping = [&]{
        if constexpr(single_block_) {
            return group;
        } else {
            return combinations;
        }
    }();

    if constexpr(single_block_) {
        assert(ngroups == ncombos);
        assert(nblocks == 1);
    }

    tatami::parallelize([&](size_t, Index_ start, Index_ length) -> void {
        std::vector<Value_> vbuffer(NC);

        size_t offset = static_cast<size_t>(start) * ncombos; // cast to avoid overflow.
        auto mean_ptr = combo_means.data() + offset;
        auto var_ptr = combo_vars.data() + offset;
        auto det_ptr = combo_detected.data() + offset;

        // A vast array of AUC-related bits and pieces.
        size_t effect_shift = ngroups * ngroups;
        AucWorkspace<Value_, Group_, Index_, Stat_, Weight_> auc_work;
        auto auc_ptr = auc;
        if (auc_ptr) {
            auc_ptr += static_cast<size_t>(start) * effect_shift;
            initialize_auc_workspace(auc_work, ngroups, nblocks, combo_size, combo_weights);
        }

        if (matrix.is_sparse()) {
            std::vector<Index_> ibuffer(NC);
            auto ext = tatami::consecutive_extractor<true>(&matrix, true, start, length);
            std::vector<Index_> tmp_index(ncombos);

            for (size_t r = start, end = start + length; r < end; ++r) {
                auto range = ext->fetch(vbuffer.data(), ibuffer.data());
                tatami_stats::variances::direct(
                    range.value,
                    range.index,
                    range.number,
                    grouping,
                    ncombos,
                    combo_size.data(),
                    mean_ptr,
                    var_ptr,
                    tmp_index.data(),
                    /* skip_nan = */ false,
                    /* invalid_count = */ NULL
                );
                mean_ptr += ncombos;
                var_ptr += ncombos;

                for (Index_ i = 0; i < range.number; ++i) {
                    det_ptr[grouping[range.index[i]]] += (range.value[i] != 0);
                }
                det_ptr += ncombos;

                if (auc_ptr) {
                    for (auto& z : auc_work.num_zeros) {
                        std::fill(z.begin(), z.end(), 0);
                    }
                    for (auto& p : auc_work.paired) {
                        p.clear();
                    }

                    for (size_t c = 0; c < NC; ++c) {
                        auto b = [&]{
                            if constexpr(single_block_) {
                                return 0;
                            } else {
                                return block[c];
                            }
                        }();
                        auto g = group[c];
                        if (ptr[c]) {
                            auc_work.paired[b].emplace_back(ptr[c], g);
                        } else {
                            ++(auc_work.num_zeros[b][g]);
                        }
                    }

                    process_auc_for_rows(auc_work, ngroups, nblocks, threshold, auc_ptr);
                    auc_ptr += effect_shift;
                }
            }

        } else {
            auto ext = tatami::consecutive_extractor<false>(&matrix, true, start, length);

            for (size_t r = start, end = start + length; r < end; ++r) {
                auto ptr = ext->fetch(vbuffer.data());
                tatami_stats::variances::direct(
                    ptr,
                    NC,
                    grouping,
                    ncombos,
                    combo_size.data(),
                    mean_ptr,
                    var_ptr,
                    /* skip_nan = */ false,
                    /* invalid_count = */ NULL
                );
                mean_ptr += ncombos;
                var_ptr += ncombos;

                for (Index_ c = 0; c < NC; ++c) {
                    det_ptr[grouping[c]] += (ptr[c] != 0);
                }
                det_ptr += ncombos;

                if (auc_ptr) {
                    auto nzIt = auc_work.num_zeros.begin();
                    for (const auto& t : auc_work.totals) {
                        std::copy(t.begin(), t.end(), nzIt->begin());
                        ++nzIt;
                    }
                    for (auto& p : auc_work.paired) {
                        p.clear();
                    }

                    for (Index_ j = 0; j < range.number; ++j) {
                        if (range.value[j]) {
                            size_t c = range.index[j];
                            auto b = [&]{
                                if constexpr(single_block_) {
                                    return 0;
                                } else {
                                    return block[c];
                                }
                            }();
                            auto g = group[c];
                            auc_work.paired[b].emplace_back(range.value[j], g);
                            --(auc_work.num_zeros[b][g]);
                        }
                    }

                    process_auc_for_rows(auc_work, ngroups, nblocks, threshold, auc_ptr);
                    auc_ptr += effect_shift;
                }
            }
        }
    }, matrix.nrow(), num_threads);
}

template<typename Value_, typename Index_, typename Combination_, typename Stat_>
void scan_matrix_by_column(
    const tatami::Matrix<Value_, Index_>& matrix, 
    size_t ncombos,
    const Combination_* combinations,
    std::vector<Stat_>& combo_means,
    std::vector<Stat_>& combo_vars,
    std::vector<Stat_>& combo_detected)
{
    Idnex_ NC = matrix.ncol();
    tatami::parallelize([&](size_t, Index_ start, Index_ length) -> void {
        std::vector<Value_> vbuffer(length);

        // Using local buffers to avoid problems with false sharing.
        std::vector<std::vector<Stat_> > tmp_means, tmp_vars, tmp_det;
        tmp_means.reserve(ncombos);
        tmp_vars.reserve(ncombos);
        tmp_dets.reserve(ncombos);
        for (size_t co = 0; co < ncombos; ++co) {
            tmp_means.emplace_back(length);
            tmp_vars.emplace_back(length);
            tmp_dets.emplace_back(length);
        }

        if (matrix.is_sparse()) {
            std::vector<Index_> ibuffer(NC);
            auto ext = tatami::consecutive_extractor<true>(&matrix, false, static_cast<Index_>(0), NC, start, length);
            std::vector<tatami_stats::RunningSparse<Stat_, Value_, Index_> > runners;
            runners.reserve(ncombos);
            for (size_t co = 0; co < ncombos; ++co) {
                runners.emplace_back(length, tmp_means[co].data(), tmp_vars[co].data(), /* skip_nan = */ false, start);
            }

            for (Index_ c = 0; c < NC; ++c) {
                auto range = ext->fetch(vbuffer.data(), ibuffer.data());
                auto co = combinations[c];
                runners[co].add(range.value, range.index, range.number);

                auto& curdet = tmp_dets[co];
                for (Index_ i = 0; i < range.number; ++i) {
                    curdet[range.index[i] - start] += (range.value[i] != 0);
                }
            }

            for (auto& run : runners) {
                run.finish();
            }

        } else {
            auto ext = tatami::consecutive_extractor<false>(&matrix, false, static_cast<Index_>(0), NC, start, length);
            std::vector<tatami_stats::RunningDense<Stat_, Value_, Index_> > runners;
            runners.reserve(ncombos);
            for (size_t co = 0; co < ncombos; ++co) {
                runners.emplace_back(length, tmp_means[co].data(), tmp_vars[co].data(), /* skip_nan = */ false);
            }

            for (Index_ c = 0; c < NC; ++c) {
                auto ptr = ext->fetch(vbuffer.data());
                auto co = combinations[c];
                runners[co].add(ptr);

                auto& curdet = tmp_dets[co];
                for (Index_ r = 0; r < length; ++r) {
                    curdet[r] += (ptr[r] != 0);
                }
            }

            for (auto& run : runners) {
                run.finish();
            }
        }

        // Moving it all into the output buffers at the end.
        size_t offset = ncombos * static_cast<size_t>(start);
        auto mean_ptr = combo_means.data() + offset;
        auto var_ptr = combo_vars.data() + offset;
        auto det_ptr = combo_detected.data() + offset;

        for (Index_ r = 0; r < length; ++r) {
            for (size_t co = 0; co < ncombos; ++co) {
                *mean_ptr = tmp_means[co][r];
                *var_ptr = tmp_vars[co][r];
                *det_ptr = tmp_detected[co][r];
                ++mean_ptr;
                ++var_ptr;
                ++det_ptr;
            }
        }

    }, matrix.nrow(), num_threads);
}

}

}

#endif
