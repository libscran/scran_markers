#ifndef SCRAN_MARKERS_UTILS_HPP
#define SCRAN_MARKERS_UTILS_HPP

#include <type_traits>

#include "sanisizer/sanisizer.hpp"

namespace scran_markers {

template<typename Input_>
using I = std::remove_cv_t<std::remove_reference_t<Input_> >;

template<typename Group_>
std::size_t total_groups(const Group_* group, const std::size_t n) {
    if (n == 0) {
        return 0;
    } else {
        return sanisizer::sum<std::size_t>(*std::max_element(group, group + n), 1);
    }
}

template<typename Group_, typename Count_>
std::vector<Count_> tabulate_groups(const Group_* group, const Count_ n) {
    const auto ngroups = total_groups(group, n);
    auto output = sanisizer::create<std::vector<Count_> >(ngroups);
    for (Count_ i = 0; i < n; ++i) {
        output[group[i]] += 1;
    }
    return output;
}

}

#endif
