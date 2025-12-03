#ifndef SCRAN_MARKERS_UTILS_HPP
#define SCRAN_MARKERS_UTILS_HPP

#include <type_traits>

namespace scran_markers {

template<typename Input_>
using I = std::remove_cv_t<std::remove_reference_t<Input_> >;

}

#endif
