#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <algorithm>
#include <sstream>
#include <numeric>
#include <string>
#include <limits>
#include <cmath>

#include <aff3ct.hpp>

#include "Decoder_polar_SCL_fast_sys.hpp"

namespace aff3ct
{
namespace module
{
template <typename R> inline R              sat_m(const R           m) { return m; }
template <          > inline signed char    sat_m(const signed char m) { return tools::saturate<signed char>(m, -128, 63); }

template <typename R>
inline void normalize_scl_metrics(std::vector<R> &metrics, const int L){}

template <>
inline void normalize_scl_metrics(std::vector<short> &metrics, const int L)
{
    auto min = *std::min_element(metrics.begin(), metrics.begin() + L);

    const auto norm = std::numeric_limits<short>::min() - min;

    for (auto i = 0; i < L; i++)
        metrics[i] += norm;
}

template <>
inline void normalize_scl_metrics(std::vector<signed char> &metrics, const int L)
{
    
}
}
}