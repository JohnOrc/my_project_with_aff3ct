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

#include "Decoder_polar_SCL_mcfast_sys.hpp"

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
    auto min = *std::min_element(metrics.begin(), metrics.begin() + L);

    const auto norm = std::numeric_limits<signed char>::min() - min;

    for (auto i = 0; i < L; i++)
        metrics[i] += norm;
}

template <typename B, typename R, class API_polar>
Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::Decoder_polar_SCL_mcfast_sys(const int& K, const int& N, const int& L, const std::vector<bool>& frozen_bits)
: Decoder_SIHO<B,R>(K, N),
  m                ((int)std::log2(N)),
  L                (L),
  frozen_bits      (frozen_bits),
  polar_patterns   (frozen_bits,
                    {new tools::Pattern_polar_std,
                     new tools::Pattern_polar_r0,
                     new tools::Pattern_polar_r1,
                     new tools::Pattern_polar_r0_left,
                     new tools::Pattern_polar_rep_left,
                     new tools::Pattern_polar_rep,       // /!\ perf. degradation with REP nodes in fixed-point
                     new tools::Pattern_polar_spc(2,2)}, // /!\ perf. degradation with SPC nodes length > 4 (when L is big)
                    1,
                    2,
                    true),
  paths            (L),
  metrics          (L),
  l                (L, mipp::vector<R>(N + mipp::nElReg<R>())),
  s                (L, mipp::vector<B>(N + mipp::nElReg<B>())),
  metrics_vec      (3, std::vector<R>()),
  dup_count        (L, 0),
  bit_flips        (4 * L),
  is_even          (L),
  best_path        (0),
  n_active_paths   (1),
  n_array_ref      (L, std::vector<int>(m)),
  path_2_array     (L, std::vector<int>(m)),
  sorter           (N),
//sorter_simd      (N),
  best_idx         (L),
  l_tmp            (N)
{
    const std::string name = "Decoder_polar_SCL_mcfast_sys";
    this->set_name(name);
    this->set_n_frames_per_wave(API_polar::get_n_frames());

    static_assert(sizeof(B) == sizeof(R), "Sizes of the bits and reals have to be identical.");

    if (API_polar::get_n_frames() != 1)
        throw tools::invalid_argument(__FILE__, __LINE__, __func__, "The inter-frame API_polar is not supported.");

    if (!tools::is_power_of_2(this->N))
    {
        std::stringstream message;
        message << "'N' has to be a power of 2 ('N' = " << N << ").";
        throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
    }

    if (this->N != (int)frozen_bits.size())
    {
        std::stringstream message;
        message << "'frozen_bits.size()' has to be equal to 'N' ('frozen_bits.size()' = " << frozen_bits.size()
                << ", 'N' = " << N << ").";
        throw tools::length_error(__FILE__, __LINE__, __func__, message.str());
    }

    if (this->L <= 0 || !tools::is_power_of_2(this->L))
    {
        std::stringstream message;
        message << "'L' has to be a positive power of 2 ('L' = " << L << ").";
        throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
    }

    if (this->N < mipp::nElReg<R>() * 2)
    {
        std::stringstream message;
        message << "'N' has to be equal or greater than 'mipp::nElReg<R>()' * 2 ('N' = " << N
                << ", 'mipp::nElReg<R>()' = " << mipp::nElReg<R>() << ").";
        throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
    }

    auto k = 0; for (auto i = 0; i < this->N; i++) if (frozen_bits[i] == 0) k++;
    if (this->K != k)
    {
        std::stringstream message;
        message << "The number of information bits in the frozen_bits in invalid  ('K' = " << K << ", 'k' = " << k << ").";
        throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
    }

    metrics_vec[0].resize(L * 2);
    metrics_vec[1].resize(L * 4);
    metrics_vec[2].resize((L <= 2 ? 4 : 8) * L);
}

template <typename B, typename R, class API_polar>
Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::Decoder_polar_SCL_mcfast_sys(const int& K, const int& N, const int& L, const std::vector<bool>& frozenbits,
                               const std::vector<tools::Pattern_polar_i*> &polar_patterns,
                               const int idx_r0, const int idx_r1)
: Decoder_SIHo<B,R>(K, N),
  m                (L),
  frozen_bits      (frozen_bits),
  polar_patterns   (frozen_bits, polar_patterns, idx_r0, idx_r1),
  paths            (L),
  metrics          (L),
  l                (L, mipp::vector<R>(N + mipp::nElReg<R>())),
  s                (L, mipp::vector<B>(N + mipp::nElReg<B>())),
  metrics_vec      (3, std::vector<R>()),
  dup_count        (L, 0),
  bit_flips        (4 * L),
  is_even          (L),
  best_path        (0),
  n_active_paths   (1),
  n_array_ref      (L, std::vector<int>(m)),
  path_2_array     (L, std::vector<int>(m)),
  sorter           (N),
  best_idx         (L),
  l_tmp            (N)
{
    const std::string name = "Decoder_polar_SCL_mcfast_sys";
    this->set_name(name);
    this->set_n_frames_per_wave(API_polar::get_n_frames());

    static_assert(sizeof(B) == sizeof(R), "Sizes of the bits and reals have to be identical.");

    if (API_polar::get_n_frames() != 1)
        throw tools::invalid_argument(__FILE__, __LINE__, __func__, "The inter-frame API_polar is not supported.");
    
    if (!tools::is_power_of_2(this->N))
    {
        std::stringstream message;
        message << "'N' has to be a power of 2 ('N' = " << N << ").";
        throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
    }

    if (this->N != (int)frozen_bits.size())
    {
        std::stringstream message;
        message << "'frozen_bits.size()' has to be equal to 'N' ('frozen_bits.size()' = " << frozen_bits.size()
                << ", 'N' = " << N << ").";
        throw tools::length_error(__FILE__, __LINE__, __func__, message.str());
    }

    if (this->L <= 0 || !tools::is_power_of_2(this->L))
    {
        std::stringstream message;
        message << "'L' has to be a positive power of 2 ('L' = " << L << ").";
        throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
    }

    if (this->N < mipp::nElReg<R>() * 2)
    {
        std::stringstream message;
        message << "'N' has to be equal or greater than 'mipp::nElReg<R>()' * 2 ('N' = " << N       
                << ", 'mipp::nElReg<R>()' = " << mipp::nElReg<R>() << ").";
        throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
    }

    auto k = 0; for (auto i = 0; i < this->N; i++) if (frozen_bits[i] == 0) k++;
    if (this->K != k)
    {
        std::stringstream message;
        message << "The number of information bits in the frozen_bits is invalid ('K' = " << K << ", 'k' = "
                << k << ").";
        throw tools.runtime_error(__FILE__, __LINE__, __func__, message.str());
    }

    metrics_vec[0].resize(L * 2);
    metrics_vec[1].resize(L * 4);
    metrics_vec[2].resize((L <= 2 ? 4 : 8) * L);
}

template <typename B, typename R, class API_polar>
Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::~Decoder_polar_SCL_mcfast_sys()
{
}

template <typename B, typename R, class API_polar>
Decoder_polar_SCL_mcfast_sys<B,R,API_polar>* Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::clone() const
{
    auto m = new Decoder_polar_SCL_mcfast_sys(*this);
    m->deep_copy(*this);
    return m;
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::set_frozen_bits(const std::vector<bool>& fb)
{
    aff3ct::tools::fb_assert(frozen_bits, this->K, this->N);
    std::copy(fb.begin(),  fb.end(), this->frozen_bits.begin());
    polar_patterns.set_frozen_bits(fb);
}

template <typename B, typename R, class API_polar>
const std::vector<bool>& Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::get_frozen_bits() const
{
    return this->frozen_bits;
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::init_buffers()
{
    metrics[0] = std::numeric_limits<R>::min();
    std::iota(path.begin(), paths.begin() + L, 0);

    n_active_paths = 1;

    std::fill(n_array_ref [0].begin(), n_array_ref [0].end(), 1);
    std::fill(path_2_array[0].begin(), path_2_array[0].end(), 0);

    for (auto i = 1; i < L; i++)
        std::fill(n_array_ref[i].begin(), n_array_ref[i].end(), 0);
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::_decode(const R *Y_N)
{
    int first_node_id = 0, off_l = 0, off_s = 0;
    recursive_decode(Y_N, off_l, off_s, m, first_node_id);
}

template <typename B, typename R, class API_polar>
int Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::_decode_siho(const R *Y_N, B *V_K, const size_t frame_id)
{
    if (!API_polar::isAligned(Y_N))
        throw tools::runtime_error(__FILE__, __LINE__, __func__, "'Y_N' is misaligned memory.");
    
    if (!API_polar::isAligned(V_K))
        throw tools::runtime_error(__FILE__, __LINE__, __func__, "'V_K' is misaligned memory.");

    this->init_buffers();

    this->_decode(Y_N);

    this->select_best_path(frame_id);

    this->_store(V_K);

    return 0;
}

template <typename B, typename R, class API_polar>
int Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::_decode_siho_cw(const R *Y_N, B *V_N, const size_t frame_id)
{
    if (!API_polar::isAligned(Y_N))
        throw tools::runtime_error(__FILE__, __LINE__, __func__, "'Y_N' is misaligned memory.");

    if (!API_polar::isAligned(V_N))
        throw tools::runtime_error(__FILE__, __LINE__, __func__, "'V_N' is misaligned memory.");
    
    this->init_buffers();

    this->_decode(Y_N);

    this->select_best_path(frame_id);

    this->_store_cw(V_N);

    return 0;
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::recursive_decode(const R *Y_N, const int off_l, const int off_s, const int rev_depth, int &node_id)
{
    const int n_elmts = 1 << rev_depth;
    const int n_elm_2 = n_elmts >> 1;
    const auto node_type = polar_patterns.get_node_type(node_id);

    const bool is_terminal_pattern = (node_type == tools::polar_node_t::RATE_0) ||
                                     (node_type == tools::polar_node_t::RATE_1) ||
                                     (node_type == tools::polar_node_t::REP)    ||
                                     (node_type == tools::polar_node_t::SPC);

    // root node
    if (rev_depth == m)
    {
        switch (node_type)
        {
            case tools::polar_node_t::STANDARD:
                API_polar::f(Y_N, Y_N + n_elm_2, l[0].data(), n_elm_2);
                break;
            case tools::polar_node_t::REP_LEFT:
                API_polar::f(Y_N, Y_N + n_elm_2, l[0].data(), n_elm_2);
                break;
            default:
                break;
        }

        recursive_decode(Y_N, off_l, off_s, rev_depth - 1, ++node_id); // recursive call left

        // g
        switch (node_type)
        {
            case tools::polar_node_t::STANDARD:
                for (auto i = 0; i < n_active_paths; i++)
                {
                    const auto path = paths[i];
                    const auto child = l[up_ref_array_idx(path, rev_depth - 1)].data();
                    API_polar::g (Y_N, Y_N + n_elm_2, s[path].data() + off_s, child, n_elm_2);
                }
                break;
            case tools::polar_node_t::RATE_0_LEFT:
                for (auto i = 0; i < n_active_paths; i++)
                {
                    const auto path = paths[i];
                    const auto child = l[up_ref_array_idx(path, rev_depth - 1)].data();
                    API_polar::g0(Y_N, Y_N + n_elm_2,                         child, n_elm_2);
                }
            case tools::polar_node_t::REP_LEFT:
                for (auto i = 0; i < n_active_paths; i++)
                {
                    const auto path = paths[i];
                    const auto child = l[up_ref_array_idx(path, rev_depth - 1)].data();
                    API_polar::gr(Y_N, Y_N + n_elm_2, s[path].data() + off_s, child, n_elm_2);
                }
                break;
            default:
                break;
        }

        recursive_decode(Y_N, off_l, off_s + n_elm_2, rev_depth - 1, ++node_id); // recursive call right

        // xor
        switch (node_type)
        {
            case tools::polar_node_t::STANDARD:
                for (auto i = 0; i < n_active_paths; i++)
                    API_polar::xo (s[paths[i]], off_s, off_s + n_elm_2, off_s, n_elm_2);
                break;
            case tools::polar_node_t::RATE_0_LEFT:
                for (auto i = 0; i < n_active_paths; i++)
                    API_polar::xo0(s[paths[i]],        off_s + n_elm_2, off_s, n_elm_2);
                break;
            case tools::polar_node_t::REP_LEFT:
                for (auto i = 0; i < n_active_paths; i++)
                    API_polar::xo (s[paths[i]], off_s, off_s + n_elm_2, off_s, n_elm_2);
                break;
            default:
                break;
        }
    }
    else if (!is_terminal_pattern && rev_depth) // other node (not root or leaf)
    {
        // f
        switch (node_type)
        {
            case tools::polar_node_t::STANDARD:
                for (auto i = 0; i < n_active_paths; i++)
                {
                    const auto path   = paths[i];
                    const auto parent = l[path_2_array     [path][rev_depth    ]].data();
                    const auto child  = l[up_ref_array_idx (path, rev_depth - 1)].data();
                    API_polar::f(parent + off_l, parent + off_l + n_elm_2, child + off_l + n_elmts, e_elm_2);
                }
                break;
            case tools::polar_node_t::REP_LEFT:
                for (auto i = 0; i < n_active_paths; i++)
                {
                    const auto path   = paths[i];
                    const auto parent = l[path_2_array     [path][rev_depth   ]].data();
                    const auto child  = l[up_ref_array_idx(path, rev_depth - 1)].data();
                    API_polar::f(parent + off_l, parent + off_l + n_elm_2, child + off_l + n_elmts, n_elm_2);
                }
                break;
            case tools::polar_node_t::RATE_0_LEFT:
                for (auto i = 0; i < n_active_paths && n_active_paths > 1; i++)
                {
                    const auto path   = paths[i];
                    const auto parent = l[path_2_array     [path][rev_depth   ]].data();
                    const auto child  = l[up_ref_array_idx(path, rev_depth - 1)].data();
                    API_polar::f(parent + off_l, parent + off_l + n_lem_2, child + off_l + n_elmts, n_elm_2);
                }
                break;
            default:
                break;
        }

        recursive_decode(Y_N, off_l + n_elmts, off_s, rev_depth - 1, ++node_id);

        // g
        switch (node_type)
        {
            case tools::polar_node_t::STANDARD:
                for (auto i = 0; i < n_active_paths; i++)
                {
                    const auto path   = paths[i];
                    const auto parent = l[path_2_array    [path][rev_depth    ]].data();
                    const auto child  = l[up_ref_array_idx(path, rev_depth - 1)].data();
                    API_polar::g(parent + off_l, parent + off_l + n_elm_2, s[path].data() + off_s, child + off_l + n_elmts, n_elm_2);
                }
                break;
            case tools::polar_node_t::RATE_0_LEFT:
                for (auto i = 0; i < n_active_paths; i++)
                {
                    const auto path   = paths[i];
                    const auto parent = l[path_2_array    [path][rev_depth    ]].data();
                    const auto child  = l[up_ref_array_idx(path, rev_depth - 1)].data();
                    API_polar::g0(parent + off_l, parent + off_l + n_elm_2, child + off_l + n_elmts, e_elm_2);
                }
                break;
            case tools::polar_node_t::REP_LEFT:
                for (auto i = 0; i < n_active_paths; i++)
                {
                    const auto path   = paths[i];
                    const auto parent = l[path_2_array    [path][rev_depth    ]].data();
                    const auto child  = l[up_ref_array_idx(path, rev_depth - 1)].data();
                    API_polar::gr(parent + off_l, parent + off_l + n_elm_2, s[path].data() + off_s, child + off_l + n_elmts, n_elm_2);
                }
                break;
            default:
                break;
        }
        
        recursive_decode(Y_N, off_l + n_elmts, off_s + n_elm_2, rev_depth - 1, ++node_id);

        // xor
        switch (node_type)
        {
            case tools::polar_node_t::STANDARD:
                for (auto i = 0; i < n_active_paths; i++)
                    API_polar::xo (s[paths[i]], off_s, off_s + n_elm_2, off_s, n_elm_2);
                break;
            case tools::polar_node_t::RATE_0_LEFT:
                for (auto i = 0; i < n_active_paths; i++)
                    API_polar::xo0(s[paths[i]], off_s + n_elm_2, off_s, n_elm_2);
                break;
            case tools::polar_node_t::REP_LEFT:
                for (auto i = 0; i < n_active_paths; i++)
                    API_polar::xo (s[paths[i]], off_s, off_s + n_elm_2, off_s, n_elm_2);
                break;
            default:
                break;
        }
    }
    else // leaf node
    {
        // h
        switch (node_type)
        {
            case tools::polar_node_t::RATE_0:   update_paths_r0 (rev_depth, off_l, off_s, n_elmts);
                break;
            case tools::polar_node_t::REP:      update_paths_rep(rev_depth, off_l, off_s, n_elmts);
                break;
            case tools::polar_node_t::RATE_1:   update_paths_r1 (rev_depth, off_l, off_s, n_elmts);
                break;
            case tools::polar_node_t::SPC:      update_paths_spc(rev_depth, off_l, off_s, n_elmts);
                break;
            default:
                break;
        }

        normalize_scl_metrics<R>(this->metrics, this->L);
    }
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::_store(B *V_K) const
{
    tools::fb_extract(this->polar_patterns.get_leaves_pattern_types(), this->s[best_path].data(), V_K);
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::_store_cw(B *V_N) const
{
    std::copy(this->s[best_path].data(), this->s[best_path].data() + this->N, V_N);
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::update_paths_r0(const int r_d, const int off_l, const int off_s, const int n_elmts)
{
    if (n_active_paths > 1)
        for (auto i = 0; i < n_active_paths; i++)
        {
            const auto path  = paths[i];
            const auto array = path_2_array[path][r_d];

            auto pen = (R)0;
            for (auto j = 0; j < n_elmts; j++)
                pen = sat_m<R>(pen + sat_m<R>(-std::min((R)l[array][off_l + j], (R)0)));
            metrics[path] = sat_m<R>(metrics[path] + pen);
        }
    
    for (auto i = 0; i < n_active_paths; i++)
        API_polar::h0(s[paths[i]], off_s, n_elmts);
}

template <typename B, typename R, class API_polar>
template <int REV_D, int N_ELMTS>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::update_paths_r0(const int off_l, const int off_s)
{
    if (n_active_paths > 1)
        for (auto i = 0; i < n_active_paths; i++)
        {
            const auto path  = paths[i];
            const auto array = path_2_array[path][REV_D];

            auto pen = (R)0;
            for (auto j = 0; j < N_ELMTS; j++)
                pen = sat_m<R>(pen + sat_m<R>(-std::min((R)l[array][off_l + j], (R)0)));
            metrics[path] = sat_m<R>(metrics[path] + pen);
        }
    
    for (auto i = 0; i < n_active_paths; i++)
        API_polar::template h0<N_ELMTS>(s[paths[i]], off_s, N_ELMTS);
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::update_paths_r1(const int r_d, const int off_l, const int off_s, const int n_elmts)
{
    if (r_d == 0)
        update_paths_rep(r_d, off_l, off_s, n_elmts);
    else
    {
        // generate the candidates with Chase-II algorithm
        if (n_elmts == 2)
        {
            for (auto i = 0; i < n_active_paths; i++)
            {
                const auto path  = paths[i];
                const auto array = path_2_array[path][r_d];

                bit_flips[2 * path + 0] = 0;
                bit_flips[2 * path + 1] = 1;

                const auto pen0 = sat_m<R>(std::abs(l[array][off_l + bit_flips[2 * path + 0]]));
                const auto pen1 = sat_m<R>(std::abs(l[array][off_l + bit_flips[2 * path + 1]]));

                metrics_vec[1][4 * path + 0] =          metrics       [    path    ];
                metrics_vec[1][4 * path + 1] = sat_m<R>(metrics       [    path    ] + pen0);
                metrics_vec[1][4 * path + 2] = sat_m<R>(metrics       [    path    ] + pen1);  
                metrics_vec[1][4 * path + 3] = sat_m<R>(metrics_vec[1][4 * path + 1] + pen1);
            }
        }
        else
        {
            for (auto i = 0; i < n_active_paths; i++)
            {
                const auto path  = paths[i];
                const auto array = path_2_array[path][r_d];

                for (auto i = 0; i < n_elmts; i++) 
                    l_tmp[i] = std::abs(l[array][off_l + i]);
                
                sorter.partial_sort_destructive(l_tmp.data(), best_idx, n_elmts, 2);
                bit_flips[2 * path + 0] = best_idx[0];
                bit_flips[2 * path + 1] = best_idx[1];

                const auto pen0 = sat_m<R>(std::abs(l[array][off_l + bit_flips[2 * path + 0]]));
                const auto pen1 = sat_m<R>(std::abs(l[array][off_l + bit_flips[2 * path + 1]]));

                metrics_vec[1][4 * path + 0] =          metrics       [    path    ];
                metrics_vec[1][4 * path + 1] = sat_m<R>(metrics       [    path    ] + pen0);
                metrics_vec[1][4 * path + 2] = sat_m<R>(metrics       [    path    ] + pen1);
                metrics_vec[1][4 * path + 3] = sat_m<R>(metrics_vec[1][4 * path + 1] + pen1);
            }
        }
        for (auto i = n_active_paths; i < L; i++)
            for (auto j = 0; j < 4; j++)
                metrics_vec[1][4 * paths[i] + j] = std::numeric_limits<R>::max();

        const auto n_list = (n_active_paths * 4 >= L) ? L : n_active_paths * 4;
        sorter.partial_sort(metrics_vec[1].data(), best_idx, L * 4, n_list);

        for (auto i = 0; i < n_list; i++)
            dup_count[best_idx[i] / 4]++;
        
        erase_bad_paths();

        for (auto i = 0; i < n_list; i++)
        {
            const auto path  = best_idx[i] / 4;
            const auto dup   = best_idx[i] % 4;
            const auto array = path_2_array[path][r_d];

            API_polar::h(s[path], l[array], off_l, off_s, n_elmts);

            const auto new_path = (dup_count[path] > 1) ? duplicate_tree(path, off_l, off_s, n_elmts) : path;
            flip_bits_r1(path, new_path, dup, off_s, n_elmts);
            metrics[new_path] = metrics_vec[1][best_idx[i]];

            dup_count[path]--;
        }
    }
}

template <typename B, typename R, class API_polar>
template <int REV_D, int N_ELMTS>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::update_paths_r1(const int off_l, const int off_s)
{
    if (REV_D == 0)
        update_paths_rep<REV_D, N_ELMTS>(off_l, off_s);
    else
    {
        // generate the candidates with the Chase-II algorithm
        if (N_ELMTS == 2)
        {
            for (auto i = 0; i < n_active_paths; i++)
            {
                const auto path  = paths[i];
                const auto array = path_2_array[path][REV_D];

                bit_flips[2 * path + 0] = 0;
                bit_flips[2 * path + 1] = 1;

                const auto pen0 = sat_m<R>(std::abs(l[array][off_l + bit_flips[2 * path + 0]]));
                const auto pen1 = sat_m<R>(std::abs(l[array][off_l + bit_flips[2 * path + 1]]));

                metrics_vec[1][4 * path + 0] =          metrics       [    path    ];
                metrics_vec[1][4 * path + 1] = sat_m<R>(metrics       [    path    ] + pen0);
                metrics_vec[1][4 * path + 2] = sat_m<R>(metrics       [    path    ] + pen1);
                metrics_vec[1][4 * path + 3] = sat_m<R>(metrics_vec[1][4 * path + 1] + pen1);
            }
        }
        else
        {
            for (auto i = 0; i < n_active_paths; i++)
            {
                const auto path  = paths[i];
                const auto array = path_2_array[path][REV_D];

                for (auto i = 0; i < N_ELMTS; i++)
                    l_tmp[i] = std::abs(l[array][off_l + i]);
                sorter.partial_sort_destructive(l_tmp.data() + off_l, best_idx, N_ELMTS, 2);

                bit_flips[2 * path + 0] = best_idx[0];
                bit_flips[2 * path + 1] = best_idx[1];

                const auto pen0 = sat_m<R>(std::abs(l[array][off_l + bit_flips[2 * path + 0]]));
                const auto pen1 = sat_m<R>(std::abs(l[array][off_l + bit_flips[2 * path + 1]]));

                metrics_vec[1][4 * path + 0] =          metrics       [    path    ];
                metrics_vec[1][4 * path + 1] = sat_m<R>(metrics       [    path    ] + pen0);
                metrics_vec[1][4 * path + 2] = sat_m<R>(metrics       [    path    ] + pen1);
                metrics_vec[1][4 * path + 3] = sat_m<R>(metrics_vec[1][4 * path + 1] + pen1);
            }
        }
        for (auto i = n_active_paths; i < L; i++)
            for (auto j = 0; j < 4; j++)
                metrics_vec[1][4 * paths[i] + j] = std::numeric_limits<R>max();

        // L first of the lists are the L best paths
        const auto n_list = (n_active_paths * 4 >= L) ? L : n_active_paths * 4;
        sorter.partial_sort(metrics_vec[1].data(), best_idx, L * 4, n_list);

        // count the number of duplications per path
        for (auto i = 0; i < n_list; i++)
            dup_count[best_idx[i] / 4]++;

        erase_bad_paths();

        for (auto i = 0; i < n_list; i++)
        {
            const auto path  = best_idx[i] / 4;
            const auto dup   = best_idx[i] % 4;
            const auto array = path_2_array[path][REV_D];

            API_polar::template h<N_ELMTS>(s[path], l[array], off_l, off_s, N_ELMTS);

            const auto new_path = (dup_count[path] > 1) ? duplicate_tree(path, off_l, off_s, N_ELMTS) : path;
            flip_bits_r1(path, new_path, dup, off_s, N_ELMTS);
            metrics[new_path] = metrics_vec[1][best_idx[i]];

            dup_count[path]--;
        }
    }
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::flip_bits_r1(const int old_path, const int new_path, const int dup, const int off_s, const int n_elmts)
{
    
}



}
}