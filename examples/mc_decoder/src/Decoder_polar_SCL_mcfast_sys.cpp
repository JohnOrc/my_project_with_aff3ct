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

    if (rev_depth == m)
    {
        
    }
}



}
}