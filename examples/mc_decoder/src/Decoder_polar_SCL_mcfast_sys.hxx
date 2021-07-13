#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <algorithm>
#include <sstream>
#include <numeric>
#include <string>
#include <limits>
#include <cmath>
#include <queue>
#include <glog/logging.h>

#include "Decoder_polar_SCL_mcfast_sys.hpp"
#include "mc_func.hpp"


namespace aff3ct
{
namespace module
{

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
                    //  new tools::Pattern_polar_spc(2,2)}, // /!\ perf. degradation with SPC nodes length > 4 (when L is big)
					new tools::Pattern_polar_spc}, // /!\ perf. degradation with SPC nodes length > 4 (when L is big)
                    1,
                    2,
                    true),
  paths            (L),
  metrics          (L),
  l                (L, mipp::vector<R>(N + mipp::nElReg<R>())),
  s                (L, mipp::vector<B>(N + mipp::nElReg<B>())),
  metrics_vec      (3, std::vector<R>()),
  dup_count        (L, 0),
  bit_flips	       (L * L), // bits indices to be flipped * L paths
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

	std::cout << "mc decoder" << std::endl;

	static_assert(sizeof(B) == sizeof(R), "Sizes of the bits and reals have to be identical.");
//	static_assert(API_polar::get_n_frames() == 1, "The inter-frame API_polar is not supported.");

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
		throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
	}

	metrics_vec[0].resize(L * 2);
	switch (L)
	{
		case 2:
			r1_mc_size  = {2};
			spc_mc_size = {2};
			metrics_vec[1].resize(L * r1_mc_size[0]); // L * |C|
			metrics_vec[2].resize(L * spc_mc_size[0]);
			break;
		case 4:
			r1_mc_size  = {4, 5};
			spc_mc_size = {2, 5};
			metrics_vec[1].resize(L * r1_mc_size[1]); // L * |C|
			metrics_vec[2].resize(L * spc_mc_size[1]);
			break;
		case 8:
			r1_mc_size  = {4, 10, 13};
			spc_mc_size = {2, 8,  13};
			metrics_vec[1].resize(L * r1_mc_size[2]); // L * |C|
			metrics_vec[2].resize(L * spc_mc_size[2]);
			break;
		case 16:
			r1_mc_size  = {4, 16, 29, 36};
			spc_mc_size = {2, 8,  27, 36};
			metrics_vec[1].resize(L * r1_mc_size[3]); // L * |C|
			metrics_vec[2].resize(L * spc_mc_size[3]);
			break;
		case 32:
			r1_mc_size  = {4, 16, 59, 80, 95};
			spc_mc_size = {2, 8,  53, 78, 95};
			metrics_vec[1].resize(L * r1_mc_size[4]); // L * |C|
			metrics_vec[2].resize(L * spc_mc_size[4]);
			break;
		case 64:
			r1_mc_size  = {4, 16,  112, 186, 226, 258};
			spc_mc_size = {2, 8,   93,  182, 225, 258};
			metrics_vec[1].resize(L * r1_mc_size[5]); // L * |C|
			metrics_vec[2].resize(L * spc_mc_size[5]);
			break;
		default:
			metrics_vec[1].resize(L * 4); // Chase-II
			metrics_vec[2].resize((L <= 2 ? 4 : 8) * L);
			break;
	}	
}

template <typename B, typename R, class API_polar>
Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::Decoder_polar_SCL_mcfast_sys(const int& K, const int& N, const int& L, const std::vector<bool>& frozen_bits,
                             const std::vector<tools::Pattern_polar_i*> &polar_patterns,
                             const int idx_r0, const int idx_r1)
: Decoder_SIHO<B,R>(K, N),
  m                ((int)std::log2(N)),
  L                (L),
  frozen_bits      (frozen_bits),
  polar_patterns   (frozen_bits, polar_patterns, idx_r0, idx_r1),
  paths            (L),
  metrics          (L),
  l                (L, mipp::vector<R>(N + mipp::nElReg<R>())),
  s                (L, mipp::vector<B>(N + mipp::nElReg<B>())),
  metrics_vec      (3, std::vector<R>()),
  dup_count        (L, 0),
  bit_flips        (L * L),
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
//	static_assert(API_polar::get_n_frames() == 1, "The inter-frame API_polar is not supported.");

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
		throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
	}

	metrics_vec[0].resize(L * 2);
	switch (L)
	{
		case 2:
			r1_mc_size  = {2};
			spc_mc_size = {2};
			metrics_vec[1].resize(L * r1_mc_size[0]); // L * |C|
			metrics_vec[2].resize(L * spc_mc_size[0]);
			break;
		case 4:
			r1_mc_size  = {4, 5};
			spc_mc_size = {2, 5};
			metrics_vec[1].resize(L * r1_mc_size[1]); // L * |C|
			metrics_vec[2].resize(L * spc_mc_size[1]);
			break;
		case 8:
			r1_mc_size  = {4, 10, 13};
			spc_mc_size = {2, 8,  13};
			metrics_vec[1].resize(L * r1_mc_size[2]); // L * |C|
			metrics_vec[2].resize(L * spc_mc_size[2]);
			break;
		case 16:
			r1_mc_size  = {4, 16, 29, 36};
			spc_mc_size = {2, 8,  27, 36};
			metrics_vec[1].resize(L * r1_mc_size[3]); // L * |C|
			metrics_vec[2].resize(L * spc_mc_size[3]);
			break;
		case 32:
			r1_mc_size  = {4, 16, 59, 80, 95};
			spc_mc_size = {2, 8,  53, 78, 95};
			metrics_vec[1].resize(L * r1_mc_size[4]); // L * |C|
			metrics_vec[2].resize(L * spc_mc_size[4]);
			break;
		case 64:
			r1_mc_size  = {4, 16,  112, 186, 226, 258};
			spc_mc_size = {2, 8,   93,  182, 225, 258};
			metrics_vec[1].resize(L * r1_mc_size[5]); // L * |C|
			metrics_vec[2].resize(L * spc_mc_size[5]);
			break;
		default:
			metrics_vec[1].resize(L * 4); // Chase-II
			metrics_vec[2].resize((L <= 2 ? 4 : 8) * L);
			break;
	}
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
	std::copy(fb.begin(), fb.end(), this->frozen_bits.begin());
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
	std::iota(paths.begin(), paths.begin() + L, 0);

	n_active_paths = 1;

	// at the beginning, path 0 points to array 0
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

//	auto t_load = std::chrono::steady_clock::now(); // ----------------------------------------------------------- LOAD
	this->init_buffers();
//	auto d_load = std::chrono::steady_clock::now() - t_load;

//	auto t_decod = std::chrono::steady_clock::now(); // -------------------------------------------------------- DECODE
	this->_decode(Y_N);
	this->select_best_path(frame_id);
//	auto d_decod = std::chrono::steady_clock::now() - t_decod;

//	auto t_store = std::chrono::steady_clock::now(); // --------------------------------------------------------- STORE
	this->_store(V_K);
//	auto d_store = std::chrono::steady_clock::now() - t_store;

//	(*this)[dec::tsk::decode_siho].update_timer(dec::tm::decode_siho::load,   d_load);
//	(*this)[dec::tsk::decode_siho].update_timer(dec::tm::decode_siho::decode, d_decod);
//	(*this)[dec::tsk::decode_siho].update_timer(dec::tm::decode_siho::store,  d_store);

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

//	auto t_load = std::chrono::steady_clock::now(); // ----------------------------------------------------------- LOAD
	this->init_buffers();
//	auto d_load = std::chrono::steady_clock::now() - t_load;

//	auto t_decod = std::chrono::steady_clock::now(); // -------------------------------------------------------- DECODE
	this->_decode(Y_N);
	this->select_best_path(frame_id);
//	auto d_decod = std::chrono::steady_clock::now() - t_decod;

//	auto t_store = std::chrono::steady_clock::now(); // --------------------------------------------------------- STORE
	this->_store_cw(V_N);
//	auto d_store = std::chrono::steady_clock::now() - t_store;

//	(*this)[dec::tsk::decode_siho_cw].update_timer(dec::tm::decode_siho_cw::load,   d_load);
//	(*this)[dec::tsk::decode_siho_cw].update_timer(dec::tm::decode_siho_cw::decode, d_decod);
//	(*this)[dec::tsk::decode_siho_cw].update_timer(dec::tm::decode_siho_cw::store,  d_store);

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
		// f
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

		recursive_decode(Y_N, off_l, off_s, rev_depth -1, ++node_id); // recursive call left

		// g
		switch (node_type)
		{
			case tools::polar_node_t::STANDARD:
				for (auto i = 0; i < n_active_paths; i++)
				{
					const auto path  = paths[i];
					const auto child = l[up_ref_array_idx(path, rev_depth -1)].data();
					API_polar::g (Y_N, Y_N + n_elm_2, s[path].data() + off_s, child, n_elm_2);
				}
				break;
			case tools::polar_node_t::RATE_0_LEFT:
				for (auto i = 0; i < n_active_paths; i++)
				{
					const auto path  = paths[i];
					const auto child = l[up_ref_array_idx(path, rev_depth -1)].data();
					API_polar::g0(Y_N, Y_N + n_elm_2,                         child, n_elm_2);
				}
				break;
			case tools::polar_node_t::REP_LEFT:
				for (auto i = 0; i < n_active_paths; i++)
				{
					const auto path  = paths[i];
					const auto child = l[up_ref_array_idx(path, rev_depth -1)].data();
					API_polar::gr(Y_N, Y_N + n_elm_2, s[path].data() + off_s, child, n_elm_2);
				}
				break;
			default:
				break;
		}

		recursive_decode(Y_N, off_l, off_s + n_elm_2, rev_depth -1, ++node_id); // recursive call right

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
					const auto parent = l[path_2_array    [path][rev_depth   ]].data();
					const auto child  = l[up_ref_array_idx(path, rev_depth -1)].data();
					API_polar::f(parent + off_l, parent + off_l + n_elm_2, child + off_l + n_elmts, n_elm_2);
				}
				break;
			case tools::polar_node_t::REP_LEFT:
				for (auto i = 0; i < n_active_paths; i++)
				{
					const auto path   = paths[i];
					const auto parent = l[path_2_array    [path][rev_depth   ]].data();
					const auto child  = l[up_ref_array_idx(path, rev_depth -1)].data();
					API_polar::f(parent + off_l, parent + off_l + n_elm_2, child + off_l + n_elmts, n_elm_2);
				}
				break;
			case tools::polar_node_t::RATE_0_LEFT:
				for (auto i = 0; i < n_active_paths && n_active_paths > 1; i++)
				{
					const auto path   = paths[i];
					const auto parent = l[path_2_array    [path][rev_depth   ]].data();
					const auto child  = l[up_ref_array_idx(path, rev_depth -1)].data();
					API_polar::f(parent + off_l, parent + off_l + n_elm_2, child + off_l + n_elmts, n_elm_2);
				}
				break;
			default:
				break;
		}

		recursive_decode(Y_N, off_l + n_elmts, off_s, rev_depth -1, ++node_id); // recursive call left

		// g
		switch (node_type)
		{
			case tools::polar_node_t::STANDARD:
				for (auto i = 0; i < n_active_paths; i++)
				{
					const auto path   = paths[i];
					const auto parent = l[path_2_array    [path][rev_depth   ]].data();
					const auto child  = l[up_ref_array_idx(path, rev_depth -1)].data();
					API_polar::g (parent + off_l, parent + off_l + n_elm_2, s[path].data() + off_s, child + off_l + n_elmts, n_elm_2);
				}
				break;
			case tools::polar_node_t::RATE_0_LEFT:
				for (auto i = 0; i < n_active_paths; i++)
				{
					const auto path   = paths[i];
					const auto parent = l[path_2_array    [path][rev_depth   ]].data();
					const auto child  = l[up_ref_array_idx(path, rev_depth -1)].data();
					API_polar::g0(parent + off_l, parent + off_l + n_elm_2,                         child + off_l + n_elmts, n_elm_2);
				}
				break;
			case tools::polar_node_t::REP_LEFT:
				for (auto i = 0; i < n_active_paths; i++)
				{
					const auto path   = paths[i];
					const auto parent = l[path_2_array    [path][rev_depth   ]].data();
					const auto child  = l[up_ref_array_idx(path, rev_depth -1)].data();
					API_polar::gr(parent + off_l, parent + off_l + n_elm_2, s[path].data() + off_s, child + off_l + n_elmts, n_elm_2);
				}
				break;
			default:
				break;
		}

		recursive_decode(Y_N, off_l + n_elmts, off_s + n_elm_2, rev_depth -1, ++node_id); // recursive call right

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
	else // leaf node
	{
		// h
		switch (node_type)
		{
			case tools::polar_node_t::RATE_0: update_paths_r0 (rev_depth, off_l, off_s, n_elmts); break;
			case tools::polar_node_t::REP:    update_paths_rep(rev_depth, off_l, off_s, n_elmts); break;
			case tools::polar_node_t::RATE_1: update_paths_r1 (rev_depth, off_l, off_s, n_elmts); break;
			case tools::polar_node_t::SPC:    update_paths_spc(rev_depth, off_l, off_s, n_elmts); break;
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
				pen = sat_m<R>(pen + sat_m<R>(-std::min((R)l[array][off_l +j], (R)0)));
			metrics[path] = sat_m<R>(metrics[path] + pen); // add a penalty to the current path metric
		}

//	if (!polar_patterns.exist_node_type(polar_node_t::RATE_0_LEFT, r_d +1))
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
				pen = sat_m<R>(pen + sat_m<R>(-std::min((R)l[array][off_l +j], (R)0)));
			metrics[path] = sat_m<R>(metrics[path] + pen); // add a penalty to the current path metric
		}

//	if (!polar_patterns.exist_node_type(polar_node_t::RATE_0_LEFT, REV_D +1))
		for (auto i = 0; i < n_active_paths; i++)
			API_polar::template h0<N_ELMTS>(s[paths[i]], off_s, N_ELMTS);
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::update_paths_r1(const int r_d, const int off_l, const int off_s, const int n_elmts)
{

	if (r_d == 0)
	{
		update_paths_rep(r_d, off_l, off_s, n_elmts);
	}
	else
	{
		auto bits_num = (n_elmts >= L) ? L -1                  : n_elmts;
		auto indice_r = (n_elmts >= L) ? (int)std::log2(L) - 1 : (int)std::log2(n_elmts) - 1;

		DLOG(INFO) << "N_elmts = " << n_elmts << " bits_num = " << bits_num << " indice_r = " << indice_r;

		for (auto i = 0; i < n_active_paths; i++)
		{
			const auto path  = paths[i];
			const auto array = path_2_array[path][r_d];			

			for (auto j = 0; j < n_elmts; j++) l_tmp[j] = std::abs(l[array][off_l + j]);
			_partial_sort(l_tmp.data(), best_idx, n_elmts, bits_num);
			_update_r1_vec(path, array, n_elmts, bits_num, off_l);	
		}
		const auto n_list = (n_active_paths * r1_mc_size[indice_r] >= L) ? L : n_active_paths * r1_mc_size[indice_r];

		DLOG(INFO) << "n_list = " << n_list << " r1 size = " << r1_mc_size[indice_r];

		_partial_sort(metrics_vec[1].data(), best_idx, n_active_paths * r1_mc_size[indice_r], n_list);

		// count the number of duplications per path, count which old_path survive
		for (auto i = 0; i < n_list; i++)
			dup_count[best_idx[i] / r1_mc_size[indice_r]]++;

		// erase bad paths
		erase_bad_paths();

		for (auto i = 0; i < n_list; i++)
		{
			const auto path  = best_idx[i] / r1_mc_size[indice_r];
			const auto dup   = best_idx[i] % r1_mc_size[indice_r];
			const auto array = path_2_array[path][r_d];
			API_polar::h(s[path], l[array], off_l, off_s, n_elmts);
			const auto new_path = (dup_count[path] > 1) ? duplicate_tree(path, off_l, off_s, n_elmts) : path;
			flip_bits_r1(path, new_path, dup, off_s, n_elmts, bits_num);
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
	{
		update_paths_rep(REV_D, off_l, off_s, N_ELMTS);
	}
	else
	{
		int bits_num = (N_ELMTS >= L) ? L -1                  : N_ELMTS;
		int indice_r = (N_ELMTS >= L) ? (int)std::log2(L) - 1 : (int)std::log2(N_ELMTS) - 1;

		DLOG(INFO) << "N_elmts = " << N_ELMTS << " bits_num = " << bits_num << " indice_r = " << indice_r;

		for (auto i = 0; i < n_active_paths; i++)
		{
			const auto path  = paths[i];
			const auto array = path_2_array[path][REV_D];			

			for (auto j = 0; j < N_ELMTS; j++) l_tmp[j] = std::abs(l[array][off_l + j]);
			_partial_sort(l_tmp.data(), best_idx, N_ELMTS, bits_num);
			_update_r1_vec(path, array, N_ELMTS, bits_num, off_l);	
		}
		const auto n_list = (n_active_paths * r1_mc_size[indice_r] >= L) ? L : n_active_paths * r1_mc_size[indice_r];

		DLOG(INFO) << "n_list = " << n_list << " r1 size = " << r1_mc_size[indice_r];

		_partial_sort(metrics_vec[1].data(), best_idx, n_active_paths * r1_mc_size[indice_r], n_list);

		// count the number of duplications per path, count which old_path survive
		for (auto i = 0; i < n_list; i++)
			dup_count[best_idx[i] / r1_mc_size[indice_r]]++;

		// erase bad paths
		erase_bad_paths();

		for (auto i = 0; i < n_list; i++)
		{
			const auto path  = best_idx[i] / r1_mc_size[indice_r];
			const auto dup   = best_idx[i] % r1_mc_size[indice_r];
			const auto array = path_2_array[path][REV_D];
			API_polar::template h<N_ELMTS>(s[path], l[array], off_l, off_s, N_ELMTS);
			const auto new_path = (dup_count[path] > 1) ? duplicate_tree(path, off_l, off_s, N_ELMTS) : path;
			flip_bits_r1(path, new_path, dup, off_s, N_ELMTS, bits_num);
			metrics[new_path] = metrics_vec[1][best_idx[i]];
			dup_count[path]--;
		}
	}
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::update_paths_rep(const int r_d, const int off_l, const int off_s, const int n_elmts)
{
	constexpr B b = tools::bit_init<B>();

	// generate the two possible candidates
	for (auto i = 0; i < n_active_paths; i++)
	{
		const auto path  = paths[i];
		const auto array = path_2_array[path][r_d];

		auto pen0 = (R)0;
		auto pen1 = (R)0;
		for (auto j = 0; j < n_elmts; j++)
		{
			pen0 = sat_m<R>(pen0 + sat_m<R>(-std::min(l[array][off_l +j], (R)0)));
			pen1 = sat_m<R>(pen1 + sat_m<R>(+std::max(l[array][off_l +j], (R)0)));
		}
		metrics_vec[0][2 * path +0] = sat_m<R>(metrics[path] + pen0);
		metrics_vec[0][2 * path +1] = sat_m<R>(metrics[path] + pen1);
	}

	if (n_active_paths <= L / 2)
	{
		const auto n_active_paths_cpy = n_active_paths;
		for (auto i = 0; i < n_active_paths_cpy; i++)
		{
			const auto path = paths[i];
			const auto new_path = duplicate_tree(path, off_l, off_s, n_elmts);

			std::fill(s[    path].begin() + off_s, s[    path].begin() + off_s + n_elmts, 0);
			std::fill(s[new_path].begin() + off_s, s[new_path].begin() + off_s + n_elmts, b);

			metrics[    path] = metrics_vec[0][2 * path +0];
			metrics[new_path] = metrics_vec[0][2 * path +1];
		}
	}
	else // n_active_paths == L
	{
		// sort hypothetic metrics
		sorter.partial_sort(metrics_vec[0].data(), best_idx, L * 2, L);

		// count the number of duplications per path
		for (auto i = 0; i < L; i++)
			dup_count[best_idx[i] / 2]++;

		// erase bad paths
		erase_bad_paths();

		// duplicate paths
		for (auto path = 0; path < L; path++)
		{
			if (dup_count[path] == 1)
			{
				const auto comp = metrics_vec[0][2 * path +0] > metrics_vec[0][2 * path +1];
				std::fill(s[path].begin() + off_s, s[path].begin() + off_s + n_elmts, comp ? b : 0);

				metrics[path] = metrics_vec[0][2 * path + (comp ? 1 : 0)];
			}
			else if (dup_count[path] == 2)
			{
				const auto new_path = duplicate_tree(path, off_l, off_s, n_elmts);
				std::fill(s[    path].begin() + off_s, s[    path].begin() + off_s + n_elmts, 0);
				std::fill(s[new_path].begin() + off_s, s[new_path].begin() + off_s + n_elmts, b);

				metrics[    path] = metrics_vec[0][2 * path +0];
				metrics[new_path] = metrics_vec[0][2 * path +1];
			}

			dup_count[path] = 0;
		}
	}
}

template <typename B, typename R, class API_polar>
template <int REV_D, int N_ELMTS>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::update_paths_rep(const int off_l, const int off_s)
{
	constexpr B b = tools::bit_init<B>();

	// generate the two possible candidates
	for (auto i = 0; i < n_active_paths; i++)
	{
		const auto path  = paths[i];
		const auto array = path_2_array[path][REV_D];

		auto pen0 = (R)0;
		auto pen1 = (R)0;
		for (auto j = 0; j < N_ELMTS; j++)
		{
			pen0 = sat_m<R>(pen0 + sat_m<R>(-std::min(l[array][off_l +j], (R)0)));
			pen1 = sat_m<R>(pen1 + sat_m<R>(+std::max(l[array][off_l +j], (R)0)));
		}
		metrics_vec[0][2 * path +0] = sat_m<R>(metrics[path] + pen0);
		metrics_vec[0][2 * path +1] = sat_m<R>(metrics[path] + pen1);
	}

	if (n_active_paths <= L / 2)
	{
		const auto n_active_paths_cpy = n_active_paths;
		for (auto i = 0; i < n_active_paths_cpy; i++)
		{
			const auto path = paths[i];
			const auto new_path = duplicate_tree(path, off_l, off_s, N_ELMTS);

			std::fill(s[    path].begin() + off_s, s[    path].begin() + off_s + N_ELMTS, 0);
			std::fill(s[new_path].begin() + off_s, s[new_path].begin() + off_s + N_ELMTS, b);

			metrics[    path] = metrics_vec[0][2 * path +0];
			metrics[new_path] = metrics_vec[0][2 * path +1];
		}
	}
	else // n_active_paths == L
	{
		// sort hypothetic metrics
		sorter.partial_sort(metrics_vec[0].data(), best_idx, L * 2, L);

		// count the number of duplications per path
		for (auto i = 0; i < L; i++)
			dup_count[best_idx[i] / 2]++;

		// erase bad paths
		erase_bad_paths();

		// duplicate paths
		for (auto path = 0; path < L; path++)
		{
			if (dup_count[path] == 1)
			{
				const auto comp = metrics_vec[0][2 * path +0] > metrics_vec[0][2 * path +1];
				std::fill(s[path].begin() + off_s, s[path].begin() + off_s + N_ELMTS, comp ? b : 0);

				metrics[path] = metrics_vec[0][2 * path + (comp ? 1 : 0)];
			}
			else if (dup_count[path] == 2)
			{
				const auto new_path = duplicate_tree(path, off_l, off_s, N_ELMTS);
				std::fill(s[    path].begin() + off_s, s[    path].begin() + off_s + N_ELMTS, 0);
				std::fill(s[new_path].begin() + off_s, s[new_path].begin() + off_s + N_ELMTS, b);

				metrics[    path] = metrics_vec[0][2 * path +0];
				metrics[new_path] = metrics_vec[0][2 * path +1];
			}

			dup_count[path] = 0;
		}
	}
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::update_paths_spc(const int r_d, const int off_l, const int off_s, const int n_elmts)
{
	auto bits_num = (n_elmts >= L) ? L                     : n_elmts;
	auto indice_r = (n_elmts >= L) ? (int)std::log2(L) - 1 : (int)std::log2(n_elmts) - 1;

	DLOG(INFO) << "spclength = " << n_elmts;

	for (auto i = 0; i < n_active_paths; i++)
	{
		const auto path  = paths[i];
		const auto array = path_2_array[paths[i]][r_d];

		for (auto j = 0; j < n_elmts; j++) l_tmp[j] = std::abs(l[array][off_l + j]);
		_partial_sort(l_tmp.data(), best_idx, n_elmts, bits_num);
		_update_spc_vec(path, array, n_elmts, bits_num, off_l);
	}
	const auto n_list = (n_active_paths * spc_mc_size[indice_r] >= L) ? L : n_active_paths * spc_mc_size[indice_r];

	_partial_sort(metrics_vec[2].data(), best_idx, n_active_paths * spc_mc_size[indice_r], n_list);

	// count the number of duplications per path, count which old_path survive
	for (auto i = 0; i < n_list; i++)
		dup_count[best_idx[i] / spc_mc_size[indice_r]]++;

	// erase bad paths
	erase_bad_paths();

	for (auto i = 0; i < n_list; i++)
	{
		const auto path  = best_idx[i] / spc_mc_size[indice_r];
		const auto dup   = best_idx[i] % spc_mc_size[indice_r];
		const auto array = path_2_array[path][r_d];
		API_polar::h(s[path], l[array], off_l, off_s, n_elmts);
		const auto new_path = (dup_count[path] > 1) ? duplicate_tree(path, off_l, off_s, n_elmts) : path;
		flip_bits_spc(path, new_path, dup, off_s, n_elmts, bits_num);
		metrics[new_path] = metrics_vec[2][best_idx[i]];
		dup_count[path]--;
	}
}

template <typename B, typename R, class API_polar>
template <int REV_D, int N_ELMTS>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::update_paths_spc(const int off_l, const int off_s)
{
	auto bits_num = (N_ELMTS >= L) ? L                     : N_ELMTS;
	auto indice_r = (N_ELMTS >= L) ? (int)std::log2(L) - 1 : (int)std::log2(N_ELMTS) - 1;


	for (auto i = 0; i < n_active_paths; i++)
	{
		const auto path  = paths[i];
		const auto array = path_2_array[path][REV_D];			

		for (auto j = 0; j < N_ELMTS; j++) l_tmp[j] = std::abs(l[array][off_l + j]);
		_partial_sort(l_tmp.data(), best_idx, N_ELMTS, bits_num);
		_update_spc_vec(path, array, N_ELMTS, bits_num, off_l);	
	}
	const auto n_list = (n_active_paths * spc_mc_size[indice_r] >= L) ? L : n_active_paths * spc_mc_size[indice_r];


	_partial_sort(metrics_vec[2].data(), best_idx, n_active_paths * spc_mc_size[indice_r], n_list);

	// count the number of duplications per path, count which old_path survive
	for (auto i = 0; i < n_list; i++)
		dup_count[best_idx[i] / spc_mc_size[indice_r]]++;

	// erase bad paths
	erase_bad_paths();

	for (auto i = 0; i < n_list; i++)
	{
		const auto path  = best_idx[i] / spc_mc_size[indice_r];
		const auto dup   = best_idx[i] % spc_mc_size[indice_r];
		const auto array = path_2_array[path][REV_D];

		API_polar::template h<N_ELMTS>(s[path], l[array], off_l, off_s, N_ELMTS);

		const auto new_path = (dup_count[path] > 1) ? duplicate_tree(path, off_l, off_s, N_ELMTS) : path;
		flip_bits_spc(path, new_path, dup, off_s, N_ELMTS, bits_num);
		metrics[new_path] = metrics_vec[2][best_idx[i]];
		dup_count[path]--;
	}
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::flip_bits_spc(const int old_path, const int new_path, const int dup, const int off_s, const int n_elmts, const int bits_num)
{
	constexpr B b = tools::bit_init<B>();

	if (!is_even[old_path])
	{
		s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
	}

	switch(L)
	{
	case 2:
		switch( dup )
		{
		case 0:
			// nothing to do
			break;
		case 1:// e + a0 + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		default:
			throw tools::runtime_error(__FILE__, __LINE__, __func__, "Flip bits error on spc node.");
			break;
		}
		break;
	case 4:
		switch( dup )
		{
		case 0:
			// nothing to do
			break;
		case 1:// e + a0 + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		case 2:// e + a0 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 3:// e + a0 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 4:// e + a1 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		default:
			throw tools::runtime_error(__FILE__, __LINE__, __func__, "Flip bits error on spc node.");
			break;
		}
		break;
	case 8:
		switch( dup )
		{
		case 0:
			// nothing to do
			break;
		case 1:// e + a0 + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		case 2:// e + a0 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 3:// e + a0 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 4:// e + a1 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 5:// e + a1 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 6:// e + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 7:// e + a0 + a1 + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 8:// e + a0 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 9:// e + a0 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 10:// e + a0 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 11:// e + a0 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 12:// e + a1 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		default:
			throw tools::runtime_error(__FILE__, __LINE__, __func__, "Flip bits error on spc node.");
			break;
		}
		break;
	case 16:
		switch( dup )
		{
		case 0:
			// nothing to do
			break;
		case 1:// e + a0 + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		case 2:// e + a0 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 3:// e + a0 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 4:// e + a1 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 5:// e + a1 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 6:// e + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 7:// e + a0 + a1 + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 8:// e + a0 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 9:// e + a0 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 10:// e + a0 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 11:// e + a0 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 12:// e + a1 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 13:// e + a1 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 14:// e + a1 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 15:// e + a1 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 16:// e + a2 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 17:// e + a2 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 18:// e + a2 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 19:// e + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 20:// e + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 21:// e + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 22:// e + a0 + a1 + a2 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 23:// e + a0 + a1 + a2 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 24:// e + a0 + a1 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 25:// e + a0 + a2 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 26:// e + a1 + a2 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 27:// e + a0 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 28:// e + a0 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 29:// e + a0 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 30:// e + a0 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 31:// e + a0 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 32:// e + a0 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 33:// e + a0 + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 34:// e + a0 + a15 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +15]] = s[old_path][off_s + bit_flips[bits_num * old_path +15]] ? 0 : b;
			break;
		case 35:// e + a1 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		default:
			throw tools::runtime_error(__FILE__, __LINE__, __func__, "Flip bits error on spc node.");
			break;
		}
		break;
	case 32:
		switch( dup )
		{
		case 0:
			// nothing to do
			break;
		case 1:// e + a0 + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		case 2:// e + a0 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 3:// e + a0 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 4:// e + a1 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 5:// e + a1 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 6:// e + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 7:// e + a0 + a1 + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 8:// e + a0 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 9:// e + a0 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 10:// e + a0 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 11:// e + a0 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 12:// e + a1 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 13:// e + a1 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 14:// e + a1 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 15:// e + a1 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 16:// e + a2 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 17:// e + a2 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 18:// e + a2 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 19:// e + a2 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 20:// e + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 21:// e + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 22:// e + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 23:// e + a3 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 24:// e + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 25:// e + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 26:// e + a4 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 27:// e + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 28:// e + a5 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 29:// e + a6 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 30:// e + a0 + a1 + a2 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 31:// e + a0 + a1 + a2 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 32:// e + a0 + a1 + a2 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 33:// e + a0 + a1 + a2 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 34:// e + a0 + a1 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 35:// e + a0 + a1 + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 36:// e + a0 + a1 + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 37:// e + a0 + a1 + a3 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 38:// e + a0 + a1 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 39:// e + a0 + a1 + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 40:// e + a0 + a1 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 41:// e + a0 + a2 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 42:// e + a0 + a2 + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 43:// e + a0 + a2 + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 44:// e + a0 + a2 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 45:// e + a0 + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 46:// e + a1 + a2 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 47:// e + a1 + a2 + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 48:// e + a1 + a2 + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 49:// e + a1 + a2 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 50:// e + a1 + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 51:// e + a2 + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 52:// e + a0 + a1 + a2 + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 53:// e + a0 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 54:// e + a0 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 55:// e + a0 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 56:// e + a0 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 57:// e + a0 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 58:// e + a0 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 59:// e + a0 + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 60:// e + a0 + a15 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +15]] = s[old_path][off_s + bit_flips[bits_num * old_path +15]] ? 0 : b;
			break;
		case 61:// e + a1 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 62:// e + a1 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 63:// e + a1 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 64:// e + a1 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 65:// e + a1 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 66:// e + a1 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 67:// e + a1 + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 68:// e + a1 + a15 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +15]] = s[old_path][off_s + bit_flips[bits_num * old_path +15]] ? 0 : b;
			break;
		case 69:// e + a2 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 70:// e + a2 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 71:// e + a2 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 72:// e + a2 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 73:// e + a3 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 74:// e + a3 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 75:// e + a4 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 76:// e + a0 + a1 + a2 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 77:// e + a0 + a1 + a2 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 78:// e + a0 + a16 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +16]] = s[old_path][off_s + bit_flips[bits_num * old_path +16]] ? 0 : b;
			break;
		case 79:// e + a0 + a17 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +17]] = s[old_path][off_s + bit_flips[bits_num * old_path +17]] ? 0 : b;
			break;
		case 80:// e + a0 + a18 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +18]] = s[old_path][off_s + bit_flips[bits_num * old_path +18]] ? 0 : b;
			break;
		case 81:// e + a0 + a19 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +19]] = s[old_path][off_s + bit_flips[bits_num * old_path +19]] ? 0 : b;
			break;
		case 82:// e + a0 + a20 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +20]] = s[old_path][off_s + bit_flips[bits_num * old_path +20]] ? 0 : b;
			break;
		case 83:// e + a0 + a21 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +21]] = s[old_path][off_s + bit_flips[bits_num * old_path +21]] ? 0 : b;
			break;
		case 84:// e + a0 + a22 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +22]] = s[old_path][off_s + bit_flips[bits_num * old_path +22]] ? 0 : b;
			break;
		case 85:// e + a0 + a23 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +23]] = s[old_path][off_s + bit_flips[bits_num * old_path +23]] ? 0 : b;
			break;
		case 86:// e + a0 + a24 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +24]] = s[old_path][off_s + bit_flips[bits_num * old_path +24]] ? 0 : b;
			break;
		case 87:// e + a0 + a25 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +25]] = s[old_path][off_s + bit_flips[bits_num * old_path +25]] ? 0 : b;
			break;
		case 88:// e + a0 + a26 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +26]] = s[old_path][off_s + bit_flips[bits_num * old_path +26]] ? 0 : b;
			break;
		case 89:// e + a0 + a27 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +27]] = s[old_path][off_s + bit_flips[bits_num * old_path +27]] ? 0 : b;
			break;
		case 90:// e + a0 + a28 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +28]] = s[old_path][off_s + bit_flips[bits_num * old_path +28]] ? 0 : b;
			break;
		case 91:// e + a0 + a29 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +29]] = s[old_path][off_s + bit_flips[bits_num * old_path +29]] ? 0 : b;
			break;
		case 92:// e + a0 + a30 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +30]] = s[old_path][off_s + bit_flips[bits_num * old_path +30]] ? 0 : b;
			break;
		case 93:// e + a0 + a31 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +31]] = s[old_path][off_s + bit_flips[bits_num * old_path +31]] ? 0 : b;
			break;
		case 94:// e + a1 + a16 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +16]] = s[old_path][off_s + bit_flips[bits_num * old_path +16]] ? 0 : b;
			break;
		default:
			throw tools::runtime_error(__FILE__, __LINE__, __func__, "Flip bits error on spc node.");
			break;
		}
		break;
	case 64:
		switch( dup )
		{
		case 0:
			// nothing to do
			break;
		case 1:// e + a0 + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		case 2:// e + a0 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 3:// e + a0 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 4:// e + a1 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 5:// e + a1 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 6:// e + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 7:// e + a0 + a1 + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 8:// e + a0 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 9:// e + a0 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 10:// e + a0 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 11:// e + a0 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 12:// e + a1 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 13:// e + a1 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 14:// e + a1 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 15:// e + a1 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 16:// e + a2 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 17:// e + a2 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 18:// e + a2 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 19:// e + a2 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 20:// e + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 21:// e + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 22:// e + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 23:// e + a3 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 24:// e + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 25:// e + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 26:// e + a4 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 27:// e + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 28:// e + a5 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 29:// e + a6 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 30:// e + a0 + a1 + a2 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 31:// e + a0 + a1 + a2 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 32:// e + a0 + a1 + a2 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 33:// e + a0 + a1 + a2 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 34:// e + a0 + a1 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 35:// e + a0 + a1 + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 36:// e + a0 + a1 + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 37:// e + a0 + a1 + a3 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 38:// e + a0 + a1 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 39:// e + a0 + a1 + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 40:// e + a0 + a1 + a4 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 41:// e + a0 + a1 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 42:// e + a0 + a1 + a5 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 43:// e + a0 + a1 + a6 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 44:// e + a0 + a2 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 45:// e + a0 + a2 + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 46:// e + a0 + a2 + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 47:// e + a0 + a2 + a3 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 48:// e + a0 + a2 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 49:// e + a0 + a2 + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 50:// e + a0 + a2 + a4 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 51:// e + a0 + a2 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 52:// e + a0 + a2 + a5 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 53:// e + a0 + a2 + a6 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 54:// e + a0 + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 55:// e + a0 + a3 + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 56:// e + a0 + a3 + a4 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 57:// e + a0 + a3 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 58:// e + a0 + a3 + a5 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 59:// e + a0 + a3 + a6 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 60:// e + a0 + a4 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 61:// e + a0 + a4 + a5 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 62:// e + a0 + a4 + a6 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 63:// e + a0 + a5 + a6 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 64:// e + a1 + a2 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 65:// e + a1 + a2 + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 66:// e + a1 + a2 + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 67:// e + a1 + a2 + a3 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 68:// e + a1 + a2 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 69:// e + a1 + a2 + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 70:// e + a1 + a2 + a4 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 71:// e + a1 + a2 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 72:// e + a1 + a2 + a5 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 73:// e + a1 + a2 + a6 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 74:// e + a1 + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 75:// e + a1 + a3 + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 76:// e + a1 + a3 + a4 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 77:// e + a1 + a3 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 78:// e + a1 + a4 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 79:// e + a2 + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 80:// e + a2 + a3 + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 81:// e + a2 + a3 + a4 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 82:// e + a2 + a3 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 83:// e + a2 + a4 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 84:// e + a3 + a4 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 85:// e + a0 + a1 + a2 + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 86:// e + a0 + a1 + a2 + a3 + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 87:// e + a0 + a1 + a2 + a3 + a4 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 88:// e + a0 + a1 + a2 + a3 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 89:// e + a0 + a1 + a2 + a4 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 90:// e + a0 + a1 + a3 + a4 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 91:// e + a0 + a2 + a3 + a4 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 92:// e + a1 + a2 + a3 + a4 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 93:// e + a0 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 94:// e + a0 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 95:// e + a0 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 96:// e + a0 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 97:// e + a0 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 98:// e + a0 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 99:// e + a0 + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 100:// e + a0 + a15 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +15]] = s[old_path][off_s + bit_flips[bits_num * old_path +15]] ? 0 : b;
			break;
		case 101:// e + a1 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 102:// e + a1 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 103:// e + a1 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 104:// e + a1 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 105:// e + a1 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 106:// e + a1 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 107:// e + a1 + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 108:// e + a1 + a15 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +15]] = s[old_path][off_s + bit_flips[bits_num * old_path +15]] ? 0 : b;
			break;
		case 109:// e + a2 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 110:// e + a2 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 111:// e + a2 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 112:// e + a2 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 113:// e + a2 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 114:// e + a2 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 115:// e + a2 + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 116:// e + a2 + a15 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +15]] = s[old_path][off_s + bit_flips[bits_num * old_path +15]] ? 0 : b;
			break;
		case 117:// e + a3 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 118:// e + a3 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 119:// e + a3 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 120:// e + a3 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 121:// e + a3 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 122:// e + a3 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 123:// e + a3 + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 124:// e + a3 + a15 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +15]] = s[old_path][off_s + bit_flips[bits_num * old_path +15]] ? 0 : b;
			break;
		case 125:// e + a4 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 126:// e + a4 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 127:// e + a4 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 128:// e + a4 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 129:// e + a4 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 130:// e + a4 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 131:// e + a4 + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 132:// e + a5 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 133:// e + a5 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 134:// e + a5 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 135:// e + a5 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 136:// e + a5 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 137:// e + a5 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 138:// e + a6 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 139:// e + a6 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 140:// e + a6 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 141:// e + a6 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 142:// e + a6 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 143:// e + a7 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 144:// e + a7 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 145:// e + a7 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 146:// e + a7 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 147:// e + a8 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 148:// e + a8 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 149:// e + a8 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 150:// e + a9 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 151:// e + a0 + a1 + a2 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 152:// e + a0 + a1 + a2 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 153:// e + a0 + a1 + a2 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 154:// e + a0 + a1 + a2 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 155:// e + a0 + a1 + a2 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 156:// e + a0 + a1 + a2 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 157:// e + a0 + a1 + a2 + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 158:// e + a0 + a1 + a2 + a15 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +15]] = s[old_path][off_s + bit_flips[bits_num * old_path +15]] ? 0 : b;
			break;
		case 159:// e + a0 + a1 + a3 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 160:// e + a0 + a1 + a3 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 161:// e + a0 + a1 + a3 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 162:// e + a0 + a1 + a3 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 163:// e + a0 + a1 + a3 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 164:// e + a0 + a1 + a4 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 165:// e + a0 + a1 + a4 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 166:// e + a0 + a1 + a4 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 167:// e + a0 + a1 + a5 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 168:// e + a0 + a1 + a5 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 169:// e + a0 + a1 + a6 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 170:// e + a0 + a1 + a7 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 171:// e + a0 + a2 + a3 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 172:// e + a0 + a2 + a3 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 173:// e + a0 + a2 + a3 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 174:// e + a0 + a2 + a3 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 175:// e + a0 + a2 + a4 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 176:// e + a0 + a2 + a5 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 177:// e + a0 + a3 + a4 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 178:// e + a1 + a2 + a3 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 179:// e + a1 + a2 + a3 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 180:// e + a1 + a2 + a3 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 181:// e + a1 + a2 + a4 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 182:// e + a0 + a16 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +16]] = s[old_path][off_s + bit_flips[bits_num * old_path +16]] ? 0 : b;
			break;
		case 183:// e + a0 + a17 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +17]] = s[old_path][off_s + bit_flips[bits_num * old_path +17]] ? 0 : b;
			break;
		case 184:// e + a0 + a18 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +18]] = s[old_path][off_s + bit_flips[bits_num * old_path +18]] ? 0 : b;
			break;
		case 185:// e + a0 + a19 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +19]] = s[old_path][off_s + bit_flips[bits_num * old_path +19]] ? 0 : b;
			break;
		case 186:// e + a0 + a20 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +20]] = s[old_path][off_s + bit_flips[bits_num * old_path +20]] ? 0 : b;
			break;
		case 187:// e + a0 + a21 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +21]] = s[old_path][off_s + bit_flips[bits_num * old_path +21]] ? 0 : b;
			break;
		case 188:// e + a0 + a22 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +22]] = s[old_path][off_s + bit_flips[bits_num * old_path +22]] ? 0 : b;
			break;
		case 189:// e + a0 + a23 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +23]] = s[old_path][off_s + bit_flips[bits_num * old_path +23]] ? 0 : b;
			break;
		case 190:// e + a0 + a24 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +24]] = s[old_path][off_s + bit_flips[bits_num * old_path +24]] ? 0 : b;
			break;
		case 191:// e + a0 + a25 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +25]] = s[old_path][off_s + bit_flips[bits_num * old_path +25]] ? 0 : b;
			break;
		case 192:// e + a0 + a26 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +26]] = s[old_path][off_s + bit_flips[bits_num * old_path +26]] ? 0 : b;
			break;
		case 193:// e + a0 + a27 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +27]] = s[old_path][off_s + bit_flips[bits_num * old_path +27]] ? 0 : b;
			break;
		case 194:// e + a0 + a28 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +28]] = s[old_path][off_s + bit_flips[bits_num * old_path +28]] ? 0 : b;
			break;
		case 195:// e + a0 + a29 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +29]] = s[old_path][off_s + bit_flips[bits_num * old_path +29]] ? 0 : b;
			break;
		case 196:// e + a0 + a30 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +30]] = s[old_path][off_s + bit_flips[bits_num * old_path +30]] ? 0 : b;
			break;
		case 197:// e + a0 + a31 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +31]] = s[old_path][off_s + bit_flips[bits_num * old_path +31]] ? 0 : b;
			break;
		case 198:// e + a1 + a16 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +16]] = s[old_path][off_s + bit_flips[bits_num * old_path +16]] ? 0 : b;
			break;
		case 199:// e + a1 + a17 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +17]] = s[old_path][off_s + bit_flips[bits_num * old_path +17]] ? 0 : b;
			break;
		case 200:// e + a1 + a18 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +18]] = s[old_path][off_s + bit_flips[bits_num * old_path +18]] ? 0 : b;
			break;
		case 201:// e + a1 + a19 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +19]] = s[old_path][off_s + bit_flips[bits_num * old_path +19]] ? 0 : b;
			break;
		case 202:// e + a1 + a20 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +20]] = s[old_path][off_s + bit_flips[bits_num * old_path +20]] ? 0 : b;
			break;
		case 203:// e + a1 + a21 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +21]] = s[old_path][off_s + bit_flips[bits_num * old_path +21]] ? 0 : b;
			break;
		case 204:// e + a1 + a22 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +22]] = s[old_path][off_s + bit_flips[bits_num * old_path +22]] ? 0 : b;
			break;
		case 205:// e + a1 + a23 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +23]] = s[old_path][off_s + bit_flips[bits_num * old_path +23]] ? 0 : b;
			break;
		case 206:// e + a1 + a24 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +24]] = s[old_path][off_s + bit_flips[bits_num * old_path +24]] ? 0 : b;
			break;
		case 207:// e + a1 + a25 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +25]] = s[old_path][off_s + bit_flips[bits_num * old_path +25]] ? 0 : b;
			break;
		case 208:// e + a1 + a26 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +26]] = s[old_path][off_s + bit_flips[bits_num * old_path +26]] ? 0 : b;
			break;
		case 209:// e + a1 + a27 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +27]] = s[old_path][off_s + bit_flips[bits_num * old_path +27]] ? 0 : b;
			break;
		case 210:// e + a1 + a28 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +28]] = s[old_path][off_s + bit_flips[bits_num * old_path +28]] ? 0 : b;
			break;
		case 211:// e + a1 + a29 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +29]] = s[old_path][off_s + bit_flips[bits_num * old_path +29]] ? 0 : b;
			break;
		case 212:// e + a1 + a30 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +30]] = s[old_path][off_s + bit_flips[bits_num * old_path +30]] ? 0 : b;
			break;
		case 213:// e + a1 + a31 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +31]] = s[old_path][off_s + bit_flips[bits_num * old_path +31]] ? 0 : b;
			break;
		case 214:// e + a2 + a16 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +16]] = s[old_path][off_s + bit_flips[bits_num * old_path +16]] ? 0 : b;
			break;
		case 215:// e + a2 + a17 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +17]] = s[old_path][off_s + bit_flips[bits_num * old_path +17]] ? 0 : b;
			break;
		case 216:// e + a2 + a18 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +18]] = s[old_path][off_s + bit_flips[bits_num * old_path +18]] ? 0 : b;
			break;
		case 217:// e + a2 + a19 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +19]] = s[old_path][off_s + bit_flips[bits_num * old_path +19]] ? 0 : b;
			break;
		case 218:// e + a2 + a20 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +20]] = s[old_path][off_s + bit_flips[bits_num * old_path +20]] ? 0 : b;
			break;
		case 219:// e + a2 + a21 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +21]] = s[old_path][off_s + bit_flips[bits_num * old_path +21]] ? 0 : b;
			break;
		case 220:// e + a2 + a22 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +22]] = s[old_path][off_s + bit_flips[bits_num * old_path +22]] ? 0 : b;
			break;
		case 221:// e + a3 + a16 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +16]] = s[old_path][off_s + bit_flips[bits_num * old_path +16]] ? 0 : b;
			break;
		case 222:// e + a3 + a17 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +17]] = s[old_path][off_s + bit_flips[bits_num * old_path +17]] ? 0 : b;
			break;
		case 223:// e + a0 + a1 + a2 + a16 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +16]] = s[old_path][off_s + bit_flips[bits_num * old_path +16]] ? 0 : b;
			break;
		case 224:// e + a0 + a1 + a2 + a17 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +17]] = s[old_path][off_s + bit_flips[bits_num * old_path +17]] ? 0 : b;
			break;
		case 225:// e + a0 + a32 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +32]] = s[old_path][off_s + bit_flips[bits_num * old_path +32]] ? 0 : b;
			break;
		case 226:// e + a0 + a33 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +33]] = s[old_path][off_s + bit_flips[bits_num * old_path +33]] ? 0 : b;
			break;
		case 227:// e + a0 + a34 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +34]] = s[old_path][off_s + bit_flips[bits_num * old_path +34]] ? 0 : b;
			break;
		case 228:// e + a0 + a35 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +35]] = s[old_path][off_s + bit_flips[bits_num * old_path +35]] ? 0 : b;
			break;
		case 229:// e + a0 + a36 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +36]] = s[old_path][off_s + bit_flips[bits_num * old_path +36]] ? 0 : b;
			break;
		case 230:// e + a0 + a37 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +37]] = s[old_path][off_s + bit_flips[bits_num * old_path +37]] ? 0 : b;
			break;
		case 231:// e + a0 + a38 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +38]] = s[old_path][off_s + bit_flips[bits_num * old_path +38]] ? 0 : b;
			break;
		case 232:// e + a0 + a39 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +39]] = s[old_path][off_s + bit_flips[bits_num * old_path +39]] ? 0 : b;
			break;
		case 233:// e + a0 + a40 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +40]] = s[old_path][off_s + bit_flips[bits_num * old_path +40]] ? 0 : b;
			break;
		case 234:// e + a0 + a41 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +41]] = s[old_path][off_s + bit_flips[bits_num * old_path +41]] ? 0 : b;
			break;
		case 235:// e + a0 + a42 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +42]] = s[old_path][off_s + bit_flips[bits_num * old_path +42]] ? 0 : b;
			break;
		case 236:// e + a0 + a43 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +43]] = s[old_path][off_s + bit_flips[bits_num * old_path +43]] ? 0 : b;
			break;
		case 237:// e + a0 + a44 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +44]] = s[old_path][off_s + bit_flips[bits_num * old_path +44]] ? 0 : b;
			break;
		case 238:// e + a0 + a45 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +45]] = s[old_path][off_s + bit_flips[bits_num * old_path +45]] ? 0 : b;
			break;
		case 239:// e + a0 + a46 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +46]] = s[old_path][off_s + bit_flips[bits_num * old_path +46]] ? 0 : b;
			break;
		case 240:// e + a0 + a47 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +47]] = s[old_path][off_s + bit_flips[bits_num * old_path +47]] ? 0 : b;
			break;
		case 241:// e + a0 + a48 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +48]] = s[old_path][off_s + bit_flips[bits_num * old_path +48]] ? 0 : b;
			break;
		case 242:// e + a0 + a49 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +49]] = s[old_path][off_s + bit_flips[bits_num * old_path +49]] ? 0 : b;
			break;
		case 243:// e + a0 + a50 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +50]] = s[old_path][off_s + bit_flips[bits_num * old_path +50]] ? 0 : b;
			break;
		case 244:// e + a0 + a51 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +51]] = s[old_path][off_s + bit_flips[bits_num * old_path +51]] ? 0 : b;
			break;
		case 245:// e + a0 + a52 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +52]] = s[old_path][off_s + bit_flips[bits_num * old_path +52]] ? 0 : b;
			break;
		case 246:// e + a0 + a53 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +53]] = s[old_path][off_s + bit_flips[bits_num * old_path +53]] ? 0 : b;
			break;
		case 247:// e + a0 + a54 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +54]] = s[old_path][off_s + bit_flips[bits_num * old_path +54]] ? 0 : b;
			break;
		case 248:// e + a0 + a55 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +55]] = s[old_path][off_s + bit_flips[bits_num * old_path +55]] ? 0 : b;
			break;
		case 249:// e + a0 + a56 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +56]] = s[old_path][off_s + bit_flips[bits_num * old_path +56]] ? 0 : b;
			break;
		case 250:// e + a0 + a57 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +57]] = s[old_path][off_s + bit_flips[bits_num * old_path +57]] ? 0 : b;
			break;
		case 251:// e + a0 + a58 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +58]] = s[old_path][off_s + bit_flips[bits_num * old_path +58]] ? 0 : b;
			break;
		case 252:// e + a0 + a59 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +59]] = s[old_path][off_s + bit_flips[bits_num * old_path +59]] ? 0 : b;
			break;
		case 253:// e + a0 + a60 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +60]] = s[old_path][off_s + bit_flips[bits_num * old_path +60]] ? 0 : b;
			break;
		case 254:// e + a0 + a61 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +61]] = s[old_path][off_s + bit_flips[bits_num * old_path +61]] ? 0 : b;
			break;
		case 255:// e + a0 + a62 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +62]] = s[old_path][off_s + bit_flips[bits_num * old_path +62]] ? 0 : b;
			break;
		case 256:// e + a0 + a63 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +63]] = s[old_path][off_s + bit_flips[bits_num * old_path +63]] ? 0 : b;
			break;
		case 257:// e + a1 + a32 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +32]] = s[old_path][off_s + bit_flips[bits_num * old_path +32]] ? 0 : b;
			break;
		default:
			throw tools::runtime_error(__FILE__, __LINE__, __func__, "Flip bits error on spc node.");
			break;
		}
		break;
	default:
		throw tools::runtime_error(__FILE__, __LINE__, __func__, "Flip bits error on SPC node.");
		break;
	}
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::delete_path(int path_id)
{
	const auto old_path = paths[path_id];
	for (auto i = 0; i < m; i++)
		n_array_ref[path_2_array[old_path][i]][i]--;

	paths[path_id] = paths[--n_active_paths];
	paths[n_active_paths] = old_path;
}

template <typename B, typename R, class API_polar>
int Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::select_best_path(const size_t frame_id)
{
	best_path = -1;
	for (auto i = 0; i < n_active_paths; i++)
		if (best_path == -1 || metrics[paths[i]] < metrics[best_path])
			best_path = paths[i];

	if (best_path == -1)
		best_path = 0;

	return n_active_paths;
}

template <typename B, typename R, class API_polar>
int Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::up_ref_array_idx(const int path, const int r_d)
{
	auto old_array = path_2_array[path][r_d];

	// if more than 1 path points to the array
	if (n_array_ref[old_array][r_d] > 1)
	{
		// allocate new array to given path, r_d
		n_array_ref[old_array][r_d]--;

		auto new_array = 0;
		while (n_array_ref[new_array][r_d]) new_array++;

		path_2_array[path     ][r_d] = new_array;
		n_array_ref [new_array][r_d]++;

		return new_array;
	}

	return old_array;
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::erase_bad_paths()
{
	// erase bad paths
	auto k = 0;
	auto n_active_paths_cpy = n_active_paths;
	for (auto i = 0; i < n_active_paths_cpy; i++)
		if (dup_count[paths[k]] == 0)
			delete_path(k);
		else
			k++;
}

template <typename B, typename R, class API_polar>
int Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::duplicate_tree(const int old_path, const int off_l, const int off_s, const int n_elmts)
{
	const auto new_path = paths[n_active_paths++];

	std::copy(path_2_array[old_path].begin(), path_2_array[old_path].end(), path_2_array[new_path].begin());

	for (auto i = 0; i < m; i++)
		n_array_ref[path_2_array[new_path][i]][i]++;

	std::copy(s[old_path].begin(), s[old_path].begin() + off_s + n_elmts, s[new_path].begin());

	return new_path;
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::_partial_sort(const R* values, std::vector<int> &pos, int n_elmts, int k)
{
    struct Node
    {
        float val;
        int idx;
        Node(float a=0, int b=0):
            val(a), idx(b) {}
    };

    struct cmp
    {
        bool operator()(Node a, Node b)
        {
            return a.val < b.val;
        }
    };

    std::priority_queue<Node, std::vector<Node>, cmp> p;

    for (int i = 0; i < k; ++i)
    {
        p.push(Node(values[i], i));
    } 

    for (int i = k; i < n_elmts; ++i)
    {
        if (p.top().val > values[i])
        {
            p.pop();
            p.push(Node(values[i], i));
        }
    }

    pos.resize(n_elmts);

    for (int i = k-1; i >= 0; --i)
    {
        pos[i] = p.top().idx;
        p.pop();
    }
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::flip_bits_r1(const int old_path, const int new_path, const int dup, const int off_s, const int n_elmts, const int bits_num)
{
	constexpr B b = tools::bit_init<B>();

	switch (L)
	{
	case 2:
		switch( dup )
		{
		case 0:
			// nothing to do
			break;
		case 1:// e + a0 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			break;
		default:
			throw tools::runtime_error(__FILE__, __LINE__, __func__, "Flip bits error on rate 1 node.");
			break;
		}
		break;
	case 4:
		switch( dup )
		{
		case 0:
			// nothing to do
			break;
		case 1:// e + a0 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			break;
		case 2:// e + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		case 3:// e + a0 + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		case 4:// e + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		default:
			throw tools::runtime_error(__FILE__, __LINE__, __func__, "Flip bits error on rate 1 node.");
			break;
		}
		break;
	case 8:
		switch( dup )
		{
		case 0:
			// nothing to do
			break;
		case 1:// e + a0 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			break;
		case 2:// e + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		case 3:// e + a0 + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		case 4:// e + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 5:// e + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 6:// e + a0 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 7:// e + a0 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 8:// e + a1 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 9:// e + a0 + a1 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 10:// e + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 11:// e + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 12:// e + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		default:
			throw tools::runtime_error(__FILE__, __LINE__, __func__, "Flip bits error on rate 1 node.");
			break;
		}
		break;
	case 16:
		switch( dup )
		{
		case 0:
			// nothing to do
			break;
		case 1:// e + a0 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			break;
		case 2:// e + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		case 3:// e + a0 + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		case 4:// e + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 5:// e + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 6:// e + a0 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 7:// e + a0 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 8:// e + a1 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 9:// e + a1 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 10:// e + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 11:// e + a0 + a1 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 12:// e + a0 + a1 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 13:// e + a0 + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 14:// e + a1 + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 15:// e + a0 + a1 + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 16:// e + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 17:// e + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 18:// e + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 19:// e + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 20:// e + a0 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 21:// e + a0 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 22:// e + a0 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 23:// e + a0 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 24:// e + a1 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 25:// e + a1 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 26:// e + a2 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 27:// e + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 28:// e + a0 + a1 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 29:// e + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 30:// e + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 31:// e + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 32:// e + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 33:// e + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 34:// e + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 35:// e + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		default:
			throw tools::runtime_error(__FILE__, __LINE__, __func__, "Flip bits error on rate 1 node.");
			break;
		}
		break;
	case 32:
		switch( dup )
		{
		case 0:
			// nothing to do
			break;
		case 1:// e + a0 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			break;
		case 2:// e + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		case 3:// e + a0 + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		case 4:// e + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 5:// e + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 6:// e + a0 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 7:// e + a0 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 8:// e + a1 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 9:// e + a1 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 10:// e + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 11:// e + a0 + a1 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 12:// e + a0 + a1 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 13:// e + a0 + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 14:// e + a1 + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 15:// e + a0 + a1 + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 16:// e + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 17:// e + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 18:// e + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 19:// e + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 20:// e + a0 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 21:// e + a0 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 22:// e + a0 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 23:// e + a0 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 24:// e + a1 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 25:// e + a1 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 26:// e + a1 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 27:// e + a1 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 28:// e + a2 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 29:// e + a2 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 30:// e + a2 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 31:// e + a2 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 32:// e + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 33:// e + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 34:// e + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 35:// e + a3 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 36:// e + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 37:// e + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 38:// e + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 39:// e + a0 + a1 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 40:// e + a0 + a1 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 41:// e + a0 + a1 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 42:// e + a0 + a1 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 43:// e + a0 + a2 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 44:// e + a0 + a2 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 45:// e + a0 + a2 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 46:// e + a0 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 47:// e + a0 + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 48:// e + a0 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 49:// e + a1 + a2 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 50:// e + a1 + a2 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 51:// e + a1 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 52:// e + a2 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 53:// e + a0 + a1 + a2 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 54:// e + a0 + a1 + a2 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 55:// e + a0 + a1 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 56:// e + a0 + a2 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 57:// e + a1 + a2 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 58:// e + a0 + a1 + a2 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 59:// e + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 60:// e + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 61:// e + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 62:// e + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 63:// e + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 64:// e + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 65:// e + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 66:// e + a15 
			s[new_path][off_s + bit_flips[bits_num * old_path +15]] = s[old_path][off_s + bit_flips[bits_num * old_path +15]] ? 0 : b;
			break;
		case 67:// e + a0 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 68:// e + a0 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 69:// e + a0 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 70:// e + a0 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 71:// e + a0 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 72:// e + a0 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 73:// e + a0 + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 74:// e + a0 + a15 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +15]] = s[old_path][off_s + bit_flips[bits_num * old_path +15]] ? 0 : b;
			break;
		case 75:// e + a1 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 76:// e + a1 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 77:// e + a1 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 78:// e + a2 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 79:// e + a0 + a1 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 80:// e + a16 
			s[new_path][off_s + bit_flips[bits_num * old_path +16]] = s[old_path][off_s + bit_flips[bits_num * old_path +16]] ? 0 : b;
			break;
		case 81:// e + a17 
			s[new_path][off_s + bit_flips[bits_num * old_path +17]] = s[old_path][off_s + bit_flips[bits_num * old_path +17]] ? 0 : b;
			break;
		case 82:// e + a18 
			s[new_path][off_s + bit_flips[bits_num * old_path +18]] = s[old_path][off_s + bit_flips[bits_num * old_path +18]] ? 0 : b;
			break;
		case 83:// e + a19 
			s[new_path][off_s + bit_flips[bits_num * old_path +19]] = s[old_path][off_s + bit_flips[bits_num * old_path +19]] ? 0 : b;
			break;
		case 84:// e + a20 
			s[new_path][off_s + bit_flips[bits_num * old_path +20]] = s[old_path][off_s + bit_flips[bits_num * old_path +20]] ? 0 : b;
			break;
		case 85:// e + a21 
			s[new_path][off_s + bit_flips[bits_num * old_path +21]] = s[old_path][off_s + bit_flips[bits_num * old_path +21]] ? 0 : b;
			break;
		case 86:// e + a22 
			s[new_path][off_s + bit_flips[bits_num * old_path +22]] = s[old_path][off_s + bit_flips[bits_num * old_path +22]] ? 0 : b;
			break;
		case 87:// e + a23 
			s[new_path][off_s + bit_flips[bits_num * old_path +23]] = s[old_path][off_s + bit_flips[bits_num * old_path +23]] ? 0 : b;
			break;
		case 88:// e + a24 
			s[new_path][off_s + bit_flips[bits_num * old_path +24]] = s[old_path][off_s + bit_flips[bits_num * old_path +24]] ? 0 : b;
			break;
		case 89:// e + a25 
			s[new_path][off_s + bit_flips[bits_num * old_path +25]] = s[old_path][off_s + bit_flips[bits_num * old_path +25]] ? 0 : b;
			break;
		case 90:// e + a26 
			s[new_path][off_s + bit_flips[bits_num * old_path +26]] = s[old_path][off_s + bit_flips[bits_num * old_path +26]] ? 0 : b;
			break;
		case 91:// e + a27 
			s[new_path][off_s + bit_flips[bits_num * old_path +27]] = s[old_path][off_s + bit_flips[bits_num * old_path +27]] ? 0 : b;
			break;
		case 92:// e + a28 
			s[new_path][off_s + bit_flips[bits_num * old_path +28]] = s[old_path][off_s + bit_flips[bits_num * old_path +28]] ? 0 : b;
			break;
		case 93:// e + a29 
			s[new_path][off_s + bit_flips[bits_num * old_path +29]] = s[old_path][off_s + bit_flips[bits_num * old_path +29]] ? 0 : b;
			break;
		case 94:// e + a30 
			s[new_path][off_s + bit_flips[bits_num * old_path +30]] = s[old_path][off_s + bit_flips[bits_num * old_path +30]] ? 0 : b;
			break;
		default:
			throw tools::runtime_error(__FILE__, __LINE__, __func__, "Flip bits error on rate 1 node.");
			break;
		}
		break;
	case 64:
		switch( dup )
		{
		case 0:
			// nothing to do
			break;
		case 1:// e + a0 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			break;
		case 2:// e + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		case 3:// e + a0 + a1 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			break;
		case 4:// e + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 5:// e + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 6:// e + a0 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 7:// e + a0 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 8:// e + a1 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 9:// e + a1 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 10:// e + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 11:// e + a0 + a1 + a2 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			break;
		case 12:// e + a0 + a1 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 13:// e + a0 + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 14:// e + a1 + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 15:// e + a0 + a1 + a2 + a3 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			break;
		case 16:// e + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 17:// e + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 18:// e + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 19:// e + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 20:// e + a0 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 21:// e + a0 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 22:// e + a0 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 23:// e + a0 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 24:// e + a1 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 25:// e + a1 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 26:// e + a1 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 27:// e + a1 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 28:// e + a2 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 29:// e + a2 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 30:// e + a2 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 31:// e + a2 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 32:// e + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 33:// e + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 34:// e + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 35:// e + a3 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 36:// e + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 37:// e + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 38:// e + a4 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 39:// e + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 40:// e + a5 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 41:// e + a6 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 42:// e + a0 + a1 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 43:// e + a0 + a1 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 44:// e + a0 + a1 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 45:// e + a0 + a1 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 46:// e + a0 + a2 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 47:// e + a0 + a2 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 48:// e + a0 + a2 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 49:// e + a0 + a2 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 50:// e + a0 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 51:// e + a0 + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 52:// e + a0 + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 53:// e + a0 + a3 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 54:// e + a0 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 55:// e + a0 + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 56:// e + a0 + a4 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 57:// e + a0 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 58:// e + a0 + a5 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 59:// e + a0 + a6 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 60:// e + a1 + a2 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 61:// e + a1 + a2 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 62:// e + a1 + a2 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 63:// e + a1 + a2 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 64:// e + a1 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 65:// e + a1 + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 66:// e + a1 + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 67:// e + a1 + a3 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 68:// e + a1 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 69:// e + a1 + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 70:// e + a1 + a4 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 71:// e + a1 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 72:// e + a2 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 73:// e + a2 + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 74:// e + a2 + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 75:// e + a2 + a3 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 76:// e + a2 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 77:// e + a2 + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 78:// e + a2 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 79:// e + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 80:// e + a3 + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 81:// e + a3 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 82:// e + a4 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 83:// e + a0 + a1 + a2 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 84:// e + a0 + a1 + a2 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 85:// e + a0 + a1 + a2 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 86:// e + a0 + a1 + a2 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 87:// e + a0 + a1 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 88:// e + a0 + a1 + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 89:// e + a0 + a1 + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 90:// e + a0 + a1 + a3 + a7 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			break;
		case 91:// e + a0 + a1 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 92:// e + a0 + a1 + a4 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 93:// e + a0 + a1 + a5 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 94:// e + a0 + a2 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 95:// e + a0 + a2 + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 96:// e + a0 + a2 + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 97:// e + a0 + a2 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 98:// e + a0 + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 99:// e + a1 + a2 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 100:// e + a1 + a2 + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 101:// e + a1 + a2 + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 102:// e + a1 + a2 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 103:// e + a1 + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 104:// e + a2 + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 105:// e + a0 + a1 + a2 + a3 + a4 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			break;
		case 106:// e + a0 + a1 + a2 + a3 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 107:// e + a0 + a1 + a2 + a3 + a6 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			break;
		case 108:// e + a0 + a1 + a2 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 109:// e + a0 + a1 + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 110:// e + a0 + a2 + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 111:// e + a1 + a2 + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		case 112:// e + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 113:// e + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 114:// e + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 115:// e + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 116:// e + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 117:// e + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 118:// e + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 119:// e + a15 
			s[new_path][off_s + bit_flips[bits_num * old_path +15]] = s[old_path][off_s + bit_flips[bits_num * old_path +15]] ? 0 : b;
			break;
		case 120:// e + a0 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 121:// e + a0 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 122:// e + a0 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 123:// e + a0 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 124:// e + a0 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 125:// e + a0 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 126:// e + a0 + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 127:// e + a0 + a15 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +15]] = s[old_path][off_s + bit_flips[bits_num * old_path +15]] ? 0 : b;
			break;
		case 128:// e + a1 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 129:// e + a1 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 130:// e + a1 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 131:// e + a1 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 132:// e + a1 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 133:// e + a1 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 134:// e + a1 + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 135:// e + a1 + a15 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +15]] = s[old_path][off_s + bit_flips[bits_num * old_path +15]] ? 0 : b;
			break;
		case 136:// e + a2 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 137:// e + a2 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 138:// e + a2 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 139:// e + a2 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 140:// e + a2 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 141:// e + a2 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 142:// e + a2 + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 143:// e + a2 + a15 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +15]] = s[old_path][off_s + bit_flips[bits_num * old_path +15]] ? 0 : b;
			break;
		case 144:// e + a3 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 145:// e + a3 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 146:// e + a3 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 147:// e + a3 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 148:// e + a3 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 149:// e + a3 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 150:// e + a4 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 151:// e + a4 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 152:// e + a4 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 153:// e + a4 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 154:// e + a4 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 155:// e + a5 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 156:// e + a5 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 157:// e + a5 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 158:// e + a5 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 159:// e + a6 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 160:// e + a6 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 161:// e + a6 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +6]] = s[old_path][off_s + bit_flips[bits_num * old_path +6]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 162:// e + a7 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 163:// e + a7 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 164:// e + a7 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +7]] = s[old_path][off_s + bit_flips[bits_num * old_path +7]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 165:// e + a8 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 166:// e + a0 + a1 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 167:// e + a0 + a1 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 168:// e + a0 + a1 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 169:// e + a0 + a1 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 170:// e + a0 + a1 + a12 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +12]] = s[old_path][off_s + bit_flips[bits_num * old_path +12]] ? 0 : b;
			break;
		case 171:// e + a0 + a1 + a13 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +13]] = s[old_path][off_s + bit_flips[bits_num * old_path +13]] ? 0 : b;
			break;
		case 172:// e + a0 + a1 + a14 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +14]] = s[old_path][off_s + bit_flips[bits_num * old_path +14]] ? 0 : b;
			break;
		case 173:// e + a0 + a1 + a15 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +15]] = s[old_path][off_s + bit_flips[bits_num * old_path +15]] ? 0 : b;
			break;
		case 174:// e + a0 + a2 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 175:// e + a0 + a2 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 176:// e + a0 + a2 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 177:// e + a0 + a2 + a11 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +11]] = s[old_path][off_s + bit_flips[bits_num * old_path +11]] ? 0 : b;
			break;
		case 178:// e + a0 + a3 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 179:// e + a0 + a3 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 180:// e + a0 + a4 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 181:// e + a1 + a2 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 182:// e + a1 + a2 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 183:// e + a1 + a2 + a10 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +10]] = s[old_path][off_s + bit_flips[bits_num * old_path +10]] ? 0 : b;
			break;
		case 184:// e + a0 + a1 + a2 + a8 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +8]] = s[old_path][off_s + bit_flips[bits_num * old_path +8]] ? 0 : b;
			break;
		case 185:// e + a0 + a1 + a2 + a9 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +9]] = s[old_path][off_s + bit_flips[bits_num * old_path +9]] ? 0 : b;
			break;
		case 186:// e + a16 
			s[new_path][off_s + bit_flips[bits_num * old_path +16]] = s[old_path][off_s + bit_flips[bits_num * old_path +16]] ? 0 : b;
			break;
		case 187:// e + a17 
			s[new_path][off_s + bit_flips[bits_num * old_path +17]] = s[old_path][off_s + bit_flips[bits_num * old_path +17]] ? 0 : b;
			break;
		case 188:// e + a18 
			s[new_path][off_s + bit_flips[bits_num * old_path +18]] = s[old_path][off_s + bit_flips[bits_num * old_path +18]] ? 0 : b;
			break;
		case 189:// e + a19 
			s[new_path][off_s + bit_flips[bits_num * old_path +19]] = s[old_path][off_s + bit_flips[bits_num * old_path +19]] ? 0 : b;
			break;
		case 190:// e + a20 
			s[new_path][off_s + bit_flips[bits_num * old_path +20]] = s[old_path][off_s + bit_flips[bits_num * old_path +20]] ? 0 : b;
			break;
		case 191:// e + a21 
			s[new_path][off_s + bit_flips[bits_num * old_path +21]] = s[old_path][off_s + bit_flips[bits_num * old_path +21]] ? 0 : b;
			break;
		case 192:// e + a22 
			s[new_path][off_s + bit_flips[bits_num * old_path +22]] = s[old_path][off_s + bit_flips[bits_num * old_path +22]] ? 0 : b;
			break;
		case 193:// e + a23 
			s[new_path][off_s + bit_flips[bits_num * old_path +23]] = s[old_path][off_s + bit_flips[bits_num * old_path +23]] ? 0 : b;
			break;
		case 194:// e + a24 
			s[new_path][off_s + bit_flips[bits_num * old_path +24]] = s[old_path][off_s + bit_flips[bits_num * old_path +24]] ? 0 : b;
			break;
		case 195:// e + a25 
			s[new_path][off_s + bit_flips[bits_num * old_path +25]] = s[old_path][off_s + bit_flips[bits_num * old_path +25]] ? 0 : b;
			break;
		case 196:// e + a26 
			s[new_path][off_s + bit_flips[bits_num * old_path +26]] = s[old_path][off_s + bit_flips[bits_num * old_path +26]] ? 0 : b;
			break;
		case 197:// e + a27 
			s[new_path][off_s + bit_flips[bits_num * old_path +27]] = s[old_path][off_s + bit_flips[bits_num * old_path +27]] ? 0 : b;
			break;
		case 198:// e + a28 
			s[new_path][off_s + bit_flips[bits_num * old_path +28]] = s[old_path][off_s + bit_flips[bits_num * old_path +28]] ? 0 : b;
			break;
		case 199:// e + a29 
			s[new_path][off_s + bit_flips[bits_num * old_path +29]] = s[old_path][off_s + bit_flips[bits_num * old_path +29]] ? 0 : b;
			break;
		case 200:// e + a30 
			s[new_path][off_s + bit_flips[bits_num * old_path +30]] = s[old_path][off_s + bit_flips[bits_num * old_path +30]] ? 0 : b;
			break;
		case 201:// e + a31 
			s[new_path][off_s + bit_flips[bits_num * old_path +31]] = s[old_path][off_s + bit_flips[bits_num * old_path +31]] ? 0 : b;
			break;
		case 202:// e + a0 + a16 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +16]] = s[old_path][off_s + bit_flips[bits_num * old_path +16]] ? 0 : b;
			break;
		case 203:// e + a0 + a17 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +17]] = s[old_path][off_s + bit_flips[bits_num * old_path +17]] ? 0 : b;
			break;
		case 204:// e + a0 + a18 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +18]] = s[old_path][off_s + bit_flips[bits_num * old_path +18]] ? 0 : b;
			break;
		case 205:// e + a0 + a19 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +19]] = s[old_path][off_s + bit_flips[bits_num * old_path +19]] ? 0 : b;
			break;
		case 206:// e + a0 + a20 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +20]] = s[old_path][off_s + bit_flips[bits_num * old_path +20]] ? 0 : b;
			break;
		case 207:// e + a0 + a21 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +21]] = s[old_path][off_s + bit_flips[bits_num * old_path +21]] ? 0 : b;
			break;
		case 208:// e + a0 + a22 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +22]] = s[old_path][off_s + bit_flips[bits_num * old_path +22]] ? 0 : b;
			break;
		case 209:// e + a0 + a23 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +23]] = s[old_path][off_s + bit_flips[bits_num * old_path +23]] ? 0 : b;
			break;
		case 210:// e + a0 + a24 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +24]] = s[old_path][off_s + bit_flips[bits_num * old_path +24]] ? 0 : b;
			break;
		case 211:// e + a0 + a25 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +25]] = s[old_path][off_s + bit_flips[bits_num * old_path +25]] ? 0 : b;
			break;
		case 212:// e + a0 + a26 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +26]] = s[old_path][off_s + bit_flips[bits_num * old_path +26]] ? 0 : b;
			break;
		case 213:// e + a0 + a27 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +27]] = s[old_path][off_s + bit_flips[bits_num * old_path +27]] ? 0 : b;
			break;
		case 214:// e + a0 + a28 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +28]] = s[old_path][off_s + bit_flips[bits_num * old_path +28]] ? 0 : b;
			break;
		case 215:// e + a0 + a29 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +29]] = s[old_path][off_s + bit_flips[bits_num * old_path +29]] ? 0 : b;
			break;
		case 216:// e + a0 + a30 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +30]] = s[old_path][off_s + bit_flips[bits_num * old_path +30]] ? 0 : b;
			break;
		case 217:// e + a0 + a31 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +31]] = s[old_path][off_s + bit_flips[bits_num * old_path +31]] ? 0 : b;
			break;
		case 218:// e + a1 + a16 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +16]] = s[old_path][off_s + bit_flips[bits_num * old_path +16]] ? 0 : b;
			break;
		case 219:// e + a1 + a17 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +17]] = s[old_path][off_s + bit_flips[bits_num * old_path +17]] ? 0 : b;
			break;
		case 220:// e + a1 + a18 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +18]] = s[old_path][off_s + bit_flips[bits_num * old_path +18]] ? 0 : b;
			break;
		case 221:// e + a1 + a19 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +19]] = s[old_path][off_s + bit_flips[bits_num * old_path +19]] ? 0 : b;
			break;
		case 222:// e + a1 + a20 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +20]] = s[old_path][off_s + bit_flips[bits_num * old_path +20]] ? 0 : b;
			break;
		case 223:// e + a1 + a21 
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +21]] = s[old_path][off_s + bit_flips[bits_num * old_path +21]] ? 0 : b;
			break;
		case 224:// e + a2 + a16 
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +16]] = s[old_path][off_s + bit_flips[bits_num * old_path +16]] ? 0 : b;
			break;
		case 225:// e + a0 + a1 + a16 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +16]] = s[old_path][off_s + bit_flips[bits_num * old_path +16]] ? 0 : b;
			break;
		case 226:// e + a32 
			s[new_path][off_s + bit_flips[bits_num * old_path +32]] = s[old_path][off_s + bit_flips[bits_num * old_path +32]] ? 0 : b;
			break;
		case 227:// e + a33 
			s[new_path][off_s + bit_flips[bits_num * old_path +33]] = s[old_path][off_s + bit_flips[bits_num * old_path +33]] ? 0 : b;
			break;
		case 228:// e + a34 
			s[new_path][off_s + bit_flips[bits_num * old_path +34]] = s[old_path][off_s + bit_flips[bits_num * old_path +34]] ? 0 : b;
			break;
		case 229:// e + a35 
			s[new_path][off_s + bit_flips[bits_num * old_path +35]] = s[old_path][off_s + bit_flips[bits_num * old_path +35]] ? 0 : b;
			break;
		case 230:// e + a36 
			s[new_path][off_s + bit_flips[bits_num * old_path +36]] = s[old_path][off_s + bit_flips[bits_num * old_path +36]] ? 0 : b;
			break;
		case 231:// e + a37 
			s[new_path][off_s + bit_flips[bits_num * old_path +37]] = s[old_path][off_s + bit_flips[bits_num * old_path +37]] ? 0 : b;
			break;
		case 232:// e + a38 
			s[new_path][off_s + bit_flips[bits_num * old_path +38]] = s[old_path][off_s + bit_flips[bits_num * old_path +38]] ? 0 : b;
			break;
		case 233:// e + a39 
			s[new_path][off_s + bit_flips[bits_num * old_path +39]] = s[old_path][off_s + bit_flips[bits_num * old_path +39]] ? 0 : b;
			break;
		case 234:// e + a40 
			s[new_path][off_s + bit_flips[bits_num * old_path +40]] = s[old_path][off_s + bit_flips[bits_num * old_path +40]] ? 0 : b;
			break;
		case 235:// e + a41 
			s[new_path][off_s + bit_flips[bits_num * old_path +41]] = s[old_path][off_s + bit_flips[bits_num * old_path +41]] ? 0 : b;
			break;
		case 236:// e + a42 
			s[new_path][off_s + bit_flips[bits_num * old_path +42]] = s[old_path][off_s + bit_flips[bits_num * old_path +42]] ? 0 : b;
			break;
		case 237:// e + a43 
			s[new_path][off_s + bit_flips[bits_num * old_path +43]] = s[old_path][off_s + bit_flips[bits_num * old_path +43]] ? 0 : b;
			break;
		case 238:// e + a44 
			s[new_path][off_s + bit_flips[bits_num * old_path +44]] = s[old_path][off_s + bit_flips[bits_num * old_path +44]] ? 0 : b;
			break;
		case 239:// e + a45 
			s[new_path][off_s + bit_flips[bits_num * old_path +45]] = s[old_path][off_s + bit_flips[bits_num * old_path +45]] ? 0 : b;
			break;
		case 240:// e + a46 
			s[new_path][off_s + bit_flips[bits_num * old_path +46]] = s[old_path][off_s + bit_flips[bits_num * old_path +46]] ? 0 : b;
			break;
		case 241:// e + a47 
			s[new_path][off_s + bit_flips[bits_num * old_path +47]] = s[old_path][off_s + bit_flips[bits_num * old_path +47]] ? 0 : b;
			break;
		case 242:// e + a48 
			s[new_path][off_s + bit_flips[bits_num * old_path +48]] = s[old_path][off_s + bit_flips[bits_num * old_path +48]] ? 0 : b;
			break;
		case 243:// e + a49 
			s[new_path][off_s + bit_flips[bits_num * old_path +49]] = s[old_path][off_s + bit_flips[bits_num * old_path +49]] ? 0 : b;
			break;
		case 244:// e + a50 
			s[new_path][off_s + bit_flips[bits_num * old_path +50]] = s[old_path][off_s + bit_flips[bits_num * old_path +50]] ? 0 : b;
			break;
		case 245:// e + a51 
			s[new_path][off_s + bit_flips[bits_num * old_path +51]] = s[old_path][off_s + bit_flips[bits_num * old_path +51]] ? 0 : b;
			break;
		case 246:// e + a52 
			s[new_path][off_s + bit_flips[bits_num * old_path +52]] = s[old_path][off_s + bit_flips[bits_num * old_path +52]] ? 0 : b;
			break;
		case 247:// e + a53 
			s[new_path][off_s + bit_flips[bits_num * old_path +53]] = s[old_path][off_s + bit_flips[bits_num * old_path +53]] ? 0 : b;
			break;
		case 248:// e + a54 
			s[new_path][off_s + bit_flips[bits_num * old_path +54]] = s[old_path][off_s + bit_flips[bits_num * old_path +54]] ? 0 : b;
			break;
		case 249:// e + a55 
			s[new_path][off_s + bit_flips[bits_num * old_path +55]] = s[old_path][off_s + bit_flips[bits_num * old_path +55]] ? 0 : b;
			break;
		case 250:// e + a56 
			s[new_path][off_s + bit_flips[bits_num * old_path +56]] = s[old_path][off_s + bit_flips[bits_num * old_path +56]] ? 0 : b;
			break;
		case 251:// e + a57 
			s[new_path][off_s + bit_flips[bits_num * old_path +57]] = s[old_path][off_s + bit_flips[bits_num * old_path +57]] ? 0 : b;
			break;
		case 252:// e + a58 
			s[new_path][off_s + bit_flips[bits_num * old_path +58]] = s[old_path][off_s + bit_flips[bits_num * old_path +58]] ? 0 : b;
			break;
		case 253:// e + a59 
			s[new_path][off_s + bit_flips[bits_num * old_path +59]] = s[old_path][off_s + bit_flips[bits_num * old_path +59]] ? 0 : b;
			break;
		case 254:// e + a60 
			s[new_path][off_s + bit_flips[bits_num * old_path +60]] = s[old_path][off_s + bit_flips[bits_num * old_path +60]] ? 0 : b;
			break;
		case 255:// e + a61 
			s[new_path][off_s + bit_flips[bits_num * old_path +61]] = s[old_path][off_s + bit_flips[bits_num * old_path +61]] ? 0 : b;
			break;
		case 256:// e + a62 
			s[new_path][off_s + bit_flips[bits_num * old_path +62]] = s[old_path][off_s + bit_flips[bits_num * old_path +62]] ? 0 : b;
			break;
		case 257:// e + a0 + a1 + a2 + a3 + a4 + a5 
			s[new_path][off_s + bit_flips[bits_num * old_path +0]] = s[old_path][off_s + bit_flips[bits_num * old_path +0]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +1]] = s[old_path][off_s + bit_flips[bits_num * old_path +1]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +2]] = s[old_path][off_s + bit_flips[bits_num * old_path +2]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +3]] = s[old_path][off_s + bit_flips[bits_num * old_path +3]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +4]] = s[old_path][off_s + bit_flips[bits_num * old_path +4]] ? 0 : b;
			s[new_path][off_s + bit_flips[bits_num * old_path +5]] = s[old_path][off_s + bit_flips[bits_num * old_path +5]] ? 0 : b;
			break;
		default:
			throw tools::runtime_error(__FILE__, __LINE__, __func__, "Flip bits error on rate 1 node.");
			break;
		}
		break;
	default:
		break;
	}
}

template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::_update_r1_vec(const int path, const int array, const int n_elmts, const int bits_num, const int off_l)
{
	int c_num = r1_mc_size[(int)std::log2((n_elmts > L) ? L : n_elmts) - 1];
	
	std ::vector<R> pen(bits_num);

	for (int i = 0; i < bits_num; i++)
	{
		bit_flips[bits_num * path + i] = best_idx[i];
		pen[i] = sat_m<R>(std::abs(l[array][off_l + bit_flips[bits_num * path + i]]));
	} 

	switch (L)
	{
	case 2:
		metrics_vec[1][c_num * path + 0] = metrics [path]; // empty
		metrics_vec[1][c_num * path + 1] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[0]); // e+a0
		break;
	case 4:
		metrics_vec[1][c_num * path + 0] = metrics [path]; // empty
		metrics_vec[1][c_num * path + 1] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[0]); // e+a0
		metrics_vec[1][c_num * path + 2] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[1]); // e+a1
		metrics_vec[1][c_num * path + 3] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[1]); // e+a0+a1
		if (n_elmts >= 4)
		{
			metrics_vec[1][c_num * path + 4] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[2]); // e+a2
		}
		break;
	case 8:		
		metrics_vec[1][c_num * path + 0] = metrics [path]; // empty
		metrics_vec[1][c_num * path + 1] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[0]); // e+a0
		metrics_vec[1][c_num * path + 2] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[1]); // e+a1
		metrics_vec[1][c_num * path + 3] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[1]); // e+a0+a1
		if (n_elmts >= 4)
		{
			metrics_vec[1][c_num * path + 4] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[2]); // e+a2
			metrics_vec[1][c_num * path + 5] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[3]); // e+a3
			metrics_vec[1][c_num * path + 6] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[2]); // e+a0+a2
			metrics_vec[1][c_num * path + 7] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[3]); // e+a0+a3
			metrics_vec[1][c_num * path + 8] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[2]); // e+a1+a2
			metrics_vec[1][c_num * path + 9] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[2]); // e+a0+a1+a2
			if (n_elmts >= 8)
			{
				metrics_vec[1][c_num * path + 10] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[4]); // e+a4
				metrics_vec[1][c_num * path + 11] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[5]); // e+a5
				metrics_vec[1][c_num * path + 12] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[6]); // e+a6
			}
		}
		break;
	case 16:
		metrics_vec[1][c_num * path + 0] = metrics [path]; // empty
		metrics_vec[1][c_num * path + 1] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[0]); // e+a0
		metrics_vec[1][c_num * path + 2] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[1]); // e+a1
		metrics_vec[1][c_num * path + 3] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[1]); // e+a0+a1
		if (n_elmts >= 4)
		{
			metrics_vec[1][c_num * path + 4] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[2]); // e+a2
			metrics_vec[1][c_num * path + 5] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[3]); // e+a3
			metrics_vec[1][c_num * path + 6] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[2]); // e+a0+a2
			metrics_vec[1][c_num * path + 7] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[3]); // e+a0+a3
			metrics_vec[1][c_num * path + 8] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[2]); // e+a1+a2
			metrics_vec[1][c_num * path + 9] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[3]); // e+a1+a3
			metrics_vec[1][c_num * path + 10] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[3]); // e+a2+a3
			metrics_vec[1][c_num * path + 11] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[2]); // e+a0+a1+a2
			metrics_vec[1][c_num * path + 12] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[3]); // e+a0+a1+a3
			metrics_vec[1][c_num * path + 13] = sat_m<R>(metrics_vec[1][c_num * path + 6] + pen[3]); // e+a0+a2+a3
			metrics_vec[1][c_num * path + 14] = sat_m<R>(metrics_vec[1][c_num * path + 8] + pen[3]); // e+a1+a2+a3
			metrics_vec[1][c_num * path + 15] = sat_m<R>(metrics_vec[1][c_num * path + 11] + pen[3]); // e+a0+a1+a2+a3
			if (n_elmts >= 8)
			{
				metrics_vec[1][c_num * path + 16] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[4]); // e+a4
				metrics_vec[1][c_num * path + 17] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[5]); // e+a5
				metrics_vec[1][c_num * path + 18] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[6]); // e+a6
				metrics_vec[1][c_num * path + 19] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[7]); // e+a7
				metrics_vec[1][c_num * path + 20] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[4]); // e+a0+a4
				metrics_vec[1][c_num * path + 21] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[5]); // e+a0+a5
				metrics_vec[1][c_num * path + 22] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[6]); // e+a0+a6
				metrics_vec[1][c_num * path + 23] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[7]); // e+a0+a7
				metrics_vec[1][c_num * path + 24] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[4]); // e+a1+a4
				metrics_vec[1][c_num * path + 25] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[5]); // e+a1+a5
				metrics_vec[1][c_num * path + 26] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[4]); // e+a2+a4
				metrics_vec[1][c_num * path + 27] = sat_m<R>(metrics_vec[1][c_num * path + 5] + pen[4]); // e+a3+a4
				metrics_vec[1][c_num * path + 28] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[4]); // e+a0+a1+a4
				if (n_elmts >= 16)
				{
					metrics_vec[1][c_num * path + 29] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[8]); // e+a8
					metrics_vec[1][c_num * path + 30] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[9]); // e+a9
					metrics_vec[1][c_num * path + 31] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[10]); // e+a10
					metrics_vec[1][c_num * path + 32] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[11]); // e+a11
					metrics_vec[1][c_num * path + 33] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[12]); // e+a12
					metrics_vec[1][c_num * path + 34] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[13]); // e+a13
					metrics_vec[1][c_num * path + 35] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[14]); // e+a14

				}
			}
		}
		break;
	case 32:
		metrics_vec[1][c_num * path + 0] = metrics [path]; // empty
		metrics_vec[1][c_num * path + 1] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[0]); // e+a0
		metrics_vec[1][c_num * path + 2] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[1]); // e+a1
		metrics_vec[1][c_num * path + 3] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[1]); // e+a0+a1
		if (n_elmts >= 4)
		{
			metrics_vec[1][c_num * path + 4] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[2]); // e+a2
			metrics_vec[1][c_num * path + 5] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[3]); // e+a3
			metrics_vec[1][c_num * path + 6] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[2]); // e+a0+a2
			metrics_vec[1][c_num * path + 7] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[3]); // e+a0+a3
			metrics_vec[1][c_num * path + 8] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[2]); // e+a1+a2
			metrics_vec[1][c_num * path + 9] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[3]); // e+a1+a3
			metrics_vec[1][c_num * path + 10] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[3]); // e+a2+a3
			metrics_vec[1][c_num * path + 11] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[2]); // e+a0+a1+a2
			metrics_vec[1][c_num * path + 12] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[3]); // e+a0+a1+a3
			metrics_vec[1][c_num * path + 13] = sat_m<R>(metrics_vec[1][c_num * path + 6] + pen[3]); // e+a0+a2+a3
			metrics_vec[1][c_num * path + 14] = sat_m<R>(metrics_vec[1][c_num * path + 8] + pen[3]); // e+a1+a2+a3
			metrics_vec[1][c_num * path + 15] = sat_m<R>(metrics_vec[1][c_num * path + 11] + pen[3]); // e+a0+a1+a2+a3
			if (n_elmts >= 8)
			{
				metrics_vec[1][c_num * path + 16] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[4]); // e+a4
				metrics_vec[1][c_num * path + 17] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[5]); // e+a5
				metrics_vec[1][c_num * path + 18] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[6]); // e+a6
				metrics_vec[1][c_num * path + 19] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[7]); // e+a7
				metrics_vec[1][c_num * path + 20] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[4]); // e+a0+a4
				metrics_vec[1][c_num * path + 21] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[5]); // e+a0+a5
				metrics_vec[1][c_num * path + 22] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[6]); // e+a0+a6
				metrics_vec[1][c_num * path + 23] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[7]); // e+a0+a7
				metrics_vec[1][c_num * path + 24] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[4]); // e+a1+a4
				metrics_vec[1][c_num * path + 25] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[5]); // e+a1+a5
				metrics_vec[1][c_num * path + 26] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[6]); // e+a1+a6
				metrics_vec[1][c_num * path + 27] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[7]); // e+a1+a7
				metrics_vec[1][c_num * path + 28] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[4]); // e+a2+a4
				metrics_vec[1][c_num * path + 29] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[5]); // e+a2+a5
				metrics_vec[1][c_num * path + 30] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[6]); // e+a2+a6
				metrics_vec[1][c_num * path + 31] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[7]); // e+a2+a7
				metrics_vec[1][c_num * path + 32] = sat_m<R>(metrics_vec[1][c_num * path + 5] + pen[4]); // e+a3+a4
				metrics_vec[1][c_num * path + 33] = sat_m<R>(metrics_vec[1][c_num * path + 5] + pen[5]); // e+a3+a5
				metrics_vec[1][c_num * path + 34] = sat_m<R>(metrics_vec[1][c_num * path + 5] + pen[6]); // e+a3+a6
				metrics_vec[1][c_num * path + 35] = sat_m<R>(metrics_vec[1][c_num * path + 5] + pen[7]); // e+a3+a7
				metrics_vec[1][c_num * path + 36] = sat_m<R>(metrics_vec[1][c_num * path + 16] + pen[5]); // e+a4+a5
				metrics_vec[1][c_num * path + 37] = sat_m<R>(metrics_vec[1][c_num * path + 16] + pen[6]); // e+a4+a6
				metrics_vec[1][c_num * path + 38] = sat_m<R>(metrics_vec[1][c_num * path + 17] + pen[6]); // e+a5+a6
				metrics_vec[1][c_num * path + 39] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[4]); // e+a0+a1+a4
				metrics_vec[1][c_num * path + 40] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[5]); // e+a0+a1+a5
				metrics_vec[1][c_num * path + 41] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[6]); // e+a0+a1+a6
				metrics_vec[1][c_num * path + 42] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[7]); // e+a0+a1+a7
				metrics_vec[1][c_num * path + 43] = sat_m<R>(metrics_vec[1][c_num * path + 6] + pen[4]); // e+a0+a2+a4
				metrics_vec[1][c_num * path + 44] = sat_m<R>(metrics_vec[1][c_num * path + 6] + pen[5]); // e+a0+a2+a5
				metrics_vec[1][c_num * path + 45] = sat_m<R>(metrics_vec[1][c_num * path + 6] + pen[6]); // e+a0+a2+a6
				metrics_vec[1][c_num * path + 46] = sat_m<R>(metrics_vec[1][c_num * path + 7] + pen[4]); // e+a0+a3+a4
				metrics_vec[1][c_num * path + 47] = sat_m<R>(metrics_vec[1][c_num * path + 7] + pen[5]); // e+a0+a3+a5
				metrics_vec[1][c_num * path + 48] = sat_m<R>(metrics_vec[1][c_num * path + 20] + pen[5]); // e+a0+a4+a5
				metrics_vec[1][c_num * path + 49] = sat_m<R>(metrics_vec[1][c_num * path + 8] + pen[4]); // e+a1+a2+a4
				metrics_vec[1][c_num * path + 50] = sat_m<R>(metrics_vec[1][c_num * path + 8] + pen[5]); // e+a1+a2+a5
				metrics_vec[1][c_num * path + 51] = sat_m<R>(metrics_vec[1][c_num * path + 9] + pen[4]); // e+a1+a3+a4
				metrics_vec[1][c_num * path + 52] = sat_m<R>(metrics_vec[1][c_num * path + 10] + pen[4]); // e+a2+a3+a4
				metrics_vec[1][c_num * path + 53] = sat_m<R>(metrics_vec[1][c_num * path + 11] + pen[4]); // e+a0+a1+a2+a4
				metrics_vec[1][c_num * path + 54] = sat_m<R>(metrics_vec[1][c_num * path + 11] + pen[5]); // e+a0+a1+a2+a5
				metrics_vec[1][c_num * path + 55] = sat_m<R>(metrics_vec[1][c_num * path + 12] + pen[4]); // e+a0+a1+a3+a4
				metrics_vec[1][c_num * path + 56] = sat_m<R>(metrics_vec[1][c_num * path + 13] + pen[4]); // e+a0+a2+a3+a4
				metrics_vec[1][c_num * path + 57] = sat_m<R>(metrics_vec[1][c_num * path + 14] + pen[4]); // e+a1+a2+a3+a4
				metrics_vec[1][c_num * path + 58] = sat_m<R>(metrics_vec[1][c_num * path + 15] + pen[4]); // e+a0+a1+a2+a3+a4
				if (n_elmts >= 16)
				{
					metrics_vec[1][c_num * path + 59] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[8]); // e+a8
					metrics_vec[1][c_num * path + 60] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[9]); // e+a9
					metrics_vec[1][c_num * path + 61] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[10]); // e+a10
					metrics_vec[1][c_num * path + 62] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[11]); // e+a11
					metrics_vec[1][c_num * path + 63] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[12]); // e+a12
					metrics_vec[1][c_num * path + 64] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[13]); // e+a13
					metrics_vec[1][c_num * path + 65] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[14]); // e+a14
					metrics_vec[1][c_num * path + 66] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[15]); // e+a15
					metrics_vec[1][c_num * path + 67] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[8]); // e+a0+a8
					metrics_vec[1][c_num * path + 68] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[9]); // e+a0+a9
					metrics_vec[1][c_num * path + 69] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[10]); // e+a0+a10
					metrics_vec[1][c_num * path + 70] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[11]); // e+a0+a11
					metrics_vec[1][c_num * path + 71] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[12]); // e+a0+a12
					metrics_vec[1][c_num * path + 72] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[13]); // e+a0+a13
					metrics_vec[1][c_num * path + 73] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[14]); // e+a0+a14
					metrics_vec[1][c_num * path + 74] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[15]); // e+a0+a15
					metrics_vec[1][c_num * path + 75] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[8]); // e+a1+a8
					metrics_vec[1][c_num * path + 76] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[9]); // e+a1+a9
					metrics_vec[1][c_num * path + 77] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[10]); // e+a1+a10
					metrics_vec[1][c_num * path + 78] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[8]); // e+a2+a8
					metrics_vec[1][c_num * path + 79] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[8]); // e+a0+a1+a8
					if (n_elmts >= 32)
					{
						metrics_vec[1][c_num * path + 80] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[16]); // e+a16
						metrics_vec[1][c_num * path + 81] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[17]); // e+a17
						metrics_vec[1][c_num * path + 82] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[18]); // e+a18
						metrics_vec[1][c_num * path + 83] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[19]); // e+a19
						metrics_vec[1][c_num * path + 84] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[20]); // e+a20
						metrics_vec[1][c_num * path + 85] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[21]); // e+a21
						metrics_vec[1][c_num * path + 86] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[22]); // e+a22
						metrics_vec[1][c_num * path + 87] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[23]); // e+a23
						metrics_vec[1][c_num * path + 88] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[24]); // e+a24
						metrics_vec[1][c_num * path + 89] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[25]); // e+a25
						metrics_vec[1][c_num * path + 90] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[26]); // e+a26
						metrics_vec[1][c_num * path + 91] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[27]); // e+a27
						metrics_vec[1][c_num * path + 92] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[28]); // e+a28
						metrics_vec[1][c_num * path + 93] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[29]); // e+a29
						metrics_vec[1][c_num * path + 94] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[30]); // e+a30
					}
				}
			}
		}
		break;
	case 64:
		metrics_vec[1][c_num * path + 0] = metrics [path]; // empty
		metrics_vec[1][c_num * path + 1] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[0]); // e+a0
		metrics_vec[1][c_num * path + 2] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[1]); // e+a1
		metrics_vec[1][c_num * path + 3] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[1]); // e+a0+a1
		if (n_elmts >= 4)
		{
			metrics_vec[1][c_num * path + 4] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[2]); // e+a2
			metrics_vec[1][c_num * path + 5] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[3]); // e+a3
			metrics_vec[1][c_num * path + 6] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[2]); // e+a0+a2
			metrics_vec[1][c_num * path + 7] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[3]); // e+a0+a3
			metrics_vec[1][c_num * path + 8] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[2]); // e+a1+a2
			metrics_vec[1][c_num * path + 9] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[3]); // e+a1+a3
			metrics_vec[1][c_num * path + 10] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[3]); // e+a2+a3
			metrics_vec[1][c_num * path + 11] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[2]); // e+a0+a1+a2
			metrics_vec[1][c_num * path + 12] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[3]); // e+a0+a1+a3
			metrics_vec[1][c_num * path + 13] = sat_m<R>(metrics_vec[1][c_num * path + 6] + pen[3]); // e+a0+a2+a3
			metrics_vec[1][c_num * path + 14] = sat_m<R>(metrics_vec[1][c_num * path + 8] + pen[3]); // e+a1+a2+a3
			metrics_vec[1][c_num * path + 15] = sat_m<R>(metrics_vec[1][c_num * path + 11] + pen[3]); // e+a0+a1+a2+a3
			if (n_elmts >= 8)
			{
				metrics_vec[1][c_num * path + 16] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[4]); // e+a4
				metrics_vec[1][c_num * path + 17] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[5]); // e+a5
				metrics_vec[1][c_num * path + 18] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[6]); // e+a6
				metrics_vec[1][c_num * path + 19] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[7]); // e+a7
				metrics_vec[1][c_num * path + 20] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[4]); // e+a0+a4
				metrics_vec[1][c_num * path + 21] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[5]); // e+a0+a5
				metrics_vec[1][c_num * path + 22] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[6]); // e+a0+a6
				metrics_vec[1][c_num * path + 23] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[7]); // e+a0+a7
				metrics_vec[1][c_num * path + 24] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[4]); // e+a1+a4
				metrics_vec[1][c_num * path + 25] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[5]); // e+a1+a5
				metrics_vec[1][c_num * path + 26] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[6]); // e+a1+a6
				metrics_vec[1][c_num * path + 27] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[7]); // e+a1+a7
				metrics_vec[1][c_num * path + 28] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[4]); // e+a2+a4
				metrics_vec[1][c_num * path + 29] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[5]); // e+a2+a5
				metrics_vec[1][c_num * path + 30] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[6]); // e+a2+a6
				metrics_vec[1][c_num * path + 31] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[7]); // e+a2+a7
				metrics_vec[1][c_num * path + 32] = sat_m<R>(metrics_vec[1][c_num * path + 5] + pen[4]); // e+a3+a4
				metrics_vec[1][c_num * path + 33] = sat_m<R>(metrics_vec[1][c_num * path + 5] + pen[5]); // e+a3+a5
				metrics_vec[1][c_num * path + 34] = sat_m<R>(metrics_vec[1][c_num * path + 5] + pen[6]); // e+a3+a6
				metrics_vec[1][c_num * path + 35] = sat_m<R>(metrics_vec[1][c_num * path + 5] + pen[7]); // e+a3+a7
				metrics_vec[1][c_num * path + 36] = sat_m<R>(metrics_vec[1][c_num * path + 16] + pen[5]); // e+a4+a5
				metrics_vec[1][c_num * path + 37] = sat_m<R>(metrics_vec[1][c_num * path + 16] + pen[6]); // e+a4+a6
				metrics_vec[1][c_num * path + 38] = sat_m<R>(metrics_vec[1][c_num * path + 16] + pen[7]); // e+a4+a7
				metrics_vec[1][c_num * path + 39] = sat_m<R>(metrics_vec[1][c_num * path + 17] + pen[6]); // e+a5+a6
				metrics_vec[1][c_num * path + 40] = sat_m<R>(metrics_vec[1][c_num * path + 17] + pen[7]); // e+a5+a7
				metrics_vec[1][c_num * path + 41] = sat_m<R>(metrics_vec[1][c_num * path + 18] + pen[7]); // e+a6+a7
				metrics_vec[1][c_num * path + 42] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[4]); // e+a0+a1+a4
				metrics_vec[1][c_num * path + 43] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[5]); // e+a0+a1+a5
				metrics_vec[1][c_num * path + 44] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[6]); // e+a0+a1+a6
				metrics_vec[1][c_num * path + 45] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[7]); // e+a0+a1+a7
				metrics_vec[1][c_num * path + 46] = sat_m<R>(metrics_vec[1][c_num * path + 6] + pen[4]); // e+a0+a2+a4
				metrics_vec[1][c_num * path + 47] = sat_m<R>(metrics_vec[1][c_num * path + 6] + pen[5]); // e+a0+a2+a5
				metrics_vec[1][c_num * path + 48] = sat_m<R>(metrics_vec[1][c_num * path + 6] + pen[6]); // e+a0+a2+a6
				metrics_vec[1][c_num * path + 49] = sat_m<R>(metrics_vec[1][c_num * path + 6] + pen[7]); // e+a0+a2+a7
				metrics_vec[1][c_num * path + 50] = sat_m<R>(metrics_vec[1][c_num * path + 7] + pen[4]); // e+a0+a3+a4
				metrics_vec[1][c_num * path + 51] = sat_m<R>(metrics_vec[1][c_num * path + 7] + pen[5]); // e+a0+a3+a5
				metrics_vec[1][c_num * path + 52] = sat_m<R>(metrics_vec[1][c_num * path + 7] + pen[6]); // e+a0+a3+a6
				metrics_vec[1][c_num * path + 53] = sat_m<R>(metrics_vec[1][c_num * path + 7] + pen[7]); // e+a0+a3+a7
				metrics_vec[1][c_num * path + 54] = sat_m<R>(metrics_vec[1][c_num * path + 20] + pen[5]); // e+a0+a4+a5
				metrics_vec[1][c_num * path + 55] = sat_m<R>(metrics_vec[1][c_num * path + 20] + pen[6]); // e+a0+a4+a6
				metrics_vec[1][c_num * path + 56] = sat_m<R>(metrics_vec[1][c_num * path + 20] + pen[7]); // e+a0+a4+a7
				metrics_vec[1][c_num * path + 57] = sat_m<R>(metrics_vec[1][c_num * path + 21] + pen[6]); // e+a0+a5+a6
				metrics_vec[1][c_num * path + 58] = sat_m<R>(metrics_vec[1][c_num * path + 21] + pen[7]); // e+a0+a5+a7
				metrics_vec[1][c_num * path + 59] = sat_m<R>(metrics_vec[1][c_num * path + 22] + pen[7]); // e+a0+a6+a7
				metrics_vec[1][c_num * path + 60] = sat_m<R>(metrics_vec[1][c_num * path + 8] + pen[4]); // e+a1+a2+a4
				metrics_vec[1][c_num * path + 61] = sat_m<R>(metrics_vec[1][c_num * path + 8] + pen[5]); // e+a1+a2+a5
				metrics_vec[1][c_num * path + 62] = sat_m<R>(metrics_vec[1][c_num * path + 8] + pen[6]); // e+a1+a2+a6
				metrics_vec[1][c_num * path + 63] = sat_m<R>(metrics_vec[1][c_num * path + 8] + pen[7]); // e+a1+a2+a7
				metrics_vec[1][c_num * path + 64] = sat_m<R>(metrics_vec[1][c_num * path + 9] + pen[4]); // e+a1+a3+a4
				metrics_vec[1][c_num * path + 65] = sat_m<R>(metrics_vec[1][c_num * path + 9] + pen[5]); // e+a1+a3+a5
				metrics_vec[1][c_num * path + 66] = sat_m<R>(metrics_vec[1][c_num * path + 9] + pen[6]); // e+a1+a3+a6
				metrics_vec[1][c_num * path + 67] = sat_m<R>(metrics_vec[1][c_num * path + 9] + pen[7]); // e+a1+a3+a7
				metrics_vec[1][c_num * path + 68] = sat_m<R>(metrics_vec[1][c_num * path + 24] + pen[5]); // e+a1+a4+a5
				metrics_vec[1][c_num * path + 69] = sat_m<R>(metrics_vec[1][c_num * path + 24] + pen[6]); // e+a1+a4+a6
				metrics_vec[1][c_num * path + 70] = sat_m<R>(metrics_vec[1][c_num * path + 24] + pen[7]); // e+a1+a4+a7
				metrics_vec[1][c_num * path + 71] = sat_m<R>(metrics_vec[1][c_num * path + 25] + pen[6]); // e+a1+a5+a6
				metrics_vec[1][c_num * path + 72] = sat_m<R>(metrics_vec[1][c_num * path + 10] + pen[4]); // e+a2+a3+a4
				metrics_vec[1][c_num * path + 73] = sat_m<R>(metrics_vec[1][c_num * path + 10] + pen[5]); // e+a2+a3+a5
				metrics_vec[1][c_num * path + 74] = sat_m<R>(metrics_vec[1][c_num * path + 10] + pen[6]); // e+a2+a3+a6
				metrics_vec[1][c_num * path + 75] = sat_m<R>(metrics_vec[1][c_num * path + 10] + pen[7]); // e+a2+a3+a7
				metrics_vec[1][c_num * path + 76] = sat_m<R>(metrics_vec[1][c_num * path + 28] + pen[5]); // e+a2+a4+a5
				metrics_vec[1][c_num * path + 77] = sat_m<R>(metrics_vec[1][c_num * path + 28] + pen[6]); // e+a2+a4+a6
				metrics_vec[1][c_num * path + 78] = sat_m<R>(metrics_vec[1][c_num * path + 29] + pen[6]); // e+a2+a5+a6
				metrics_vec[1][c_num * path + 79] = sat_m<R>(metrics_vec[1][c_num * path + 32] + pen[5]); // e+a3+a4+a5
				metrics_vec[1][c_num * path + 80] = sat_m<R>(metrics_vec[1][c_num * path + 32] + pen[6]); // e+a3+a4+a6
				metrics_vec[1][c_num * path + 81] = sat_m<R>(metrics_vec[1][c_num * path + 33] + pen[6]); // e+a3+a5+a6
				metrics_vec[1][c_num * path + 82] = sat_m<R>(metrics_vec[1][c_num * path + 36] + pen[6]); // e+a4+a5+a6
				metrics_vec[1][c_num * path + 83] = sat_m<R>(metrics_vec[1][c_num * path + 11] + pen[4]); // e+a0+a1+a2+a4
				metrics_vec[1][c_num * path + 84] = sat_m<R>(metrics_vec[1][c_num * path + 11] + pen[5]); // e+a0+a1+a2+a5
				metrics_vec[1][c_num * path + 85] = sat_m<R>(metrics_vec[1][c_num * path + 11] + pen[6]); // e+a0+a1+a2+a6
				metrics_vec[1][c_num * path + 86] = sat_m<R>(metrics_vec[1][c_num * path + 11] + pen[7]); // e+a0+a1+a2+a7
				metrics_vec[1][c_num * path + 87] = sat_m<R>(metrics_vec[1][c_num * path + 12] + pen[4]); // e+a0+a1+a3+a4
				metrics_vec[1][c_num * path + 88] = sat_m<R>(metrics_vec[1][c_num * path + 12] + pen[5]); // e+a0+a1+a3+a5
				metrics_vec[1][c_num * path + 89] = sat_m<R>(metrics_vec[1][c_num * path + 12] + pen[6]); // e+a0+a1+a3+a6
				metrics_vec[1][c_num * path + 90] = sat_m<R>(metrics_vec[1][c_num * path + 12] + pen[7]); // e+a0+a1+a3+a7
				metrics_vec[1][c_num * path + 91] = sat_m<R>(metrics_vec[1][c_num * path + 42] + pen[5]); // e+a0+a1+a4+a5
				metrics_vec[1][c_num * path + 92] = sat_m<R>(metrics_vec[1][c_num * path + 42] + pen[6]); // e+a0+a1+a4+a6
				metrics_vec[1][c_num * path + 93] = sat_m<R>(metrics_vec[1][c_num * path + 43] + pen[6]); // e+a0+a1+a5+a6
				metrics_vec[1][c_num * path + 94] = sat_m<R>(metrics_vec[1][c_num * path + 13] + pen[4]); // e+a0+a2+a3+a4
				metrics_vec[1][c_num * path + 95] = sat_m<R>(metrics_vec[1][c_num * path + 13] + pen[5]); // e+a0+a2+a3+a5
				metrics_vec[1][c_num * path + 96] = sat_m<R>(metrics_vec[1][c_num * path + 13] + pen[6]); // e+a0+a2+a3+a6
				metrics_vec[1][c_num * path + 97] = sat_m<R>(metrics_vec[1][c_num * path + 46] + pen[5]); // e+a0+a2+a4+a5
				metrics_vec[1][c_num * path + 98] = sat_m<R>(metrics_vec[1][c_num * path + 50] + pen[5]); // e+a0+a3+a4+a5
				metrics_vec[1][c_num * path + 99] = sat_m<R>(metrics_vec[1][c_num * path + 14] + pen[4]); // e+a1+a2+a3+a4
				metrics_vec[1][c_num * path + 100] = sat_m<R>(metrics_vec[1][c_num * path + 14] + pen[5]); // e+a1+a2+a3+a5
				metrics_vec[1][c_num * path + 101] = sat_m<R>(metrics_vec[1][c_num * path + 14] + pen[6]); // e+a1+a2+a3+a6
				metrics_vec[1][c_num * path + 102] = sat_m<R>(metrics_vec[1][c_num * path + 60] + pen[5]); // e+a1+a2+a4+a5
				metrics_vec[1][c_num * path + 103] = sat_m<R>(metrics_vec[1][c_num * path + 64] + pen[5]); // e+a1+a3+a4+a5
				metrics_vec[1][c_num * path + 104] = sat_m<R>(metrics_vec[1][c_num * path + 72] + pen[5]); // e+a2+a3+a4+a5
				metrics_vec[1][c_num * path + 105] = sat_m<R>(metrics_vec[1][c_num * path + 15] + pen[4]); // e+a0+a1+a2+a3+a4
				metrics_vec[1][c_num * path + 106] = sat_m<R>(metrics_vec[1][c_num * path + 15] + pen[5]); // e+a0+a1+a2+a3+a5
				metrics_vec[1][c_num * path + 107] = sat_m<R>(metrics_vec[1][c_num * path + 15] + pen[6]); // e+a0+a1+a2+a3+a6
				metrics_vec[1][c_num * path + 108] = sat_m<R>(metrics_vec[1][c_num * path + 83] + pen[5]); // e+a0+a1+a2+a4+a5
				metrics_vec[1][c_num * path + 109] = sat_m<R>(metrics_vec[1][c_num * path + 87] + pen[5]); // e+a0+a1+a3+a4+a5
				metrics_vec[1][c_num * path + 110] = sat_m<R>(metrics_vec[1][c_num * path + 94] + pen[5]); // e+a0+a2+a3+a4+a5
				metrics_vec[1][c_num * path + 111] = sat_m<R>(metrics_vec[1][c_num * path + 99] + pen[5]); // e+a1+a2+a3+a4+a5
				if (n_elmts >= 16)
				{
					metrics_vec[1][c_num * path + 112] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[8]); // e+a8
					metrics_vec[1][c_num * path + 113] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[9]); // e+a9
					metrics_vec[1][c_num * path + 114] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[10]); // e+a10
					metrics_vec[1][c_num * path + 115] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[11]); // e+a11
					metrics_vec[1][c_num * path + 116] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[12]); // e+a12
					metrics_vec[1][c_num * path + 117] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[13]); // e+a13
					metrics_vec[1][c_num * path + 118] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[14]); // e+a14
					metrics_vec[1][c_num * path + 119] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[15]); // e+a15
					metrics_vec[1][c_num * path + 120] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[8]); // e+a0+a8
					metrics_vec[1][c_num * path + 121] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[9]); // e+a0+a9
					metrics_vec[1][c_num * path + 122] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[10]); // e+a0+a10
					metrics_vec[1][c_num * path + 123] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[11]); // e+a0+a11
					metrics_vec[1][c_num * path + 124] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[12]); // e+a0+a12
					metrics_vec[1][c_num * path + 125] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[13]); // e+a0+a13
					metrics_vec[1][c_num * path + 126] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[14]); // e+a0+a14
					metrics_vec[1][c_num * path + 127] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[15]); // e+a0+a15
					metrics_vec[1][c_num * path + 128] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[8]); // e+a1+a8
					metrics_vec[1][c_num * path + 129] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[9]); // e+a1+a9
					metrics_vec[1][c_num * path + 130] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[10]); // e+a1+a10
					metrics_vec[1][c_num * path + 131] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[11]); // e+a1+a11
					metrics_vec[1][c_num * path + 132] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[12]); // e+a1+a12
					metrics_vec[1][c_num * path + 133] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[13]); // e+a1+a13
					metrics_vec[1][c_num * path + 134] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[14]); // e+a1+a14
					metrics_vec[1][c_num * path + 135] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[15]); // e+a1+a15
					metrics_vec[1][c_num * path + 136] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[8]); // e+a2+a8
					metrics_vec[1][c_num * path + 137] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[9]); // e+a2+a9
					metrics_vec[1][c_num * path + 138] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[10]); // e+a2+a10
					metrics_vec[1][c_num * path + 139] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[11]); // e+a2+a11
					metrics_vec[1][c_num * path + 140] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[12]); // e+a2+a12
					metrics_vec[1][c_num * path + 141] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[13]); // e+a2+a13
					metrics_vec[1][c_num * path + 142] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[14]); // e+a2+a14
					metrics_vec[1][c_num * path + 143] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[15]); // e+a2+a15
					metrics_vec[1][c_num * path + 144] = sat_m<R>(metrics_vec[1][c_num * path + 5] + pen[8]); // e+a3+a8
					metrics_vec[1][c_num * path + 145] = sat_m<R>(metrics_vec[1][c_num * path + 5] + pen[9]); // e+a3+a9
					metrics_vec[1][c_num * path + 146] = sat_m<R>(metrics_vec[1][c_num * path + 5] + pen[10]); // e+a3+a10
					metrics_vec[1][c_num * path + 147] = sat_m<R>(metrics_vec[1][c_num * path + 5] + pen[11]); // e+a3+a11
					metrics_vec[1][c_num * path + 148] = sat_m<R>(metrics_vec[1][c_num * path + 5] + pen[12]); // e+a3+a12
					metrics_vec[1][c_num * path + 149] = sat_m<R>(metrics_vec[1][c_num * path + 5] + pen[13]); // e+a3+a13
					metrics_vec[1][c_num * path + 150] = sat_m<R>(metrics_vec[1][c_num * path + 16] + pen[8]); // e+a4+a8
					metrics_vec[1][c_num * path + 151] = sat_m<R>(metrics_vec[1][c_num * path + 16] + pen[9]); // e+a4+a9
					metrics_vec[1][c_num * path + 152] = sat_m<R>(metrics_vec[1][c_num * path + 16] + pen[10]); // e+a4+a10
					metrics_vec[1][c_num * path + 153] = sat_m<R>(metrics_vec[1][c_num * path + 16] + pen[11]); // e+a4+a11
					metrics_vec[1][c_num * path + 154] = sat_m<R>(metrics_vec[1][c_num * path + 16] + pen[12]); // e+a4+a12
					metrics_vec[1][c_num * path + 155] = sat_m<R>(metrics_vec[1][c_num * path + 17] + pen[8]); // e+a5+a8
					metrics_vec[1][c_num * path + 156] = sat_m<R>(metrics_vec[1][c_num * path + 17] + pen[9]); // e+a5+a9
					metrics_vec[1][c_num * path + 157] = sat_m<R>(metrics_vec[1][c_num * path + 17] + pen[10]); // e+a5+a10
					metrics_vec[1][c_num * path + 158] = sat_m<R>(metrics_vec[1][c_num * path + 17] + pen[11]); // e+a5+a11
					metrics_vec[1][c_num * path + 159] = sat_m<R>(metrics_vec[1][c_num * path + 18] + pen[8]); // e+a6+a8
					metrics_vec[1][c_num * path + 160] = sat_m<R>(metrics_vec[1][c_num * path + 18] + pen[9]); // e+a6+a9
					metrics_vec[1][c_num * path + 161] = sat_m<R>(metrics_vec[1][c_num * path + 18] + pen[10]); // e+a6+a10
					metrics_vec[1][c_num * path + 162] = sat_m<R>(metrics_vec[1][c_num * path + 19] + pen[8]); // e+a7+a8
					metrics_vec[1][c_num * path + 163] = sat_m<R>(metrics_vec[1][c_num * path + 19] + pen[9]); // e+a7+a9
					metrics_vec[1][c_num * path + 164] = sat_m<R>(metrics_vec[1][c_num * path + 19] + pen[10]); // e+a7+a10
					metrics_vec[1][c_num * path + 165] = sat_m<R>(metrics_vec[1][c_num * path + 112] + pen[9]); // e+a8+a9
					metrics_vec[1][c_num * path + 166] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[8]); // e+a0+a1+a8
					metrics_vec[1][c_num * path + 167] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[9]); // e+a0+a1+a9
					metrics_vec[1][c_num * path + 168] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[10]); // e+a0+a1+a10
					metrics_vec[1][c_num * path + 169] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[11]); // e+a0+a1+a11
					metrics_vec[1][c_num * path + 170] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[12]); // e+a0+a1+a12
					metrics_vec[1][c_num * path + 171] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[13]); // e+a0+a1+a13
					metrics_vec[1][c_num * path + 172] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[14]); // e+a0+a1+a14
					metrics_vec[1][c_num * path + 173] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[15]); // e+a0+a1+a15
					metrics_vec[1][c_num * path + 174] = sat_m<R>(metrics_vec[1][c_num * path + 6] + pen[8]); // e+a0+a2+a8
					metrics_vec[1][c_num * path + 175] = sat_m<R>(metrics_vec[1][c_num * path + 6] + pen[9]); // e+a0+a2+a9
					metrics_vec[1][c_num * path + 176] = sat_m<R>(metrics_vec[1][c_num * path + 6] + pen[10]); // e+a0+a2+a10
					metrics_vec[1][c_num * path + 177] = sat_m<R>(metrics_vec[1][c_num * path + 6] + pen[11]); // e+a0+a2+a11
					metrics_vec[1][c_num * path + 178] = sat_m<R>(metrics_vec[1][c_num * path + 7] + pen[8]); // e+a0+a3+a8
					metrics_vec[1][c_num * path + 179] = sat_m<R>(metrics_vec[1][c_num * path + 7] + pen[9]); // e+a0+a3+a9
					metrics_vec[1][c_num * path + 180] = sat_m<R>(metrics_vec[1][c_num * path + 20] + pen[8]); // e+a0+a4+a8
					metrics_vec[1][c_num * path + 181] = sat_m<R>(metrics_vec[1][c_num * path + 8] + pen[8]); // e+a1+a2+a8
					metrics_vec[1][c_num * path + 182] = sat_m<R>(metrics_vec[1][c_num * path + 8] + pen[9]); // e+a1+a2+a9
					metrics_vec[1][c_num * path + 183] = sat_m<R>(metrics_vec[1][c_num * path + 8] + pen[10]); // e+a1+a2+a10
					metrics_vec[1][c_num * path + 184] = sat_m<R>(metrics_vec[1][c_num * path + 11] + pen[8]); // e+a0+a1+a2+a8
					metrics_vec[1][c_num * path + 185] = sat_m<R>(metrics_vec[1][c_num * path + 11] + pen[9]); // e+a0+a1+a2+a9
					if (n_elmts >= 32)
					{
						metrics_vec[1][c_num * path + 186] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[16]); // e+a16
						metrics_vec[1][c_num * path + 187] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[17]); // e+a17
						metrics_vec[1][c_num * path + 188] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[18]); // e+a18
						metrics_vec[1][c_num * path + 189] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[19]); // e+a19
						metrics_vec[1][c_num * path + 190] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[20]); // e+a20
						metrics_vec[1][c_num * path + 191] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[21]); // e+a21
						metrics_vec[1][c_num * path + 192] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[22]); // e+a22
						metrics_vec[1][c_num * path + 193] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[23]); // e+a23
						metrics_vec[1][c_num * path + 194] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[24]); // e+a24
						metrics_vec[1][c_num * path + 195] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[25]); // e+a25
						metrics_vec[1][c_num * path + 196] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[26]); // e+a26
						metrics_vec[1][c_num * path + 197] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[27]); // e+a27
						metrics_vec[1][c_num * path + 198] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[28]); // e+a28
						metrics_vec[1][c_num * path + 199] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[29]); // e+a29
						metrics_vec[1][c_num * path + 200] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[30]); // e+a30
						metrics_vec[1][c_num * path + 201] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[31]); // e+a31
						metrics_vec[1][c_num * path + 202] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[16]); // e+a0+a16
						metrics_vec[1][c_num * path + 203] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[17]); // e+a0+a17
						metrics_vec[1][c_num * path + 204] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[18]); // e+a0+a18
						metrics_vec[1][c_num * path + 205] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[19]); // e+a0+a19
						metrics_vec[1][c_num * path + 206] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[20]); // e+a0+a20
						metrics_vec[1][c_num * path + 207] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[21]); // e+a0+a21
						metrics_vec[1][c_num * path + 208] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[22]); // e+a0+a22
						metrics_vec[1][c_num * path + 209] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[23]); // e+a0+a23
						metrics_vec[1][c_num * path + 210] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[24]); // e+a0+a24
						metrics_vec[1][c_num * path + 211] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[25]); // e+a0+a25
						metrics_vec[1][c_num * path + 212] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[26]); // e+a0+a26
						metrics_vec[1][c_num * path + 213] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[27]); // e+a0+a27
						metrics_vec[1][c_num * path + 214] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[28]); // e+a0+a28
						metrics_vec[1][c_num * path + 215] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[29]); // e+a0+a29
						metrics_vec[1][c_num * path + 216] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[30]); // e+a0+a30
						metrics_vec[1][c_num * path + 217] = sat_m<R>(metrics_vec[1][c_num * path + 1] + pen[31]); // e+a0+a31
						metrics_vec[1][c_num * path + 218] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[16]); // e+a1+a16
						metrics_vec[1][c_num * path + 219] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[17]); // e+a1+a17
						metrics_vec[1][c_num * path + 220] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[18]); // e+a1+a18
						metrics_vec[1][c_num * path + 221] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[19]); // e+a1+a19
						metrics_vec[1][c_num * path + 222] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[20]); // e+a1+a20
						metrics_vec[1][c_num * path + 223] = sat_m<R>(metrics_vec[1][c_num * path + 2] + pen[21]); // e+a1+a21
						metrics_vec[1][c_num * path + 224] = sat_m<R>(metrics_vec[1][c_num * path + 4] + pen[16]); // e+a2+a16
						metrics_vec[1][c_num * path + 225] = sat_m<R>(metrics_vec[1][c_num * path + 3] + pen[16]); // e+a0+a1+a16
						if (n_elmts >= 64)
						{
							metrics_vec[1][c_num * path + 226] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[32]); // e+a32
							metrics_vec[1][c_num * path + 227] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[33]); // e+a33
							metrics_vec[1][c_num * path + 228] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[34]); // e+a34
							metrics_vec[1][c_num * path + 229] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[35]); // e+a35
							metrics_vec[1][c_num * path + 230] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[36]); // e+a36
							metrics_vec[1][c_num * path + 231] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[37]); // e+a37
							metrics_vec[1][c_num * path + 232] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[38]); // e+a38
							metrics_vec[1][c_num * path + 233] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[39]); // e+a39
							metrics_vec[1][c_num * path + 234] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[40]); // e+a40
							metrics_vec[1][c_num * path + 235] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[41]); // e+a41
							metrics_vec[1][c_num * path + 236] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[42]); // e+a42
							metrics_vec[1][c_num * path + 237] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[43]); // e+a43
							metrics_vec[1][c_num * path + 238] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[44]); // e+a44
							metrics_vec[1][c_num * path + 239] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[45]); // e+a45
							metrics_vec[1][c_num * path + 240] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[46]); // e+a46
							metrics_vec[1][c_num * path + 241] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[47]); // e+a47
							metrics_vec[1][c_num * path + 242] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[48]); // e+a48
							metrics_vec[1][c_num * path + 243] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[49]); // e+a49
							metrics_vec[1][c_num * path + 244] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[50]); // e+a50
							metrics_vec[1][c_num * path + 245] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[51]); // e+a51
							metrics_vec[1][c_num * path + 246] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[52]); // e+a52
							metrics_vec[1][c_num * path + 247] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[53]); // e+a53
							metrics_vec[1][c_num * path + 248] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[54]); // e+a54
							metrics_vec[1][c_num * path + 249] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[55]); // e+a55
							metrics_vec[1][c_num * path + 250] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[56]); // e+a56
							metrics_vec[1][c_num * path + 251] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[57]); // e+a57
							metrics_vec[1][c_num * path + 252] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[58]); // e+a58
							metrics_vec[1][c_num * path + 253] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[59]); // e+a59
							metrics_vec[1][c_num * path + 254] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[60]); // e+a60
							metrics_vec[1][c_num * path + 255] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[61]); // e+a61
							metrics_vec[1][c_num * path + 256] = sat_m<R>(metrics_vec[1][c_num * path + 0] + pen[62]); // e+a62
							metrics_vec[1][c_num * path + 257] = sat_m<R>(metrics_vec[1][c_num * path + 105] + pen[5]); // e+a0+a1+a2+a3+a4+a5
						}
					}
				}
			}
		}
		break;
	default:
		// Chase-II
		break;
	}				
}
template <typename B, typename R, class API_polar>
void Decoder_polar_SCL_mcfast_sys<B,R,API_polar>
::_update_spc_vec(const int path, const int array, const int n_elmts, const int bits_num, const int off_l)
{
	int c_num = spc_mc_size[(int)std::log2((n_elmts > L) ? L : n_elmts) - 1];
	std ::vector<R> pen(bits_num);
		
	for (auto i = 0; i < bits_num; i++)
	{
		bit_flips[bits_num * path + i] = best_idx[i];
		pen[i] = sat_m<R>(std::abs(l[array][off_l + bit_flips[bits_num * path + i]]));
	}
	
	auto sum = 0;

	for (auto i = 0; i < n_elmts; i++)
	{
		sum ^= (l[array][off_l + i] < 0);
	}	

	is_even[path] = (sum == 0);
	if (!is_even[path])
	{
		metrics[path] = metrics[path] + pen[0];
		pen[0] = -pen[0];
	}

	switch (L)
	{
	case 2:
		metrics_vec[2][c_num * path + 0] = metrics[path]; // empty
		metrics_vec[2][c_num * path + 1] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[1]); // e+a0+a1
		break;
	case 4:
		metrics_vec[2][c_num * path + 0] = metrics[path]; // empty
		metrics_vec[2][c_num * path + 1] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[1]); // e+a0+a1
		if (n_elmts >= 4)
		{
			metrics_vec[2][c_num * path + 2] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[2]); // e+a0+a2
			metrics_vec[2][c_num * path + 3] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[3]); // e+a0+a3
			metrics_vec[2][c_num * path + 4] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[2]); // e+a1+a2
		}
		break;
	case 8:
		metrics_vec[2][c_num * path + 0] = metrics[path]; // empty
		metrics_vec[2][c_num * path + 1] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[1]); // e+a0+a1
		if (n_elmts >= 4)
		{
			metrics_vec[2][c_num * path + 2] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[2]); // e+a0+a2
			metrics_vec[2][c_num * path + 3] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[3]); // e+a0+a3
			metrics_vec[2][c_num * path + 4] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[2]); // e+a1+a2
			metrics_vec[2][c_num * path + 5] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[3]); // e+a1+a3
			metrics_vec[2][c_num * path + 6] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[3]); // e+a2+a3
			metrics_vec[2][c_num * path + 7] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[3]); // e+a0+a1+a2+a3
			if (n_elmts >= 8)
			{
				metrics_vec[2][c_num * path + 8] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[4]); // e+a0+a4
				metrics_vec[2][c_num * path + 9] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[5]); // e+a0+a5
				metrics_vec[2][c_num * path + 10] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[6]); // e+a0+a6
				metrics_vec[2][c_num * path + 11] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[7]); // e+a0+a7
				metrics_vec[2][c_num * path + 12] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[4]); // e+a1+a4
			}
		}
		break;
	case 16:
		metrics_vec[2][c_num * path + 0] = metrics[path]; // empty
		metrics_vec[2][c_num * path + 1] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[1]); // e+a0+a1
		if (n_elmts >= 4)
		{
			metrics_vec[2][c_num * path + 2] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[2]); // e+a0+a2
			metrics_vec[2][c_num * path + 3] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[3]); // e+a0+a3
			metrics_vec[2][c_num * path + 4] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[2]); // e+a1+a2
			metrics_vec[2][c_num * path + 5] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[3]); // e+a1+a3
			metrics_vec[2][c_num * path + 6] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[3]); // e+a2+a3
			metrics_vec[2][c_num * path + 7] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[3]); // e+a0+a1+a2+a3
			if (n_elmts >= 8)
			{
				metrics_vec[2][c_num * path + 8] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[4]); // e+a0+a4
				metrics_vec[2][c_num * path + 9] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[5]); // e+a0+a5
				metrics_vec[2][c_num * path + 10] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[6]); // e+a0+a6
				metrics_vec[2][c_num * path + 11] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[7]); // e+a0+a7
				metrics_vec[2][c_num * path + 12] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[4]); // e+a1+a4
				metrics_vec[2][c_num * path + 13] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[5]); // e+a1+a5
				metrics_vec[2][c_num * path + 14] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[6]); // e+a1+a6
				metrics_vec[2][c_num * path + 15] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[7]); // e+a1+a7
				metrics_vec[2][c_num * path + 16] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[4]); // e+a2+a4
				metrics_vec[2][c_num * path + 17] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[5]); // e+a2+a5
				metrics_vec[2][c_num * path + 18] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[6]); // e+a2+a6
				metrics_vec[2][c_num * path + 19] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[4]); // e+a3+a4
				metrics_vec[2][c_num * path + 20] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[5]); // e+a3+a5
				metrics_vec[2][c_num * path + 21] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[4] + pen[5]); // e+a4+a5
				metrics_vec[2][c_num * path + 22] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[4]); // e+a0+a1+a2+a4
				metrics_vec[2][c_num * path + 23] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[5]); // e+a0+a1+a2+a5
				metrics_vec[2][c_num * path + 24] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[3] + pen[4]); // e+a0+a1+a3+a4
				metrics_vec[2][c_num * path + 25] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[3] + pen[4]); // e+a0+a2+a3+a4
				metrics_vec[2][c_num * path + 26] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[3] + pen[4]); // e+a1+a2+a3+a4
				if (n_elmts >= 16)
				{
					metrics_vec[2][c_num * path + 27] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[8]); // e+a0+a8
					metrics_vec[2][c_num * path + 28] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[9]); // e+a0+a9
					metrics_vec[2][c_num * path + 29] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[10]); // e+a0+a10
					metrics_vec[2][c_num * path + 30] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[11]); // e+a0+a11
					metrics_vec[2][c_num * path + 31] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[12]); // e+a0+a12
					metrics_vec[2][c_num * path + 32] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[13]); // e+a0+a13
					metrics_vec[2][c_num * path + 33] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[14]); // e+a0+a14
					metrics_vec[2][c_num * path + 34] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[15]); // e+a0+a15
					metrics_vec[2][c_num * path + 35] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[8]); // e+a1+a8
				}
			}
		}
		break;
	case 32:
		metrics_vec[2][c_num * path + 0] = metrics [path]; // empty
		metrics_vec[2][c_num * path + 1] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[1]); // e+a0+a1
		if (n_elmts >= 4)
		{
			metrics_vec[2][c_num * path + 2] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[2]); // e+a0+a2
			metrics_vec[2][c_num * path + 3] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[3]); // e+a0+a3
			metrics_vec[2][c_num * path + 4] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[2]); // e+a1+a2
			metrics_vec[2][c_num * path + 5] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[3]); // e+a1+a3
			metrics_vec[2][c_num * path + 6] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[3]); // e+a2+a3
			metrics_vec[2][c_num * path + 7] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[3]); // e+a0+a1+a2+a3
			if (n_elmts >= 8)
			{
				metrics_vec[2][c_num * path + 8] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[4]); // e+a0+a4
				metrics_vec[2][c_num * path + 9] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[5]); // e+a0+a5
				metrics_vec[2][c_num * path + 10] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[6]); // e+a0+a6
				metrics_vec[2][c_num * path + 11] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[7]); // e+a0+a7
				metrics_vec[2][c_num * path + 12] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[4]); // e+a1+a4
				metrics_vec[2][c_num * path + 13] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[5]); // e+a1+a5
				metrics_vec[2][c_num * path + 14] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[6]); // e+a1+a6
				metrics_vec[2][c_num * path + 15] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[7]); // e+a1+a7
				metrics_vec[2][c_num * path + 16] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[4]); // e+a2+a4
				metrics_vec[2][c_num * path + 17] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[5]); // e+a2+a5
				metrics_vec[2][c_num * path + 18] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[6]); // e+a2+a6
				metrics_vec[2][c_num * path + 19] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[7]); // e+a2+a7
				metrics_vec[2][c_num * path + 20] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[4]); // e+a3+a4
				metrics_vec[2][c_num * path + 21] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[5]); // e+a3+a5
				metrics_vec[2][c_num * path + 22] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[6]); // e+a3+a6
				metrics_vec[2][c_num * path + 23] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[7]); // e+a3+a7
				metrics_vec[2][c_num * path + 24] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[4] + pen[5]); // e+a4+a5
				metrics_vec[2][c_num * path + 25] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[4] + pen[6]); // e+a4+a6
				metrics_vec[2][c_num * path + 26] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[4] + pen[7]); // e+a4+a7
				metrics_vec[2][c_num * path + 27] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[5] + pen[6]); // e+a5+a6
				metrics_vec[2][c_num * path + 28] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[5] + pen[7]); // e+a5+a7
				metrics_vec[2][c_num * path + 29] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[6] + pen[7]); // e+a6+a7
				metrics_vec[2][c_num * path + 30] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[4]); // e+a0+a1+a2+a4
				metrics_vec[2][c_num * path + 31] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[5]); // e+a0+a1+a2+a5
				metrics_vec[2][c_num * path + 32] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[6]); // e+a0+a1+a2+a6
				metrics_vec[2][c_num * path + 33] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[7]); // e+a0+a1+a2+a7
				metrics_vec[2][c_num * path + 34] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[3] + pen[4]); // e+a0+a1+a3+a4
				metrics_vec[2][c_num * path + 35] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[3] + pen[5]); // e+a0+a1+a3+a5
				metrics_vec[2][c_num * path + 36] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[3] + pen[6]); // e+a0+a1+a3+a6
				metrics_vec[2][c_num * path + 37] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[3] + pen[7]); // e+a0+a1+a3+a7
				metrics_vec[2][c_num * path + 38] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[4] + pen[5]); // e+a0+a1+a4+a5
				metrics_vec[2][c_num * path + 39] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[4] + pen[6]); // e+a0+a1+a4+a6
				metrics_vec[2][c_num * path + 40] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[5] + pen[6]); // e+a0+a1+a5+a6
				metrics_vec[2][c_num * path + 41] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[3] + pen[4]); // e+a0+a2+a3+a4
				metrics_vec[2][c_num * path + 42] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[3] + pen[5]); // e+a0+a2+a3+a5
				metrics_vec[2][c_num * path + 43] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[3] + pen[6]); // e+a0+a2+a3+a6
				metrics_vec[2][c_num * path + 44] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[4] + pen[5]); // e+a0+a2+a4+a5
				metrics_vec[2][c_num * path + 45] = sat_m<R>(metrics_vec[2][c_num * path + 3] + pen[4] + pen[5]); // e+a0+a3+a4+a5
				metrics_vec[2][c_num * path + 46] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[3] + pen[4]); // e+a1+a2+a3+a4
				metrics_vec[2][c_num * path + 47] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[3] + pen[5]); // e+a1+a2+a3+a5
				metrics_vec[2][c_num * path + 48] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[3] + pen[6]); // e+a1+a2+a3+a6
				metrics_vec[2][c_num * path + 49] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[4] + pen[5]); // e+a1+a2+a4+a5
				metrics_vec[2][c_num * path + 50] = sat_m<R>(metrics_vec[2][c_num * path + 5] + pen[4] + pen[5]); // e+a1+a3+a4+a5
				metrics_vec[2][c_num * path + 51] = sat_m<R>(metrics_vec[2][c_num * path + 6] + pen[4] + pen[5]); // e+a2+a3+a4+a5
				metrics_vec[2][c_num * path + 52] = sat_m<R>(metrics_vec[2][c_num * path + 7] + pen[4] + pen[5]); // e+a0+a1+a2+a3+a4+a5
				if (n_elmts >= 16)
				{
					metrics_vec[2][c_num * path + 53] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[8]); // e+a0+a8
					metrics_vec[2][c_num * path + 54] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[9]); // e+a0+a9
					metrics_vec[2][c_num * path + 55] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[10]); // e+a0+a10
					metrics_vec[2][c_num * path + 56] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[11]); // e+a0+a11
					metrics_vec[2][c_num * path + 57] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[12]); // e+a0+a12
					metrics_vec[2][c_num * path + 58] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[13]); // e+a0+a13
					metrics_vec[2][c_num * path + 59] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[14]); // e+a0+a14
					metrics_vec[2][c_num * path + 60] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[15]); // e+a0+a15
					metrics_vec[2][c_num * path + 61] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[8]); // e+a1+a8
					metrics_vec[2][c_num * path + 62] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[9]); // e+a1+a9
					metrics_vec[2][c_num * path + 63] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[10]); // e+a1+a10
					metrics_vec[2][c_num * path + 64] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[11]); // e+a1+a11
					metrics_vec[2][c_num * path + 65] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[12]); // e+a1+a12
					metrics_vec[2][c_num * path + 66] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[13]); // e+a1+a13
					metrics_vec[2][c_num * path + 67] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[14]); // e+a1+a14
					metrics_vec[2][c_num * path + 68] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[15]); // e+a1+a15
					metrics_vec[2][c_num * path + 69] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[8]); // e+a2+a8
					metrics_vec[2][c_num * path + 70] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[9]); // e+a2+a9
					metrics_vec[2][c_num * path + 71] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[10]); // e+a2+a10
					metrics_vec[2][c_num * path + 72] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[11]); // e+a2+a11
					metrics_vec[2][c_num * path + 73] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[8]); // e+a3+a8
					metrics_vec[2][c_num * path + 74] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[9]); // e+a3+a9
					metrics_vec[2][c_num * path + 75] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[4] + pen[8]); // e+a4+a8
					metrics_vec[2][c_num * path + 76] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[8]); // e+a0+a1+a2+a8
					metrics_vec[2][c_num * path + 77] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[9]); // e+a0+a1+a2+a9
					if (n_elmts >= 32)
					{
						metrics_vec[2][c_num * path + 78] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[16]); // e+a0+a16
						metrics_vec[2][c_num * path + 79] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[17]); // e+a0+a17
						metrics_vec[2][c_num * path + 80] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[18]); // e+a0+a18
						metrics_vec[2][c_num * path + 81] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[19]); // e+a0+a19
						metrics_vec[2][c_num * path + 82] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[20]); // e+a0+a20
						metrics_vec[2][c_num * path + 83] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[21]); // e+a0+a21
						metrics_vec[2][c_num * path + 84] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[22]); // e+a0+a22
						metrics_vec[2][c_num * path + 85] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[23]); // e+a0+a23
						metrics_vec[2][c_num * path + 86] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[24]); // e+a0+a24
						metrics_vec[2][c_num * path + 87] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[25]); // e+a0+a25
						metrics_vec[2][c_num * path + 88] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[26]); // e+a0+a26
						metrics_vec[2][c_num * path + 89] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[27]); // e+a0+a27
						metrics_vec[2][c_num * path + 90] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[28]); // e+a0+a28
						metrics_vec[2][c_num * path + 91] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[29]); // e+a0+a29
						metrics_vec[2][c_num * path + 92] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[30]); // e+a0+a30
						metrics_vec[2][c_num * path + 93] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[31]); // e+a0+a31
						metrics_vec[2][c_num * path + 94] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[16]); // e+a1+a16
					}
				}
			}
		}
		break;
	case 64:
		metrics_vec[2][c_num * path + 0] = metrics [path]; // empty
		metrics_vec[2][c_num * path + 1] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[1]); // e+a0+a1
		if (n_elmts >= 4)
		{
			metrics_vec[2][c_num * path + 2] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[2]); // e+a0+a2
			metrics_vec[2][c_num * path + 3] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[3]); // e+a0+a3
			metrics_vec[2][c_num * path + 4] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[2]); // e+a1+a2
			metrics_vec[2][c_num * path + 5] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[3]); // e+a1+a3
			metrics_vec[2][c_num * path + 6] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[3]); // e+a2+a3
			metrics_vec[2][c_num * path + 7] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[3]); // e+a0+a1+a2+a3
			if (n_elmts >= 8)
			{
				metrics_vec[2][c_num * path + 8] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[4]); // e+a0+a4
				metrics_vec[2][c_num * path + 9] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[5]); // e+a0+a5
				metrics_vec[2][c_num * path + 10] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[6]); // e+a0+a6
				metrics_vec[2][c_num * path + 11] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[7]); // e+a0+a7
				metrics_vec[2][c_num * path + 12] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[4]); // e+a1+a4
				metrics_vec[2][c_num * path + 13] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[5]); // e+a1+a5
				metrics_vec[2][c_num * path + 14] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[6]); // e+a1+a6
				metrics_vec[2][c_num * path + 15] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[7]); // e+a1+a7
				metrics_vec[2][c_num * path + 16] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[4]); // e+a2+a4
				metrics_vec[2][c_num * path + 17] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[5]); // e+a2+a5
				metrics_vec[2][c_num * path + 18] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[6]); // e+a2+a6
				metrics_vec[2][c_num * path + 19] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[7]); // e+a2+a7
				metrics_vec[2][c_num * path + 20] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[4]); // e+a3+a4
				metrics_vec[2][c_num * path + 21] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[5]); // e+a3+a5
				metrics_vec[2][c_num * path + 22] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[6]); // e+a3+a6
				metrics_vec[2][c_num * path + 23] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[7]); // e+a3+a7
				metrics_vec[2][c_num * path + 24] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[4] + pen[5]); // e+a4+a5
				metrics_vec[2][c_num * path + 25] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[4] + pen[6]); // e+a4+a6
				metrics_vec[2][c_num * path + 26] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[4] + pen[7]); // e+a4+a7
				metrics_vec[2][c_num * path + 27] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[5] + pen[6]); // e+a5+a6
				metrics_vec[2][c_num * path + 28] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[5] + pen[7]); // e+a5+a7
				metrics_vec[2][c_num * path + 29] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[6] + pen[7]); // e+a6+a7
				metrics_vec[2][c_num * path + 30] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[4]); // e+a0+a1+a2+a4
				metrics_vec[2][c_num * path + 31] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[5]); // e+a0+a1+a2+a5
				metrics_vec[2][c_num * path + 32] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[6]); // e+a0+a1+a2+a6
				metrics_vec[2][c_num * path + 33] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[7]); // e+a0+a1+a2+a7
				metrics_vec[2][c_num * path + 34] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[3] + pen[4]); // e+a0+a1+a3+a4
				metrics_vec[2][c_num * path + 35] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[3] + pen[5]); // e+a0+a1+a3+a5
				metrics_vec[2][c_num * path + 36] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[3] + pen[6]); // e+a0+a1+a3+a6
				metrics_vec[2][c_num * path + 37] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[3] + pen[7]); // e+a0+a1+a3+a7
				metrics_vec[2][c_num * path + 38] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[4] + pen[5]); // e+a0+a1+a4+a5
				metrics_vec[2][c_num * path + 39] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[4] + pen[6]); // e+a0+a1+a4+a6
				metrics_vec[2][c_num * path + 40] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[4] + pen[7]); // e+a0+a1+a4+a7
				metrics_vec[2][c_num * path + 41] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[5] + pen[6]); // e+a0+a1+a5+a6
				metrics_vec[2][c_num * path + 42] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[5] + pen[7]); // e+a0+a1+a5+a7
				metrics_vec[2][c_num * path + 43] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[6] + pen[7]); // e+a0+a1+a6+a7
				metrics_vec[2][c_num * path + 44] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[3] + pen[4]); // e+a0+a2+a3+a4
				metrics_vec[2][c_num * path + 45] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[3] + pen[5]); // e+a0+a2+a3+a5
				metrics_vec[2][c_num * path + 46] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[3] + pen[6]); // e+a0+a2+a3+a6
				metrics_vec[2][c_num * path + 47] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[3] + pen[7]); // e+a0+a2+a3+a7
				metrics_vec[2][c_num * path + 48] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[4] + pen[5]); // e+a0+a2+a4+a5
				metrics_vec[2][c_num * path + 49] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[4] + pen[6]); // e+a0+a2+a4+a6
				metrics_vec[2][c_num * path + 50] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[4] + pen[7]); // e+a0+a2+a4+a7
				metrics_vec[2][c_num * path + 51] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[5] + pen[6]); // e+a0+a2+a5+a6
				metrics_vec[2][c_num * path + 52] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[5] + pen[7]); // e+a0+a2+a5+a7
				metrics_vec[2][c_num * path + 53] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[6] + pen[7]); // e+a0+a2+a6+a7
				metrics_vec[2][c_num * path + 54] = sat_m<R>(metrics_vec[2][c_num * path + 3] + pen[4] + pen[5]); // e+a0+a3+a4+a5
				metrics_vec[2][c_num * path + 55] = sat_m<R>(metrics_vec[2][c_num * path + 3] + pen[4] + pen[6]); // e+a0+a3+a4+a6
				metrics_vec[2][c_num * path + 56] = sat_m<R>(metrics_vec[2][c_num * path + 3] + pen[4] + pen[7]); // e+a0+a3+a4+a7
				metrics_vec[2][c_num * path + 57] = sat_m<R>(metrics_vec[2][c_num * path + 3] + pen[5] + pen[6]); // e+a0+a3+a5+a6
				metrics_vec[2][c_num * path + 58] = sat_m<R>(metrics_vec[2][c_num * path + 3] + pen[5] + pen[7]); // e+a0+a3+a5+a7
				metrics_vec[2][c_num * path + 59] = sat_m<R>(metrics_vec[2][c_num * path + 3] + pen[6] + pen[7]); // e+a0+a3+a6+a7
				metrics_vec[2][c_num * path + 60] = sat_m<R>(metrics_vec[2][c_num * path + 8] + pen[5] + pen[6]); // e+a0+a4+a5+a6
				metrics_vec[2][c_num * path + 61] = sat_m<R>(metrics_vec[2][c_num * path + 8] + pen[5] + pen[7]); // e+a0+a4+a5+a7
				metrics_vec[2][c_num * path + 62] = sat_m<R>(metrics_vec[2][c_num * path + 8] + pen[6] + pen[7]); // e+a0+a4+a6+a7
				metrics_vec[2][c_num * path + 63] = sat_m<R>(metrics_vec[2][c_num * path + 9] + pen[6] + pen[7]); // e+a0+a5+a6+a7
				metrics_vec[2][c_num * path + 64] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[3] + pen[4]); // e+a1+a2+a3+a4
				metrics_vec[2][c_num * path + 65] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[3] + pen[5]); // e+a1+a2+a3+a5
				metrics_vec[2][c_num * path + 66] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[3] + pen[6]); // e+a1+a2+a3+a6
				metrics_vec[2][c_num * path + 67] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[3] + pen[7]); // e+a1+a2+a3+a7
				metrics_vec[2][c_num * path + 68] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[4] + pen[5]); // e+a1+a2+a4+a5
				metrics_vec[2][c_num * path + 69] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[4] + pen[6]); // e+a1+a2+a4+a6
				metrics_vec[2][c_num * path + 70] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[4] + pen[7]); // e+a1+a2+a4+a7
				metrics_vec[2][c_num * path + 71] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[5] + pen[6]); // e+a1+a2+a5+a6
				metrics_vec[2][c_num * path + 72] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[5] + pen[7]); // e+a1+a2+a5+a7
				metrics_vec[2][c_num * path + 73] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[6] + pen[7]); // e+a1+a2+a6+a7
				metrics_vec[2][c_num * path + 74] = sat_m<R>(metrics_vec[2][c_num * path + 5] + pen[4] + pen[5]); // e+a1+a3+a4+a5
				metrics_vec[2][c_num * path + 75] = sat_m<R>(metrics_vec[2][c_num * path + 5] + pen[4] + pen[6]); // e+a1+a3+a4+a6
				metrics_vec[2][c_num * path + 76] = sat_m<R>(metrics_vec[2][c_num * path + 5] + pen[4] + pen[7]); // e+a1+a3+a4+a7
				metrics_vec[2][c_num * path + 77] = sat_m<R>(metrics_vec[2][c_num * path + 5] + pen[5] + pen[6]); // e+a1+a3+a5+a6
				metrics_vec[2][c_num * path + 78] = sat_m<R>(metrics_vec[2][c_num * path + 12] + pen[5] + pen[6]); // e+a1+a4+a5+a6
				metrics_vec[2][c_num * path + 79] = sat_m<R>(metrics_vec[2][c_num * path + 6] + pen[4] + pen[5]); // e+a2+a3+a4+a5
				metrics_vec[2][c_num * path + 80] = sat_m<R>(metrics_vec[2][c_num * path + 6] + pen[4] + pen[6]); // e+a2+a3+a4+a6
				metrics_vec[2][c_num * path + 81] = sat_m<R>(metrics_vec[2][c_num * path + 6] + pen[4] + pen[7]); // e+a2+a3+a4+a7
				metrics_vec[2][c_num * path + 82] = sat_m<R>(metrics_vec[2][c_num * path + 6] + pen[5] + pen[6]); // e+a2+a3+a5+a6
				metrics_vec[2][c_num * path + 83] = sat_m<R>(metrics_vec[2][c_num * path + 16] + pen[5] + pen[6]); // e+a2+a4+a5+a6
				metrics_vec[2][c_num * path + 84] = sat_m<R>(metrics_vec[2][c_num * path + 20] + pen[5] + pen[6]); // e+a3+a4+a5+a6
				metrics_vec[2][c_num * path + 85] = sat_m<R>(metrics_vec[2][c_num * path + 7] + pen[4] + pen[5]); // e+a0+a1+a2+a3+a4+a5
				metrics_vec[2][c_num * path + 86] = sat_m<R>(metrics_vec[2][c_num * path + 7] + pen[4] + pen[6]); // e+a0+a1+a2+a3+a4+a6
				metrics_vec[2][c_num * path + 87] = sat_m<R>(metrics_vec[2][c_num * path + 7] + pen[4] + pen[7]); // e+a0+a1+a2+a3+a4+a7
				metrics_vec[2][c_num * path + 88] = sat_m<R>(metrics_vec[2][c_num * path + 7] + pen[5] + pen[6]); // e+a0+a1+a2+a3+a5+a6
				metrics_vec[2][c_num * path + 89] = sat_m<R>(metrics_vec[2][c_num * path + 30] + pen[5] + pen[6]); // e+a0+a1+a2+a4+a5+a6
				metrics_vec[2][c_num * path + 90] = sat_m<R>(metrics_vec[2][c_num * path + 34] + pen[5] + pen[6]); // e+a0+a1+a3+a4+a5+a6
				metrics_vec[2][c_num * path + 91] = sat_m<R>(metrics_vec[2][c_num * path + 44] + pen[5] + pen[6]); // e+a0+a2+a3+a4+a5+a6
				metrics_vec[2][c_num * path + 92] = sat_m<R>(metrics_vec[2][c_num * path + 64] + pen[5] + pen[6]); // e+a1+a2+a3+a4+a5+a6
				if (n_elmts >= 16)
				{
					metrics_vec[2][c_num * path + 93] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[8]); // e+a0+a8
					metrics_vec[2][c_num * path + 94] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[9]); // e+a0+a9
					metrics_vec[2][c_num * path + 95] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[10]); // e+a0+a10
					metrics_vec[2][c_num * path + 96] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[11]); // e+a0+a11
					metrics_vec[2][c_num * path + 97] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[12]); // e+a0+a12
					metrics_vec[2][c_num * path + 98] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[13]); // e+a0+a13
					metrics_vec[2][c_num * path + 99] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[14]); // e+a0+a14
					metrics_vec[2][c_num * path + 100] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[15]); // e+a0+a15
					metrics_vec[2][c_num * path + 101] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[8]); // e+a1+a8
					metrics_vec[2][c_num * path + 102] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[9]); // e+a1+a9
					metrics_vec[2][c_num * path + 103] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[10]); // e+a1+a10
					metrics_vec[2][c_num * path + 104] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[11]); // e+a1+a11
					metrics_vec[2][c_num * path + 105] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[12]); // e+a1+a12
					metrics_vec[2][c_num * path + 106] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[13]); // e+a1+a13
					metrics_vec[2][c_num * path + 107] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[14]); // e+a1+a14
					metrics_vec[2][c_num * path + 108] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[15]); // e+a1+a15
					metrics_vec[2][c_num * path + 109] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[8]); // e+a2+a8
					metrics_vec[2][c_num * path + 110] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[9]); // e+a2+a9
					metrics_vec[2][c_num * path + 111] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[10]); // e+a2+a10
					metrics_vec[2][c_num * path + 112] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[11]); // e+a2+a11
					metrics_vec[2][c_num * path + 113] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[12]); // e+a2+a12
					metrics_vec[2][c_num * path + 114] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[13]); // e+a2+a13
					metrics_vec[2][c_num * path + 115] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[14]); // e+a2+a14
					metrics_vec[2][c_num * path + 116] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[15]); // e+a2+a15
					metrics_vec[2][c_num * path + 117] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[8]); // e+a3+a8
					metrics_vec[2][c_num * path + 118] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[9]); // e+a3+a9
					metrics_vec[2][c_num * path + 119] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[10]); // e+a3+a10
					metrics_vec[2][c_num * path + 120] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[11]); // e+a3+a11
					metrics_vec[2][c_num * path + 121] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[12]); // e+a3+a12
					metrics_vec[2][c_num * path + 122] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[13]); // e+a3+a13
					metrics_vec[2][c_num * path + 123] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[14]); // e+a3+a14
					metrics_vec[2][c_num * path + 124] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[15]); // e+a3+a15
					metrics_vec[2][c_num * path + 125] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[4] + pen[8]); // e+a4+a8
					metrics_vec[2][c_num * path + 126] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[4] + pen[9]); // e+a4+a9
					metrics_vec[2][c_num * path + 127] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[4] + pen[10]); // e+a4+a10
					metrics_vec[2][c_num * path + 128] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[4] + pen[11]); // e+a4+a11
					metrics_vec[2][c_num * path + 129] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[4] + pen[12]); // e+a4+a12
					metrics_vec[2][c_num * path + 130] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[4] + pen[13]); // e+a4+a13
					metrics_vec[2][c_num * path + 131] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[4] + pen[14]); // e+a4+a14
					metrics_vec[2][c_num * path + 132] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[5] + pen[8]); // e+a5+a8
					metrics_vec[2][c_num * path + 133] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[5] + pen[9]); // e+a5+a9
					metrics_vec[2][c_num * path + 134] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[5] + pen[10]); // e+a5+a10
					metrics_vec[2][c_num * path + 135] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[5] + pen[11]); // e+a5+a11
					metrics_vec[2][c_num * path + 136] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[5] + pen[12]); // e+a5+a12
					metrics_vec[2][c_num * path + 137] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[5] + pen[13]); // e+a5+a13
					metrics_vec[2][c_num * path + 138] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[6] + pen[8]); // e+a6+a8
					metrics_vec[2][c_num * path + 139] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[6] + pen[9]); // e+a6+a9
					metrics_vec[2][c_num * path + 140] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[6] + pen[10]); // e+a6+a10
					metrics_vec[2][c_num * path + 141] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[6] + pen[11]); // e+a6+a11
					metrics_vec[2][c_num * path + 142] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[6] + pen[12]); // e+a6+a12
					metrics_vec[2][c_num * path + 143] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[7] + pen[8]); // e+a7+a8
					metrics_vec[2][c_num * path + 144] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[7] + pen[9]); // e+a7+a9
					metrics_vec[2][c_num * path + 145] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[7] + pen[10]); // e+a7+a10
					metrics_vec[2][c_num * path + 146] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[7] + pen[11]); // e+a7+a11
					metrics_vec[2][c_num * path + 147] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[8] + pen[9]); // e+a8+a9
					metrics_vec[2][c_num * path + 148] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[8] + pen[10]); // e+a8+a10
					metrics_vec[2][c_num * path + 149] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[8] + pen[11]); // e+a8+a11
					metrics_vec[2][c_num * path + 150] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[9] + pen[10]); // e+a9+a10
					metrics_vec[2][c_num * path + 151] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[8]); // e+a0+a1+a2+a8
					metrics_vec[2][c_num * path + 152] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[9]); // e+a0+a1+a2+a9
					metrics_vec[2][c_num * path + 153] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[10]); // e+a0+a1+a2+a10
					metrics_vec[2][c_num * path + 154] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[11]); // e+a0+a1+a2+a11
					metrics_vec[2][c_num * path + 155] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[12]); // e+a0+a1+a2+a12
					metrics_vec[2][c_num * path + 156] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[13]); // e+a0+a1+a2+a13
					metrics_vec[2][c_num * path + 157] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[14]); // e+a0+a1+a2+a14
					metrics_vec[2][c_num * path + 158] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[15]); // e+a0+a1+a2+a15
					metrics_vec[2][c_num * path + 159] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[3] + pen[8]); // e+a0+a1+a3+a8
					metrics_vec[2][c_num * path + 160] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[3] + pen[9]); // e+a0+a1+a3+a9
					metrics_vec[2][c_num * path + 161] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[3] + pen[10]); // e+a0+a1+a3+a10
					metrics_vec[2][c_num * path + 162] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[3] + pen[11]); // e+a0+a1+a3+a11
					metrics_vec[2][c_num * path + 163] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[3] + pen[12]); // e+a0+a1+a3+a12
					metrics_vec[2][c_num * path + 164] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[4] + pen[8]); // e+a0+a1+a4+a8
					metrics_vec[2][c_num * path + 165] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[4] + pen[9]); // e+a0+a1+a4+a9
					metrics_vec[2][c_num * path + 166] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[4] + pen[10]); // e+a0+a1+a4+a10
					metrics_vec[2][c_num * path + 167] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[5] + pen[8]); // e+a0+a1+a5+a8
					metrics_vec[2][c_num * path + 168] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[5] + pen[9]); // e+a0+a1+a5+a9
					metrics_vec[2][c_num * path + 169] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[6] + pen[8]); // e+a0+a1+a6+a8
					metrics_vec[2][c_num * path + 170] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[7] + pen[8]); // e+a0+a1+a7+a8
					metrics_vec[2][c_num * path + 171] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[3] + pen[8]); // e+a0+a2+a3+a8
					metrics_vec[2][c_num * path + 172] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[3] + pen[9]); // e+a0+a2+a3+a9
					metrics_vec[2][c_num * path + 173] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[3] + pen[10]); // e+a0+a2+a3+a10
					metrics_vec[2][c_num * path + 174] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[3] + pen[11]); // e+a0+a2+a3+a11
					metrics_vec[2][c_num * path + 175] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[4] + pen[8]); // e+a0+a2+a4+a8
					metrics_vec[2][c_num * path + 176] = sat_m<R>(metrics_vec[2][c_num * path + 2] + pen[5] + pen[8]); // e+a0+a2+a5+a8
					metrics_vec[2][c_num * path + 177] = sat_m<R>(metrics_vec[2][c_num * path + 3] + pen[4] + pen[8]); // e+a0+a3+a4+a8
					metrics_vec[2][c_num * path + 178] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[3] + pen[8]); // e+a1+a2+a3+a8
					metrics_vec[2][c_num * path + 179] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[3] + pen[9]); // e+a1+a2+a3+a9
					metrics_vec[2][c_num * path + 180] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[3] + pen[10]); // e+a1+a2+a3+a10
					metrics_vec[2][c_num * path + 181] = sat_m<R>(metrics_vec[2][c_num * path + 4] + pen[4] + pen[8]); // e+a1+a2+a4+a8
					if (n_elmts >= 32)
					{
						metrics_vec[2][c_num * path + 182] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[16]); // e+a0+a16
						metrics_vec[2][c_num * path + 183] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[17]); // e+a0+a17
						metrics_vec[2][c_num * path + 184] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[18]); // e+a0+a18
						metrics_vec[2][c_num * path + 185] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[19]); // e+a0+a19
						metrics_vec[2][c_num * path + 186] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[20]); // e+a0+a20
						metrics_vec[2][c_num * path + 187] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[21]); // e+a0+a21
						metrics_vec[2][c_num * path + 188] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[22]); // e+a0+a22
						metrics_vec[2][c_num * path + 189] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[23]); // e+a0+a23
						metrics_vec[2][c_num * path + 190] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[24]); // e+a0+a24
						metrics_vec[2][c_num * path + 191] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[25]); // e+a0+a25
						metrics_vec[2][c_num * path + 192] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[26]); // e+a0+a26
						metrics_vec[2][c_num * path + 193] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[27]); // e+a0+a27
						metrics_vec[2][c_num * path + 194] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[28]); // e+a0+a28
						metrics_vec[2][c_num * path + 195] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[29]); // e+a0+a29
						metrics_vec[2][c_num * path + 196] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[30]); // e+a0+a30
						metrics_vec[2][c_num * path + 197] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[31]); // e+a0+a31
						metrics_vec[2][c_num * path + 198] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[16]); // e+a1+a16
						metrics_vec[2][c_num * path + 199] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[17]); // e+a1+a17
						metrics_vec[2][c_num * path + 200] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[18]); // e+a1+a18
						metrics_vec[2][c_num * path + 201] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[19]); // e+a1+a19
						metrics_vec[2][c_num * path + 202] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[20]); // e+a1+a20
						metrics_vec[2][c_num * path + 203] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[21]); // e+a1+a21
						metrics_vec[2][c_num * path + 204] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[22]); // e+a1+a22
						metrics_vec[2][c_num * path + 205] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[23]); // e+a1+a23
						metrics_vec[2][c_num * path + 206] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[24]); // e+a1+a24
						metrics_vec[2][c_num * path + 207] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[25]); // e+a1+a25
						metrics_vec[2][c_num * path + 208] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[26]); // e+a1+a26
						metrics_vec[2][c_num * path + 209] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[27]); // e+a1+a27
						metrics_vec[2][c_num * path + 210] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[28]); // e+a1+a28
						metrics_vec[2][c_num * path + 211] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[29]); // e+a1+a29
						metrics_vec[2][c_num * path + 212] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[30]); // e+a1+a30
						metrics_vec[2][c_num * path + 213] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[31]); // e+a1+a31
						metrics_vec[2][c_num * path + 214] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[16]); // e+a2+a16
						metrics_vec[2][c_num * path + 215] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[17]); // e+a2+a17
						metrics_vec[2][c_num * path + 216] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[18]); // e+a2+a18
						metrics_vec[2][c_num * path + 217] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[19]); // e+a2+a19
						metrics_vec[2][c_num * path + 218] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[20]); // e+a2+a20
						metrics_vec[2][c_num * path + 219] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[21]); // e+a2+a21
						metrics_vec[2][c_num * path + 220] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[2] + pen[22]); // e+a2+a22
						metrics_vec[2][c_num * path + 221] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[16]); // e+a3+a16
						metrics_vec[2][c_num * path + 222] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[3] + pen[17]); // e+a3+a17
						metrics_vec[2][c_num * path + 223] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[16]); // e+a0+a1+a2+a16
						metrics_vec[2][c_num * path + 224] = sat_m<R>(metrics_vec[2][c_num * path + 1] + pen[2] + pen[17]); // e+a0+a1+a2+a17
						if (n_elmts >= 64)
						{
							metrics_vec[2][c_num * path + 225] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[32]); // e+a0+a32
							metrics_vec[2][c_num * path + 226] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[33]); // e+a0+a33
							metrics_vec[2][c_num * path + 227] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[34]); // e+a0+a34
							metrics_vec[2][c_num * path + 228] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[35]); // e+a0+a35
							metrics_vec[2][c_num * path + 229] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[36]); // e+a0+a36
							metrics_vec[2][c_num * path + 230] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[37]); // e+a0+a37
							metrics_vec[2][c_num * path + 231] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[38]); // e+a0+a38
							metrics_vec[2][c_num * path + 232] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[39]); // e+a0+a39
							metrics_vec[2][c_num * path + 233] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[40]); // e+a0+a40
							metrics_vec[2][c_num * path + 234] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[41]); // e+a0+a41
							metrics_vec[2][c_num * path + 235] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[42]); // e+a0+a42
							metrics_vec[2][c_num * path + 236] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[43]); // e+a0+a43
							metrics_vec[2][c_num * path + 237] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[44]); // e+a0+a44
							metrics_vec[2][c_num * path + 238] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[45]); // e+a0+a45
							metrics_vec[2][c_num * path + 239] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[46]); // e+a0+a46
							metrics_vec[2][c_num * path + 240] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[47]); // e+a0+a47
							metrics_vec[2][c_num * path + 241] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[48]); // e+a0+a48
							metrics_vec[2][c_num * path + 242] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[49]); // e+a0+a49
							metrics_vec[2][c_num * path + 243] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[50]); // e+a0+a50
							metrics_vec[2][c_num * path + 244] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[51]); // e+a0+a51
							metrics_vec[2][c_num * path + 245] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[52]); // e+a0+a52
							metrics_vec[2][c_num * path + 246] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[53]); // e+a0+a53
							metrics_vec[2][c_num * path + 247] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[54]); // e+a0+a54
							metrics_vec[2][c_num * path + 248] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[55]); // e+a0+a55
							metrics_vec[2][c_num * path + 249] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[56]); // e+a0+a56
							metrics_vec[2][c_num * path + 250] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[57]); // e+a0+a57
							metrics_vec[2][c_num * path + 251] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[58]); // e+a0+a58
							metrics_vec[2][c_num * path + 252] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[59]); // e+a0+a59
							metrics_vec[2][c_num * path + 253] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[60]); // e+a0+a60
							metrics_vec[2][c_num * path + 254] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[61]); // e+a0+a61
							metrics_vec[2][c_num * path + 255] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[62]); // e+a0+a62
							metrics_vec[2][c_num * path + 256] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[0] + pen[63]); // e+a0+a63
							metrics_vec[2][c_num * path + 257] = sat_m<R>(metrics_vec[2][c_num * path + 0] + pen[1] + pen[32]); // e+a1+a32
						}
					}
				}
			}
		}
		break;
	default:
	// Chase-II
		break;
	}
}
}
}
