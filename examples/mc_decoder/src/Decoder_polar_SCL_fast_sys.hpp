/*!
 * \file
 * \brief Class module::Decoder_polar_SCL_fast_sys.
 */
#ifndef DECODER_POLAR_SCL_FAST_SYS
#define DECODER_POLAR_SCL_FAST_SYS

#include <vector>
#include <mipp.h>
#include <aff3ct.hpp>


namespace aff3ct
{
namespace module
{
template <typename B = int, typename R = float,
		  class API_polar = tools::API_polar_dynamic_seq<B, R, tools::f_LLR <  R>,
		  													   tools::g_LLR <B,R>,
															   tools::g0_LLR<  R>,
															   tools::h_LLR <B,R>,
															   tools::xo_STD<B  >>>
class Decoder_polar_SCL_fast_sys : public Decoder_SIHO<B,R>, public tools::Interface_notify_frozenbits_update
{
protected:
	const 	int							m;				// graph depth
			int							L;				// maximum paths number
	const	std::vector<bool>&			frozen_bits;
			tools::Pattern_polar_parser	polar_patterns;

				std	::vector<int >		paths;			// active paths
				std ::vector<R	 >		metrics;		// path metrics
	std::vector<mipp::vector<R	>>		l;				// llrs
	std::vector<mipp::vector<B  >>		s;				// partial sums
	std::vector<std ::vector<R	>>		metrics_vec;	// list of candidate metrics to be sorted
				std ::vector<int >		dup_count;		// number of duplications of a path, at updating time
				std ::vector<int >		bit_flips;		// index of the bits to be flipped
				std ::vector<bool>		is_even;		// used to store parity of a spc node

	int									best_path;
	int									n_active_paths;

	// each following 2D vector is of size L * m
	std::vector<std::vector<int>>		n_array_ref;	// number of times an array is used
	std::vector<std::vector<int>>		path_2_array;	// give array used by a path

	tools::LC_sorter<R>					sorter;
	std::vector<int>					best_idx;
	mipp::vector<R>						l_tmp;

public:
	Decoder_polar_SCL_fast_sys(const int& K, const int& N, const int& L, const std::vector<bool>& frozen_bits);

	Decoder_polar_SCL_fast_sys(const int& K, const int& N, const int& L, const std::vector<bool>& frozen_bits,
							   const std::vector<tools::Pattern_polar_i*> &polar_patterns,
							   const int idx_r0, const int idx_r1);

	virtual ~Decoder_polar_SCL_fast_sys();

	virtual Decoder_polar_SCL_fast_sys<B,R,API_polar>* clone() const;

	virtual void notify_noise_update();

protected:
	virtual void _decode		(const R *Y_N								);
			int	 _decode_siho	(const R *Y_N, B *V_K, const size_t frame_id);
			int	 _decode_siho_cw(const R *Y_N, B *V_N, const size_t frame_id);
	virtual void _store			(			   B *V_K 						) const;
	virtual void _store_cw		(			   B *V_N 						) const;

	inline void recursive_decode(const R *Y_N, const int off_l, const int off_s, const int rev_depth, int &node_id);

	inline void update_paths_r0 (const int rev_depth, const int off_l, const off_s, const int n_elmts);
	inline void update_paths_r1 (const int rev_depth, const int off_l, const off_s, const int n_elmts);
	inline void update_paths_rep(const int rev_depth, const int off_l, const off_s, const int n_elmts);
	inline void update_paths_spc(const int rev_depth, const int off_l, const off_s, const int n_elmts);

	template <int REV_D, int N_ELMTS> inline void update_paths_r0 (const int off_l, const int off_s);
	template <int REV_D, int N_ELMTS> inline void update_paths_r1 (const int off_l, const int off_s);
	template <int REV_D, int N_ELMTS> inline void update_paths_rep(const int off_l, const int off_s);
	template <int REV_D, int N_ELMTS> inline void update_paths_spc(const int off_l, const int off_s);

	virtual inline void init_buffers	(							  );
			inline void delete_path		(int path_id				  );
	virtual inline int	select_best_path(const size_t frame_id		  );
			inline int	up_ref_array_idx(const int path, const int r_d);

private:
	inline void flip_bits_r1 (const int old_path, const int new_path, const int dup, const int off_s, const int n_elmts);
	inline void flip_bits_spc(const int old_path, const int new_path, const int dup, const int off_s, const int n_elmts);

	inline void erase_bad_paths	();
	inline int	duplicate_tree	(const int old_path, const int off_l, const int off_s, const int n_elmts);
};
}
}
