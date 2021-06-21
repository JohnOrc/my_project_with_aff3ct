#include <aff3ct.hpp>

namespace aff3ct
{
namespace factory
{
extern const std::string Decoder_polar_name;
extern const std::string Decoder_polar_prefix;
class Decoder_polar_mc : public Decoder
{
public:
	// ----------------------------------------------------------------------------------------------------- PARAMETERS
	// optional parameters
	std::string simd_strategy = "";
	std::string polar_nodes   = "{R0,R0L,R1,REP,REPL,SPC}";
	bool        full_adaptive = true;
	int         n_ite         = 1;
	int         L             = 8;
	int         T             = 8;

	// -------------------------------------------------------------------------------------------------------- METHODS
	explicit Decoder_polar_mc(const std::string &p = Decoder_polar_prefix);
	virtual ~Decoder_polar_mc() = default;
	Decoder_polar_mc* clone() const;

	// parameters construction
	virtual void get_description(cli::Argument_map_info &args) const;
	virtual void store          (const cli::Argument_map_value &vals);
	virtual void get_headers    (std::map<std::string,tools::header_list>& headers, const bool full = true) const;

	// builder
	template <typename B = int, typename Q = float>
	module::Decoder_SISO<B,Q>* build_siso(const std::vector<bool> &frozen_bits,
	                                       module::Encoder<B> *encoder = nullptr) const;

	template <typename B = int, typename Q = float>
	module::Decoder_SIHO<B,Q>* build(const std::vector<bool> &frozen_bits, const module::CRC<B> *crc = nullptr,
	                                 module::Encoder<B> *encoder = nullptr) const;

	template <typename B = int, typename Q = float>
	module::Decoder_SIHO<B,Q>* build_gen(const module::CRC<B> *crc = nullptr,
	                                           module::Encoder<B> *encoder = nullptr) const;

	static const std::vector<bool>& get_frozen_bits(const std::string &implem);

private:
	template <typename B = int, typename Q = float, class API_polar>
	module::Decoder_SIHO<B,Q>* _build(const std::vector<bool> &frozen_bits, const module::CRC<B> *crc = nullptr,
	                                  module::Encoder<B> *encoder = nullptr) const;

	template <typename B = int, typename Q = float, class API_polar>
	module::Decoder_SIHO<B,Q>* _build_scl_fast(const std::vector<bool> &frozen_bits,
	                                           const module::CRC<B> *crc = nullptr,
	                                                 module::Encoder<B> *encoder = nullptr) const;

	template <typename B = int, typename Q = float, class API_polar>
	module::Decoder_SIHO<B,Q>* _build_gen(const module::CRC<B> *crc = nullptr,
	                                            module::Encoder<B> *encoder = nullptr) const;

protected:
	Decoder_polar_mc(const std::string &n, const std::string &p);
};
}
}
