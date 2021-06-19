#include <vector>
#include <string>
#include <map>
#include <cli.hpp>

#include "Tools/Factory/Header.hpp"
#include "Module/CRC/CRC.hpp"
#include "Tools/Codec/Polar/Codec_polar.hpp"
#include "Factory/Tools/Code/Polar/Frozenbits_generator.hpp"
#include "Factory/Tools/Codec/Codec_SISO.hpp"

namespace aff3ct
{
namespace factory
{
extern const std::string Codec_polar_name;
extern const std::string Codec_polar_prefix;
class Codec_polar_mc : public Codec_SISO
{
public:
	// ----------------------------------------------------------------------------------------------------- PARAMETERS
	// depending parameters
	tools::auto_cloned_unique_ptr<Frozenbits_generator> fbg;

	// -------------------------------------------------------------------------------------------------------- METHODS
	explicit Codec_polar_mc(const std::string &p = Codec_polar_prefix);
	virtual ~Codec_polar_mc() = default;
	Codec_polar_mc* clone() const;
	void enable_puncturer();

	virtual std::vector<std::string> get_names      () const;
	virtual std::vector<std::string> get_short_names() const;
	virtual std::vector<std::string> get_prefixes   () const;

	// parameters construction
	void get_description(cli::Argument_map_info &args) const;
	void store          (const cli::Argument_map_value &vals);
	void get_headers    (std::map<std::string,tools::header_list>& headers, const bool full = true) const;

	// builder
	template <typename B = int, typename Q = float>
	tools::Codec_polar<B,Q>* build(const module::CRC<B>* crc = nullptr) const;
};
}
}