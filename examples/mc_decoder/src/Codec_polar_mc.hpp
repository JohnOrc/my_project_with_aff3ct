#include <vector>
#include <string>
#include <map>
#include <cli.hpp>

#include <aff3ct.hpp>

namespace aff3ct
{
namespace factory
{
extern const std::string Codec_polar_prefix;
class Codec_polar_mc : public Codec_polar
{
public:
	// -------------------------------------------------------------------------------------------------------- METHODS
	explicit Codec_polar_mc(const std::string &p = Codec_polar_prefix);
	virtual ~Codec_polar_mc() = default;
	Codec_polar_mc* clone() const;

	// TODO
	tools::Codec_polar<B,Q>* build(const module::CRC<B>* crc = nullptr) const;
};
}
}