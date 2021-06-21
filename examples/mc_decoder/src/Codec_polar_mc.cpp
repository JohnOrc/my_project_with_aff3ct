#include <aff3ct.hpp>

#include "Codec_polar_mc.hpp"
//#include "Decoder_polar_mc.hpp"

namespace aff3ct
{
namespace factory
{
Codec_polar_mc
::Codec_polar_mc(const std::string &prefix)
: Codec_polar(Codec_polar_name, prefix)
{
	Codec::set_enc(new Encoder_polar("enc"));
	Codec::set_dec(new Decoder_polar("dec"));

	std::cout << "COdec_polar_mc loaded" << std::endl;

  	fbg = new Frozenbits_generator(enc->get_prefix()+"-fb");
}

Codec_polar_mc* Codec_polar_mc
::clone() const
{
	return new Codec_polar_mc(*this);
}


template <typename B, typename Q>
tools::Codec_polar<B,Q>* Codec_polar_mc
::build(const module::CRC<B> *crc) const
{
	return new tools::Codec_polar<B,Q>(*fbg,
	                                   dynamic_cast<const Encoder_polar  &>(*enc),
	                                   dynamic_cast<const Decoder_polar  &>(*dec),
	                                   dynamic_cast<const Puncturer_polar*>(pct.get()),
	                                   crc);
}

// ==================================================================================== explicit template instantiation
#include "Tools/types.h"
#ifdef AFF3CT_MULTI_PREC
template aff3ct::tools::Codec_polar<B_8 ,Q_8 >* aff3ct::factory::Codec_polar_mc::build<B_8 ,Q_8 >(const aff3ct::module::CRC<B_8 >*) const;
template aff3ct::tools::Codec_polar<B_16,Q_16>* aff3ct::factory::Codec_polar_mc::build<B_16,Q_16>(const aff3ct::module::CRC<B_16>*) const;
template aff3ct::tools::Codec_polar<B_32,Q_32>* aff3ct::factory::Codec_polar_mc::build<B_32,Q_32>(const aff3ct::module::CRC<B_32>*) const;
template aff3ct::tools::Codec_polar<B_64,Q_64>* aff3ct::factory::Codec_polar_mc::build<B_64,Q_64>(const aff3ct::module::CRC<B_64>*) const;
#else
template aff3ct::tools::Codec_polar<B,Q>* aff3ct::factory::Codec_polar_mc::build<B,Q>(const aff3ct::module::CRC<B>*) const;
#endif
// ==================================================================================== explicit template instantiation
}
}