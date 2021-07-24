#include <vector>
#include <aff3ct.hpp>

/**
 * @brief Generate frozen bits for a specific SNR
 *
 * @param k The number of information bits
 * @param n The codeword length
 * @param snr_max Estimated SNR (in dB)
 * @return auto std::vector<bool> of frozen bits
 */
std::vector<bool> generate_frozen_bits(const int k, const int n,
                                       const float snr_max) {
  // calculate constants
  auto r = static_cast<float>(k * 1.0 / n);
  const auto esn0 = aff3ct::tools::ebn0_to_esn0(snr_max, r);
  const auto ebn0 = aff3ct::tools::esn0_to_ebn0(esn0);
  const auto sigma = aff3ct::tools::esn0_to_sigma(esn0);

  // set noise
  aff3ct::tools::Frozenbits_generator_GA_Arikan frozen_bits_generator(k, n);
  auto noise = std::unique_ptr<aff3ct::tools::Sigma<>>(
    new aff3ct::tools::Sigma<>());
  noise->set_values(sigma, ebn0, esn0);

  // generate frozen bits
  std::vector<bool> frozen_bits(n);
  frozen_bits_generator.set_noise(*noise);
  frozen_bits_generator.generate(frozen_bits);
  return frozen_bits;
}
