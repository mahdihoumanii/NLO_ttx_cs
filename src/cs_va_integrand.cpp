#include "cs_va_integrand.hpp"

namespace ttbar_qtsub {

double va_integrand_factor_product(const VAIntegrandFactors& factors) {
    return factors.ps_factor
         * factors.ew_rescaling
         * factors.pdf_factor
         * factors.alpha_s_rescaling;
}

VAIntegrandResult assemble_va_integrand(
    double kernel,
    const VAIntegrandFactors& factors) {
    const double factor_product = va_integrand_factor_product(factors);
    return VAIntegrandResult{
        kernel,
        factors,
        factor_product,
        factor_product * kernel,
    };
}

VAIntegrandResult assemble_va_integrand(
    const VAKernelResult& kernel,
    const VAIntegrandFactors& factors) {
    return assemble_va_integrand(kernel.VA_total, factors);
}

VAIntegrandResult compute_va_integrand(
    CSChannel channel,
    const std::vector<FourVector>& momenta,
    double alpha_s,
    double muR,
    double mt,
    const VAIntegrandFactors& factors) {
    return assemble_va_integrand(
        compute_va_kernel(channel, momenta, alpha_s, muR, mt),
        factors);
}

}  // namespace ttbar_qtsub