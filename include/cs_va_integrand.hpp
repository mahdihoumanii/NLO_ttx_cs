#pragma once

#include <vector>

#include "cs_born_correlations.hpp"
#include "cs_kinematics.hpp"
#include "cs_va_kernel.hpp"

namespace ttbar_qtsub {

struct VAIntegrandFactors {
    double ps_factor = 1.0;
    double pdf_factor = 1.0;
    double alpha_s_rescaling = 1.0;
    double ew_rescaling = 1.0;
};

struct VAIntegrandResult {
    double kernel;
    VAIntegrandFactors factors;
    double factor_product;
    double integrand;
};

// Multiply the four post-kernel factors in the same order used by
// determine_integrand_VA(): ps_factor * ew_rescaling * pdf_factor * alpha_s_rescaling.
double va_integrand_factor_product(const VAIntegrandFactors& factors);

// Assemble the final VA integrand from a prevalidated VA kernel value.
VAIntegrandResult assemble_va_integrand(
    double kernel,
    const VAIntegrandFactors& factors);

// Convenience overload when the validated kernel pieces are already available.
VAIntegrandResult assemble_va_integrand(
    const VAKernelResult& kernel,
    const VAIntegrandFactors& factors);

// Compute the validated VA kernel and multiply it by explicit post-kernel
// factors. This mirrors the validated weight structure without trying to
// reproduce an external phase-space generator.
VAIntegrandResult compute_va_integrand(
    CSChannel channel,
    const std::vector<FourVector>& momenta,
    double alpha_s,
    double muR,
    double mt,
    const VAIntegrandFactors& factors);

}  // namespace ttbar_qtsub