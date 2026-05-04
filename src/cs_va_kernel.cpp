#include "cs_va_kernel.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "core/constants.hpp"
#include "cs_ct_shift.hpp"
#include "cs_ioperator.hpp"
#include "virtual.hpp"

namespace ttbar_qtsub {
namespace {

double shat_from_momenta(const std::vector<FourVector>& momenta) {
    return project_sij(momenta[0], momenta[1], true, true);
}

double costheta_from_momenta(const std::vector<FourVector>& momenta, double mt) {
    const double shat = shat_from_momenta(momenta);
    if (shat <= 4.0 * mt * mt) {
        throw std::domain_error("costheta_from_momenta: shat must be above threshold");
    }

    const double beta = std::sqrt(1.0 - 4.0 * mt * mt / shat);
    if (beta <= 0.0) {
        throw std::domain_error("costheta_from_momenta: beta must be positive");
    }

    const double t = project_sij(momenta[0], momenta[2], true, false);
    double costheta = (-2.0 * mt * mt + shat + 2.0 * t) / (beta * shat);

    if (costheta > 1.0 && costheta < 1.0 + kEpsilonCosTheta) {
        costheta = 1.0;
    }
    else if (costheta < -1.0 && costheta > -1.0 - kEpsilonCosTheta) {
        costheta = -1.0;
    }
    else if (costheta < -1.0 || costheta > 1.0) {
        throw std::domain_error("costheta_from_momenta: reconstructed cos(theta) is outside [-1,1]");
    }

    return costheta;
}

double compute_local_virtual_finite(CSChannel channel,
                                    double shat,
                                    double mt,
                                    double muR,
                                    double alpha_s,
                                    double costheta) {
    const double gs = std::sqrt(4.0 * kPi * alpha_s);
    switch (channel) {
        case CSChannel::QQbar:
            return qqbttb1l(shat, mt, muR, gs, costheta);
        case CSChannel::GG:
            return ggttb1l(shat, mt, muR, gs, costheta);
    }

    throw std::invalid_argument("compute_local_virtual_finite: unsupported channel");
}

}  // namespace

VAKernelResult compute_va_kernel(
    CSChannel channel,
    const std::vector<FourVector>& momenta,
    double alpha_s,
    double muR,
    double mt) {
    if (momenta.size() != 4) {
        throw std::invalid_argument("compute_va_kernel: expected exactly four Born momenta");
    }

    const double shat = shat_from_momenta(momenta);
    const double costheta = costheta_from_momenta(momenta, mt);
    const double virtual_with_ct = compute_local_virtual_finite(
        channel, shat, mt, muR, alpha_s, costheta);
    const auto ct_shift_result = compute_ct_shift_local_exact(channel, momenta, alpha_s, muR, mt);
    const double ct_shift = ct_shift_result.total;
    const double virtual_finite = virtual_with_ct - ct_shift;
    const double ioperator_finite = assemble_cdst_ioperator_finite(
        channel,
        shat,
        mt,
        muR,
        alpha_s,
        costheta,
        momenta[0],
        momenta[1],
        momenta[2],
        momenta[3]).finite;
    return VAKernelResult{
        virtual_finite,
        ct_shift,
        ioperator_finite,
        virtual_finite + ioperator_finite,
        virtual_finite + ct_shift + ioperator_finite,
    };
}

}  // namespace ttbar_qtsub