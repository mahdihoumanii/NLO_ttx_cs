#include "cs_ct_shift.hpp"

#include <cmath>
#include <stdexcept>

#include "born.hpp"
#include "core/constants.hpp"

namespace ttbar_qtsub {
namespace {

constexpr int kBornAlphaSPower = 2;

double shat_from_momenta(const std::vector<FourVector>& momenta) {
    if (momenta.size() != 4) {
        throw std::invalid_argument("shat_from_momenta: expected exactly four Born momenta");
    }
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

double decoupling_value(double Born, double alpha_s, double muR, double mt) {
    const double log_mu_mt = std::log((muR * muR) / (mt * mt));
    return -static_cast<double>(kBornAlphaSPower)
         * (alpha_s / (3.0 * kPi))
         * kTR
         * log_mu_mt
         * Born;
}

double born_value_for_channel(CSChannel channel,
                              double shat,
                              double mt,
                              double muR,
                              double alpha_s,
                              double costheta) {
    const double gs = std::sqrt(4.0 * kPi * alpha_s);
    switch (channel) {
        case CSChannel::QQbar:
            return born_qq(shat, mt, muR, gs, costheta);
        case CSChannel::GG:
            return born_gg(shat, mt, muR, gs, costheta);
    }

    throw std::invalid_argument("born_value_for_channel: unsupported channel");
}

double gg_top_os_counterterm_shape(double beta, double costheta) {
    const double b2 = beta * beta;
    const double c2 = costheta * costheta;
    const double b4 = b2 * b2;
    const double b6 = b4 * b2;
    const double b8 = b4 * b4;
    const double c4 = c2 * c2;
    const double c6 = c4 * c2;
    const double c8 = c4 * c4;

    const double numerator =
        (-14.0 * b2 + 23.0 * b2 * c2
         - 7.0 * b4 + 17.0 * b4 * c2 - 37.0 * b4 * c4
         + 14.0 * b6 - 92.0 * b6 * c2 + 146.0 * b6 * c4 - 41.0 * b6 * c6
         + 50.0 * b8 * c2 - 93.0 * b8 * c4 + 43.0 * b8 * c6 - 9.0 * b8 * c8) / 9.0;
    const double denominator = std::pow(1.0 - b2 * c2, 3);
    return numerator / denominator;
}

double local_exact_total_value(CSChannel channel,
                               double shat,
                               double mt,
                               double muR,
                               double alpha_s,
                               double costheta) {
    if (muR <= 0.0) {
        throw std::domain_error("local_exact_total_value: muR must be positive");
    }
    if (mt <= 0.0) {
        throw std::domain_error("local_exact_total_value: mt must be positive");
    }
    if (shat <= 4.0 * mt * mt) {
        return 0.0;
    }

    const double log_mu_mt = std::log((muR * muR) / (mt * mt));
    const double alpha_over_2pi = alpha_s / (2.0 * kPi);
    const double gs2 = 4.0 * kPi * alpha_s;
    const double gs4 = gs2 * gs2;

    switch (channel) {
        case CSChannel::QQbar: {
            const double born_value = born_value_for_channel(channel, shat, mt, muR, alpha_s, costheta);
            const double coefficient = -4.0 * kCF - (10.0 / 3.0) * log_mu_mt;
            return alpha_over_2pi * coefficient * born_value;
        }
        case CSChannel::GG: {
            const double beta = std::sqrt(1.0 - 4.0 * mt * mt / shat);
            const double scale_factor = 1.0 + 0.75 * log_mu_mt;
            return alpha_over_2pi * gs4 * scale_factor * gg_top_os_counterterm_shape(beta, costheta);
        }
    }

    throw std::invalid_argument("local_exact_total_value: unsupported channel");
}

}  // namespace

CTShiftResult compute_ct_shift_decoupling_approx(
    CSChannel channel,
    double Born,
    double alpha_s,
    double muR,
    double mt) {
    switch (channel) {
        case CSChannel::QQbar:
        case CSChannel::GG:
            break;
    }

    if (muR <= 0.0) {
        throw std::domain_error("compute_ct_shift: muR must be positive");
    }
    if (mt <= 0.0) {
        throw std::domain_error("compute_ct_shift: mt must be positive");
    }

    const double decoupling = decoupling_value(Born, alpha_s, muR, mt);

    return CTShiftResult{decoupling, decoupling};
}

CTShiftResult compute_ct_shift_local_exact(
    CSChannel channel,
    const std::vector<FourVector>& momenta,
    double alpha_s,
    double muR,
    double mt) {
    if (muR <= 0.0) {
        throw std::domain_error("compute_ct_shift_local_exact: muR must be positive");
    }
    if (mt <= 0.0) {
        throw std::domain_error("compute_ct_shift_local_exact: mt must be positive");
    }

    const double shat = shat_from_momenta(momenta);
    const double costheta = costheta_from_momenta(momenta, mt);
    const double born_value = born_value_for_channel(channel, shat, mt, muR, alpha_s, costheta);
    const double exact = local_exact_total_value(channel, shat, mt, muR, alpha_s, costheta);
    const double decoupling = decoupling_value(born_value, alpha_s, muR, mt);
    return CTShiftResult{exact, decoupling};
}

CTShiftResult compute_ct_shift(
    CSChannel channel,
    double Born,
    double alpha_s,
    double muR,
    double mt) {
    return compute_ct_shift_decoupling_approx(channel, Born, alpha_s, muR, mt);
}

}  // namespace ttbar_qtsub