#include "cs_ioperator.hpp"

#include <cmath>
#include <stdexcept>

#include "core/constants.hpp"

namespace ttbar_qtsub {
namespace {

constexpr double kDeltaIR1 = 0.0;
constexpr double kDeltaIR2 = 0.0;
constexpr int kSwitchPolenorm = 0;

double dilog_real(double x) {
    if (std::abs(x) < 1.0e-15) {
        return 0.0;
    }
    if (std::abs(x - 1.0) < 1.0e-15) {
        return kPiSq / 6.0;
    }
    if (x < 0.0) {
        const double y = x / (x - 1.0);
        const double ly = std::log(1.0 - x);
        return -dilog_real(y) - 0.5 * ly * ly;
    }
    if (x > 0.5) {
        return kPiSq / 6.0 - std::log(x) * std::log(1.0 - x) - dilog_real(1.0 - x);
    }

    double sum = 0.0;
    double xk = x;
    for (int k = 1; k <= 500000; ++k) {
        const double denom = static_cast<double>(k) * static_cast<double>(k);
        const double term = xk / denom;
        sum += term;
        if (std::abs(term) < 1.0e-15) {
            break;
        }
        xk *= x;
    }
    return sum;
}

double kallen_lambda(double a, double b, double c) {
    return a * a + b * b + c * c - 2.0 * a * b - 2.0 * a * c - 2.0 * b * c;
}

double gamma_q() {
    return 1.5 * kCF;
}

double gamma_g() {
    return (11.0 / 6.0) * kCA - (2.0 / 3.0) * kTR * static_cast<double>(kNfDefault);
}

double K_q() {
    return (3.5 - kPiSq / 6.0) * kCF;
}

double K_g() {
    return (67.0 / 18.0 - kPiSq / 6.0) * kCA - (10.0 / 9.0) * kTR * static_cast<double>(kNfDefault);
}

int pair_mass_pattern(const OrderedColorCorrelation& pair) {
    if (!pair.emitter_massive && !pair.spectator_massive) {
        return 0;
    }
    if (pair.emitter_massive && !pair.spectator_massive) {
        return 1;
    }
    if (!pair.emitter_massive && pair.spectator_massive) {
        return 2;
    }
    return 3;
}

bool emitter_is_gluon(const OrderedColorCorrelation& pair) {
    return pair.emitter_type == "gluon";
}

double finite0_for_emitter(const OrderedColorCorrelation& pair, double mt, double sij) {
    const double temp_finite_massive_full_q = -kPiSq / 3.0 + (gamma_q() + K_q()) / kCF;
    const double temp_finite_massive_full_g_in = -kPiSq / 3.0 + (gamma_g() + K_g()) / kCA;

    if (pair.emitter_massive) {
        const double gamma_q0 = kCF * (kDeltaIR1 + 0.5 * std::log((mt * mt) / sij) - 2.0);
        return temp_finite_massive_full_q + gamma_q0 / kCF;
    }
    if (emitter_is_gluon(pair)) {
        return temp_finite_massive_full_g_in;
    }
    return -kPiSq / 3.0 + (gamma_q() + K_q()) / kCF;
}

double finite1_for_emitter(const OrderedColorCorrelation& pair) {
    if (pair.emitter_massive) {
        const double gamma_q1 = -0.5 * kCF;
        return (gamma_q1 + gamma_q()) / kCF;
    }
    if (emitter_is_gluon(pair)) {
        return gamma_g() / kCA;
    }
    return gamma_q() / kCF;
}

CDSTPairCoefficients compute_pair_coefficients(const OrderedColorCorrelation& pair,
                                               double mt,
                                               double mu_reg,
                                               double alpha_s) {
    if (pair.matrix_sij <= 0.0) {
        throw std::domain_error("assemble_cdst_ioperator_finite: matrix_sij must be positive");
    }

    const double mu2IR = mu_reg * mu_reg;
    const double sij = pair.matrix_sij;
    const double mt2 = mt * mt;

    CDSTPairCoefficients result{};
    result.emitter = pair.emitter;
    result.spectator = pair.spectator;
    result.Bij = pair.Bij;
    result.matrix_sij = sij;
    result.finite0 = finite0_for_emitter(pair, mt, sij);
    result.finite1 = finite1_for_emitter(pair);

    const int mass_pattern = pair_mass_pattern(pair);
    if (mass_pattern == 0) {
        result.VS0 = kDeltaIR2;
        if (kSwitchPolenorm == 1) {
            result.VS0 += -kPiSq / 6.0;
        }
        result.VS1 = kDeltaIR1;
        result.VS2 = 1.0;
        result.VNS = 0.0;
    }
    else if (mass_pattern == 1 || mass_pattern == 2) {
        const double temp_mass2 = mt2;
        const double Q2 = temp_mass2 + sij;
        const double sqrtQ2 = std::sqrt(Q2);
        const double temp_log_s_Q2 = std::log(sij / Q2);
        const double log_mass2_over_s = std::log(temp_mass2 / sij);
        const double log_mass2_over_Q2 = std::log(temp_mass2 / Q2);

        result.VS0 = 0.5 * kDeltaIR2
                   + 0.5 * log_mass2_over_s * kDeltaIR1
                   - 0.25 * log_mass2_over_s * log_mass2_over_s
                   - kPiSq / 12.0
                   - 0.5 * log_mass2_over_s * temp_log_s_Q2
                   - 0.5 * log_mass2_over_Q2 * temp_log_s_Q2;
        result.VS1 = 0.5 * kDeltaIR1 + 0.5 * log_mass2_over_s;
        result.VS2 = 0.5;
        if (kSwitchPolenorm == 1) {
            result.VS0 += -0.5 * (kPiSq / 6.0);
        }

        if (mass_pattern == 1) {
            result.VNS = gamma_q() / kCF * std::log(sij / Q2)
                       + kPiSq / 6.0
                       - dilog_real(sij / Q2)
                       - 2.0 * std::log(sij / Q2)
                       - temp_mass2 / sij * std::log(temp_mass2 / Q2);
        }
        else {
            if (emitter_is_gluon(pair)) {
                result.VNS = gamma_g() / kCA
                           * (std::log(sij / Q2)
                              - 2.0 * std::log((sqrtQ2 - mt) / sqrtQ2)
                              - 2.0 * mt / (sqrtQ2 + mt))
                           + kPiSq / 6.0
                           - dilog_real(sij / Q2);
                // The extra heavy-flavor correction is only activated
                // for final-state massless gluon emitters. Born ttbar pairs do
                // not contain such emitters, so no extra term is added here.
            }
            else {
                result.VNS = gamma_q() / kCF
                           * (std::log(sij / Q2)
                              - 2.0 * std::log((sqrtQ2 - mt) / sqrtQ2)
                              - 2.0 * mt / (sqrtQ2 + mt))
                           + kPiSq / 6.0
                           - dilog_real(sij / Q2);
            }
        }
    }
    else {
        const double emitter_mass2 = mt2;
        const double spectator_mass2 = mt2;
        const double Q2 = emitter_mass2 + spectator_mass2 + sij;
        const double sqrtQ2 = std::sqrt(Q2);
        const double mu2_emitter = emitter_mass2 / Q2;
        const double mu2_spectator = spectator_mass2 / Q2;
        const double temp_r2jk = emitter_mass2 * spectator_mass2 / (sij * sij);

        double delta = 0.0;
        double vik = 0.0;
        if (temp_r2jk < 1.0e-3) {
            delta = 2.0 * temp_r2jk
                  + 2.0 * std::pow(temp_r2jk, 2)
                  + 4.0 * std::pow(temp_r2jk, 3)
                  + 10.0 * std::pow(temp_r2jk, 4)
                  + 28.0 * std::pow(temp_r2jk, 5);
            vik = 1.0 - delta;
        }
        else {
            vik = std::sqrt(kallen_lambda(1.0, mu2_emitter, mu2_spectator))
                / (1.0 - mu2_emitter - mu2_spectator);
            delta = 1.0 - vik;
        }

        const double rho2 = delta / (1.0 + vik);
        const double logrho2 = std::log(rho2);
        const double logrho = 0.5 * logrho2;
        const double rho2_emitter = (1.0 - vik + 2.0 * mu2_emitter / (1.0 - mu2_emitter - mu2_spectator))
                                  / (1.0 + vik + 2.0 * mu2_emitter / (1.0 - mu2_emitter - mu2_spectator));
        const double rho2_spectator = (1.0 - vik + 2.0 * mu2_spectator / (1.0 - mu2_emitter - mu2_spectator))
                                    / (1.0 + vik + 2.0 * mu2_spectator / (1.0 - mu2_emitter - mu2_spectator));

        result.VS0 = (1.0 / vik)
                   * (logrho * kDeltaIR1
                      - 0.25 * std::pow(std::log(rho2_emitter), 2)
                      - 0.25 * std::pow(std::log(rho2_spectator), 2)
                      - kPiSq / 6.0
                      + logrho * std::log(Q2 / sij));
        result.VS1 = logrho / vik;
        result.VS2 = 0.0;

        result.VNS = gamma_q() / kCF * std::log(sij / Q2)
                   + (1.0 / vik)
                     * (logrho2 * std::log(1.0 + rho2)
                        + 2.0 * dilog_real(rho2)
                        - dilog_real(1.0 - rho2_emitter)
                        - dilog_real(1.0 - rho2_spectator)
                        - kPiSq / 6.0)
                   + std::log((sqrtQ2 - mt) / sqrtQ2)
                   - 2.0 * std::log((std::pow(sqrtQ2 - mt, 2) - mt2) / Q2)
                   - 2.0 * mt2 / sij * std::log(mt / (sqrtQ2 - mt))
                   - mt / (sqrtQ2 - mt)
                   + 2.0 * mt * (2.0 * mt - sqrtQ2) / sij
                   + kPiSq / 2.0;
    }

    result.log_mu2_over_sij = std::log(mu2IR / sij);
    result.log2_half = 0.5 * result.log_mu2_over_sij * result.log_mu2_over_sij;
    result.bracket_finite = (result.VS0 + result.VNS + result.finite0)
                          + (result.VS1 + result.finite1) * result.log_mu2_over_sij
                          + result.VS2 * result.log2_half;
    result.contribution_finite = -alpha_s / (2.0 * kPi) * result.bracket_finite * result.Bij;

    return result;
}

}  // namespace

CSIOperatorResult assemble_cdst_ioperator_finite(
    CSChannel channel,
    double shat,
    double mt,
    double mu_reg,
    double alpha_s,
    double costheta,
    const FourVector& p1,
    const FourVector& p2,
    const FourVector& p3,
    const FourVector& p4) {
    const double gs = std::sqrt(4.0 * kPi * alpha_s);
    const auto ordered_pairs = build_ordered_born_correlations(
        channel, shat, mt, mu_reg, gs, costheta, p1, p2, p3, p4);

    CSIOperatorResult result{};
    result.finite = 0.0;
    result.pairs.reserve(ordered_pairs.size());

    for (const auto& pair : ordered_pairs) {
        auto coefficients = compute_pair_coefficients(pair, mt, mu_reg, alpha_s);
        result.finite += coefficients.contribution_finite;
        result.pairs.push_back(coefficients);
    }

    return result;
}

}  // namespace ttbar_qtsub