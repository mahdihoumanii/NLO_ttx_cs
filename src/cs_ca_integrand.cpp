#include "cs_ca_integrand.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <stdexcept>

#include <LHAPDF/LHAPDF.h>

#include "born.hpp"
#include "color_correlation.hpp"
#include "core/constants.hpp"

namespace ttbar_qtsub {
namespace {

constexpr int kNf = 5;
constexpr double kTRLocal = 0.5;
constexpr double kPi2Over3 = kPiSq / 3.0;
constexpr double kPi2Over6 = kPiSq / 6.0;

struct CAGroupDefinition {
    std::string name;
    int type = 0;
    int emitter = 0;
    std::array<int, 3> in_collinear{};
    std::vector<CAPDFCombination> pdf;
};

struct SplittingValues {
    std::array<std::array<double, 3>, 4> Kbar{};
    std::array<std::array<double, 3>, 4> Kt{};
    std::array<std::array<double, 3>, 4> P{};
    std::array<std::array<double, 3>, 4> P_reg{};
    std::array<std::array<double, 3>, 4> Kbar_plus{};
    std::array<std::array<double, 3>, 4> Kt_plus{};
    std::array<std::array<double, 3>, 4> P_plus{};
    std::array<std::array<double, 3>, 4> intKbar_plus{};
    std::array<std::array<double, 3>, 4> intKt_plus{};
    std::array<std::array<double, 3>, 4> intP_plus{};
    std::array<double, 4> Kbar_delta{};
    std::array<double, 4> Kt_delta{};
    std::array<double, 4> P_delta{};
};

struct MassiveValues {
    std::array<std::array<std::array<double, 5>, 3>, 4> Kit{};
    std::array<std::array<std::array<double, 5>, 3>, 4> intKit_plus{};
    std::array<std::array<std::array<double, 5>, 3>, 4> Kit_plus_x{};
    std::array<std::array<std::array<double, 5>, 3>, 4> Kit_plus_1{};
    std::array<std::array<std::array<double, 5>, 3>, 4> Kit_plus_outside_x{};
    std::array<std::array<std::array<double, 5>, 3>, 4> Kit_plus_outside_1{};
    std::array<std::array<std::array<double, 5>, 3>, 4> Kit_delta{};
    std::array<std::array<double, 5>, 3> s{};
    std::array<std::array<double, 5>, 3> sx{};
};

double li2_series(double x) {
    double sum = 0.0;
    double x_power = x;
    for (int k = 1; k <= 500000; ++k) {
        const double term = x_power / (static_cast<double>(k) * static_cast<double>(k));
        sum += term;
        if (std::abs(term) < 1.0e-15) {
            break;
        }
        x_power *= x;
    }
    return sum;
}

double li2(double x) {
    if (x == 0.0) {
        return 0.0;
    }
    if (std::abs(x - 1.0) < 1.0e-14) {
        return kPi2Over6;
    }
    if (x < 0.0) {
        const double y = x / (x - 1.0);
        const double log_1_minus_x = std::log(1.0 - x);
        return -li2(y) - 0.5 * log_1_minus_x * log_1_minus_x;
    }
    if (x > 0.5 && x < 1.0) {
        return kPi2Over6 - std::log(x) * std::log(1.0 - x) - li2(1.0 - x);
    }
    if (x >= 0.0 && x <= 0.5) {
        return li2_series(x);
    }

    throw std::domain_error("li2: unsupported real branch");
}

LHAPDF::PDF& active_pdf() {
    static const int quiet_lhapdf = []() {
        LHAPDF::setVerbosity(0);
        return 0;
    }();
    (void)quiet_lhapdf;
    static std::unique_ptr<LHAPDF::PDF> pdf(LHAPDF::mkPDF("NNPDF31_nlo_as_0118", 0));
    return *pdf;
}

double xfx(int flavor, double x, double muF) {
    return active_pdf().xfxQ(flavor, x, muF);
}

double cs_gamma_g() {
    return (11.0 / 6.0) * kCA - (2.0 / 3.0) * kTRLocal * kNf;
}

double cs_gamma_q() {
    return (3.0 / 2.0) * kCF;
}

double P_qg(double x) {
    return kCF * (1.0 + sqr(1.0 - x)) / x;
}

double P_gq(double x) {
    return kTRLocal * (sqr(x) + sqr(1.0 - x));
}

double P_qq_reg(double x) {
    return -kCF * (1.0 + x);
}

double P_qq(double x) {
    return P_qq_reg(x);
}

double P_qq_plus(double x) {
    return kCF * 2.0 / (1.0 - x);
}

double intP_qq_plus(double x) {
    return -kCF * 2.0 * std::log(1.0 - x);
}

double P_gg_reg(double x) {
    return 2.0 * kCA * ((1.0 - x) / x - 1.0 + x * (1.0 - x));
}

double P_gg(double x) {
    return P_gg_reg(x);
}

double P_gg_plus(double x) {
    return 2.0 * kCA / (1.0 - x);
}

double intP_gg_plus(double x) {
    return -2.0 * kCA * std::log(1.0 - x);
}

double Kbar_qg(double x) {
    return P_qg(x) * std::log((1.0 - x) / x) + kCF * x;
}

double Kbar_gq(double x) {
    return P_gq(x) * std::log((1.0 - x) / x) + 2.0 * kTRLocal * x * (1.0 - x);
}

double Kbar_qq(double x) {
    return kCF * (-(1.0 + x) * std::log((1.0 - x) / x) + (1.0 - x));
}

double Kbar_qq_plus(double x) {
    return kCF * ((2.0 / (1.0 - x)) * std::log((1.0 - x) / x));
}

double intKbar_qq_plus(double x) {
    return kCF * (-sqr(std::log(1.0 - x)) - 2.0 * li2(1.0 - x) + kPi2Over3);
}

double Kbar_gg(double x) {
    return 2.0 * kCA * ((1.0 - x) / x - 1.0 + x * (1.0 - x)) * std::log((1.0 - x) / x);
}

double Kbar_gg_plus(double x) {
    return kCA * ((2.0 / (1.0 - x)) * std::log((1.0 - x) / x));
}

double intKbar_gg_plus(double x) {
    return kCA * (-sqr(std::log(1.0 - x)) - 2.0 * li2(1.0 - x) + kPi2Over3);
}

double Kt_qg(double x) {
    return P_qg(x) * std::log(1.0 - x);
}

double Kt_gq(double x) {
    return P_gq(x) * std::log(1.0 - x);
}

double Kt_qq(double x) {
    return P_qq_reg(x) * std::log(1.0 - x);
}

double Kt_qq_plus(double x) {
    return kCF * (2.0 / (1.0 - x)) * std::log(1.0 - x);
}

double intKt_qq_plus(double x) {
    return -kCF * sqr(std::log(1.0 - x));
}

double Kt_gg(double x) {
    return P_gg_reg(x) * std::log(1.0 - x);
}

double Kt_gg_plus(double x) {
    return kCA * (2.0 / (1.0 - x)) * std::log(1.0 - x);
}

double intKt_gg_plus(double x) {
    return -kCA * sqr(std::log(1.0 - x));
}

double born_me2(CSChannel channel, double shat, double mt, double muR, double gs, double costheta) {
    switch (channel) {
        case CSChannel::QQbar:
            return born_qq(shat, mt, muR, gs, costheta);
        case CSChannel::GG:
            return born_gg(shat, mt, muR, gs, costheta);
    }
    throw std::invalid_argument("born_me2: unsupported channel");
}

double color_me2(CSChannel channel,
                 double shat,
                 double mt,
                 double muR,
                 double gs,
                 double costheta,
                 int emitter,
                 int spectator) {
    switch (channel) {
        case CSChannel::QQbar:
            return color_correlation_qq(shat, mt, muR, gs, costheta, emitter, spectator);
        case CSChannel::GG:
            return color_correlation_gg(shat, mt, muR, gs, costheta, emitter, spectator);
    }
    throw std::invalid_argument("color_me2: unsupported channel");
}

double shat_from_momenta(const std::vector<FourVector>& momenta) {
    return project_sij(momenta[0], momenta[1], true, true);
}

double costheta_from_momenta(const std::vector<FourVector>& momenta, double mt) {
    const double shat = shat_from_momenta(momenta);
    const double beta = std::sqrt(1.0 - 4.0 * mt * mt / shat);
    double costheta = momenta[2].pz / (momenta[2].E * beta);
    if (costheta > 1.0 && costheta < 1.0 + kEpsilonCosTheta) {
        costheta = 1.0;
    }
    if (costheta < -1.0 && costheta > -1.0 - kEpsilonCosTheta) {
        costheta = -1.0;
    }
    return costheta;
}

std::vector<CAPDFCombination> expand_wildcards(const std::vector<CAPDFCombination>& combinations) {
    std::vector<CAPDFCombination> expanded;
    for (const auto& combination : combinations) {
        if (combination.flav1 == 10 && combination.flav2 == 10) {
            for (int flav1 = -kNf; flav1 <= kNf; ++flav1) {
                if (flav1 == 0) {
                    continue;
                }
                for (int flav2 = -kNf; flav2 <= kNf; ++flav2) {
                    if (flav2 == 0) {
                        continue;
                    }
                    expanded.push_back({combination.direction, flav1, flav2});
                }
            }
        }
        else if (combination.flav1 == 10) {
            for (int flav1 = -kNf; flav1 <= kNf; ++flav1) {
                if (flav1 != 0) {
                    expanded.push_back({combination.direction, flav1, combination.flav2});
                }
            }
        }
        else if (combination.flav2 == 10) {
            for (int flav2 = -kNf; flav2 <= kNf; ++flav2) {
                if (flav2 != 0) {
                    expanded.push_back({combination.direction, combination.flav1, flav2});
                }
            }
        }
        else {
            expanded.push_back(combination);
        }
    }
    return expanded;
}

std::string emitter_label(int emitter) {
    return emitter == 1 ? "a" : "b";
}

std::vector<CAGroupDefinition> build_groups(CABornSubprocess subprocess) {
    const auto born_pdf = ca_born_pdf_combinations(subprocess);
    const CAPDFCombination representative = born_pdf.front();
    const std::array<int, 3> representative_partons{{0, representative.flav1, representative.flav2}};
    std::vector<CAGroupDefinition> groups;

    for (int emitter = 1; emitter < 3; ++emitter) {
        for (int i_t = 0; i_t < 2; ++i_t) {
            std::vector<CAPDFCombination> pdf = born_pdf;
            CAGroupDefinition group;
            group.emitter = emitter;
            const bool hard_is_gluon = representative_partons[emitter] == 0;
            if (hard_is_gluon) {
                if (i_t == 0) {
                    group.type = 0;
                    group.in_collinear = {1, 0, 0};
                    group.in_collinear[emitter] = 1;
                    group.name = "C_" + emitter_label(emitter) + "^{gg}()";
                }
                else {
                    group.type = 2;
                    group.in_collinear = {0, 0, 0};
                    group.in_collinear[emitter] = 1;
                    for (auto& combination : pdf) {
                        if (emitter == 1) {
                            combination.flav1 = 10;
                        }
                        else {
                            combination.flav2 = 10;
                        }
                    }
                    group.name = "C_" + emitter_label(emitter) + "^{qg}()";
                }
            }
            else {
                if (i_t == 0) {
                    group.type = 3;
                    group.in_collinear = {0, 0, 0};
                    group.in_collinear[emitter] = 1;
                    for (auto& combination : pdf) {
                        if (emitter == 1) {
                            combination.flav1 = 0;
                        }
                        else {
                            combination.flav2 = 0;
                        }
                    }
                    group.name = "C_" + emitter_label(emitter) + "^{gq}()";
                }
                else {
                    group.type = 1;
                    group.in_collinear = {1, 0, 0};
                    group.in_collinear[emitter] = 1;
                    group.name = "C_" + emitter_label(emitter) + "^{qq}()";
                }
            }
            group.pdf = expand_wildcards(pdf);
            groups.push_back(group);
        }
    }

    return groups;
}

std::vector<int> spectators_for_emitter(int emitter) {
    if (emitter == 1) {
        return {0, 2, 3, 4};
    }
    return {0, 1, 3, 4};
}

std::array<std::array<int, 3>, 4> dipole_splitting(const std::vector<CAGroupDefinition>& groups) {
    std::array<std::array<int, 3>, 4> splitting{};
    for (const auto& group : groups) {
        for (int emitter = 0; emitter < 3; ++emitter) {
            if (group.in_collinear[emitter] == 1) {
                splitting[group.type][emitter] = 1;
            }
        }
    }
    return splitting;
}

SplittingValues splitting_values(const std::array<std::array<int, 3>, 4>& splitting,
                                 double x1,
                                 double x2,
                                 double z1,
                                 double z2) {
    SplittingValues values;
    values.P_delta[0] = (11.0 / 6.0) * kCA - (2.0 / 3.0) * kNf * kTRLocal;
    values.P_delta[1] = kCF * (3.0 / 2.0);
    values.Kbar_delta[0] = -(kCA * ((50.0 / 9.0) - kPiSq) - kTRLocal * kNf * (16.0 / 9.0));
    values.Kbar_delta[1] = kCF * (-(5.0 - kPiSq));
    values.Kt_delta[0] = kCA * (-kPiSq / 3.0);
    values.Kt_delta[1] = kCF * (-kPiSq / 3.0);

    for (int emitter = 1; emitter < 3; ++emitter) {
        const double z = emitter == 1 ? z1 : z2;
        const double x = emitter == 1 ? x1 : x2;
        if (splitting[0][emitter] == 1) {
            values.Kbar[0][emitter] = Kbar_gg(z);
            values.Kt[0][emitter] = Kt_gg(z);
            values.P[0][emitter] = P_gg(z);
            values.P_reg[0][emitter] = P_gg_reg(z);
            values.Kbar_plus[0][emitter] = Kbar_gg_plus(z);
            values.Kt_plus[0][emitter] = Kt_gg_plus(z);
            values.P_plus[0][emitter] = P_gg_plus(z);
            values.intKbar_plus[0][emitter] = intKbar_gg_plus(x);
            values.intKt_plus[0][emitter] = intKt_gg_plus(x);
            values.intP_plus[0][emitter] = intP_gg_plus(x);
        }
        if (splitting[1][emitter] == 1) {
            values.Kbar[1][emitter] = Kbar_qq(z);
            values.Kt[1][emitter] = Kt_qq(z);
            values.P[1][emitter] = P_qq(z);
            values.P_reg[1][emitter] = P_qq_reg(z);
            values.Kbar_plus[1][emitter] = Kbar_qq_plus(z);
            values.Kt_plus[1][emitter] = Kt_qq_plus(z);
            values.P_plus[1][emitter] = P_qq_plus(z);
            values.intKbar_plus[1][emitter] = intKbar_qq_plus(x);
            values.intKt_plus[1][emitter] = intKt_qq_plus(x);
            values.intP_plus[1][emitter] = intP_qq_plus(x);
        }
        if (splitting[2][emitter] == 1) {
            values.Kbar[2][emitter] = Kbar_qg(z);
            values.Kt[2][emitter] = Kt_qg(z);
            values.P[2][emitter] = P_qg(z);
            values.P_reg[2][emitter] = values.P[2][emitter];
        }
        if (splitting[3][emitter] == 1) {
            values.Kbar[3][emitter] = Kbar_gq(z);
            values.Kt[3][emitter] = Kt_gq(z);
            values.P[3][emitter] = P_gq(z);
            values.P_reg[3][emitter] = values.P[3][emitter];
        }
    }
    return values;
}

std::vector<std::array<int, 2>> unique_pairs(const std::vector<CAGroupDefinition>& groups) {
    std::vector<std::array<int, 2>> pairs;
    for (const auto& group : groups) {
        for (const int spectator : spectators_for_emitter(group.emitter)) {
            if (spectator == 0) {
                continue;
            }
            std::array<int, 2> pair{{std::min(group.emitter, spectator), std::max(group.emitter, spectator)}};
            if (std::find(pairs.begin(), pairs.end(), pair) == pairs.end()) {
                pairs.push_back(pair);
            }
        }
    }
    return pairs;
}

double z_for_emitter(const CAInput& input, int emitter) {
    return emitter == 1 ? input.z1 : input.z2;
}

double x_for_emitter(const CAInput& input, int emitter) {
    return emitter == 1 ? input.x1 : input.x2;
}

double gz_for_emitter(const CAInput& input, int emitter) {
    return emitter == 1 ? input.g_z1 : input.g_z2;
}

MassiveValues massive_values(const CAInput& input,
                             const std::vector<CAGroupDefinition>& groups,
                             const std::array<std::array<int, 3>, 4>& splitting) {
    MassiveValues values;
    const std::array<FourVector, 5> legs{{{}, input.momenta[0], input.momenta[1], input.momenta[2], input.momenta[3]}};
    const auto pairs = unique_pairs(groups);
    const double m = input.mt;
    const double m2 = m * m;

    for (const auto& pair : pairs) {
        const int emitter = pair[0];
        const int spectator = pair[1];
        values.s[emitter][spectator] = matrix_sij(legs[emitter], legs[spectator]);
        if (spectator < 3) {
            continue;
        }
        const double z = z_for_emitter(input, emitter);
        const double x = x_for_emitter(input, emitter);
        const double s = values.s[emitter][spectator];
        values.sx[emitter][spectator] = s / z;
        const double sx = values.sx[emitter][spectator];
        const double mu2 = m2 / s;
        const double mu2_x = m2 / sx;

        for (int type = 0; type < 4; ++type) {
            if (splitting[type][1] != 1 && splitting[type][2] != 1) {
                continue;
            }
            if (type == 0) {
                values.Kit[type][emitter][spectator] =
                    -2.0 * std::log(2.0 - z) / (1.0 - z)
                    + 2.0 * m2 / (z * sx) * std::log(m2 / ((1.0 - z) * sx + m2));
                values.Kit_plus_x[type][emitter][spectator] =
                    2.0 * std::log(1.0 - z) / (1.0 - z)
                    + (1.0 - z) / (2.0 * sqr(1.0 - z + mu2_x))
                    - 2.0 / (1.0 - z) * (1.0 + std::log(1.0 - z + mu2_x));
                values.Kit_plus_1[type][emitter][spectator] =
                    2.0 * std::log(1.0 - z) / (1.0 - z)
                    + (1.0 - z) / (2.0 * sqr(1.0 - z + mu2))
                    - 2.0 / (1.0 - z) * (1.0 + std::log(1.0 - z + mu2));
                values.intKit_plus[type][emitter][spectator] =
                    -sqr(std::log(1.0 - x))
                    + 0.5 * (-mu2 / (1.0 - x + mu2) + mu2 / (1.0 + mu2) - std::log((1.0 - x + mu2) / (1.0 + mu2)))
                    + 2.0 * (li2(-1.0 / mu2) - li2(-(1.0 - x) / mu2) + std::log(1.0 - x) * (1.0 + std::log(mu2)));
                values.Kit_plus_outside_x[type][emitter][spectator] =
                    std::log(((2.0 - z) * sx) / ((2.0 - z) * sx + m2))
                    + std::log(2.0 + mu2_x - z);
                values.Kit_plus_outside_1[type][emitter][spectator] =
                    std::log(s / (s + m2)) + std::log(1.0 + mu2);
                values.Kit_delta[type][emitter][spectator] =
                    -cs_gamma_q() / kCF
                    + mu2 * std::log(m2 / (s + m2))
                    + 0.5 * m2 / (s + m2);
            }
            else if (type == 1) {
                values.Kit[type][emitter][spectator] = -2.0 * std::log(2.0 - z) / (1.0 - z);
                values.Kit_plus_x[type][emitter][spectator] =
                    2.0 * std::log(1.0 - z) / (1.0 - z)
                    + (1.0 - z) / (2.0 * sqr(1.0 - z + mu2_x))
                    - 2.0 / (1.0 - z) * (1.0 + std::log(1.0 - z + mu2_x));
                values.Kit_plus_1[type][emitter][spectator] =
                    2.0 * std::log(1.0 - z) / (1.0 - z)
                    + (1.0 - z) / (2.0 * sqr(1.0 - z + mu2))
                    - 2.0 / (1.0 - z) * (1.0 + std::log(1.0 - z + mu2));
                values.intKit_plus[type][emitter][spectator] =
                    -sqr(std::log(1.0 - x))
                    + 0.5 * (-mu2 / (1.0 - x + mu2) + mu2 / (1.0 + mu2) - std::log((1.0 - x + mu2) / (1.0 + mu2)))
                    + 2.0 * (li2(-1.0 / mu2) - li2(-(1.0 - x) / mu2) + std::log(1.0 - x) * (1.0 + std::log(mu2)));
                values.Kit_plus_outside_x[type][emitter][spectator] =
                    std::log(((2.0 - z) * sx) / ((2.0 - z) * sx + m2))
                    + std::log(2.0 + mu2_x - z);
                values.Kit_plus_outside_1[type][emitter][spectator] =
                    std::log(s / (s + m2)) + std::log(1.0 + mu2);
                values.Kit_delta[type][emitter][spectator] =
                    -cs_gamma_q() / kCF
                    + mu2 * std::log(m2 / (s + m2))
                    + 0.5 * m2 / (s + m2);
            }
            else if (type == 2) {
                values.Kit[type][emitter][spectator] =
                    2.0 * kCF / kCA * m2 / (z * sx) * std::log(m2 / ((1.0 - z) * sx + m2));
            }
        }
    }
    return values;
}

std::array<double, 3> pdf_factor_raw(const CAGroupDefinition& group,
                                     double x1,
                                     double x2,
                                     double z1,
                                     double z2,
                                     double muF) {
    std::array<double, 3> result{};
    for (int iz = 0; iz < 3; ++iz) {
        if (group.in_collinear[iz] != 1) {
            continue;
        }
        const double px1 = iz == 1 ? x1 / z1 : x1;
        const double px2 = iz == 2 ? x2 / z2 : x2;
        for (const auto& combination : group.pdf) {
            result[iz] += xfx(combination.flav1, px1, muF) * xfx(combination.flav2, px2, muF);
        }
    }
    return result;
}

std::array<double, 3> pdf_factor_cv(const CAGroupDefinition& group,
                                    double x1,
                                    double x2,
                                    double z1,
                                    double z2,
                                    double muF) {
    std::array<double, 3> result = pdf_factor_raw(group, x1, x2, z1, z2, muF);
    const double tau = x1 * x2;
    for (int iz = 0; iz < 3; ++iz) {
        if (group.in_collinear[iz] != 1) {
            continue;
        }
        result[iz] /= tau;
        if (iz == 1) {
            result[iz] *= z1;
        }
        else if (iz == 2) {
            result[iz] *= z2;
        }
    }
    return result;
}

std::array<double, 3> assemble_integrand_slots(const CAInput& input,
                                               const CAGroupDefinition& group,
                                               const std::array<double, 3>& kp,
                                               const std::array<double, 3>& pdf_factor) {
    std::array<double, 3> slots{};
    const double gz = gz_for_emitter(input, group.emitter);
    slots[0] = kp[0] / gz * pdf_factor[group.emitter];
    if (group.in_collinear[0] == 1) {
        slots[1] = (kp[1] / gz + kp[2]) * pdf_factor[0];
    }
    slots[2] = slots[0] + slots[1];
    return slots;
}

void validate_input(const CAInput& input) {
    if (input.momenta.size() != 4) {
        throw std::invalid_argument("compute_ca_integrand: expected four Born momenta");
    }
    if (!(input.x1 > 0.0 && input.x1 <= 1.0) || !(input.x2 > 0.0 && input.x2 <= 1.0)) {
        throw std::invalid_argument("compute_ca_integrand: x1/x2 must be in (0, 1]");
    }
    if (!(input.z1 >= input.x1 && input.z1 <= 1.0) || !(input.z2 >= input.x2 && input.z2 <= 1.0)) {
        throw std::invalid_argument("compute_ca_integrand: z1/z2 must be in [x, 1]");
    }
    if (!(input.alpha_s > 0.0) || !(input.muF > 0.0) || !(input.muR > 0.0) || !(input.mt > 0.0)) {
        throw std::invalid_argument("compute_ca_integrand: alpha_s, scales, and mt must be positive");
    }
}

}  // namespace

CSChannel ca_channel(CABornSubprocess subprocess) {
    switch (subprocess) {
        case CABornSubprocess::GG_TTbar:
            return CSChannel::GG;
        case CABornSubprocess::DDbar_TTbar:
        case CABornSubprocess::UUbar_TTbar:
        case CABornSubprocess::BBbar_TTbar:
            return CSChannel::QQbar;
    }
    throw std::invalid_argument("ca_channel: unsupported subprocess");
}

const char* ca_subprocess_name(CABornSubprocess subprocess) {
    switch (subprocess) {
        case CABornSubprocess::GG_TTbar:
            return "gg_tt~";
        case CABornSubprocess::DDbar_TTbar:
            return "dd~_tt~";
        case CABornSubprocess::UUbar_TTbar:
            return "uu~_tt~";
        case CABornSubprocess::BBbar_TTbar:
            return "bb~_tt~";
    }
    return "unknown";
}

std::vector<CAPDFCombination> ca_born_pdf_combinations(CABornSubprocess subprocess) {
    switch (subprocess) {
        case CABornSubprocess::GG_TTbar:
            return {{1, 0, 0}};
        case CABornSubprocess::DDbar_TTbar:
            return {{1, 1, -1}, {1, 3, -3}, {-1, 1, -1}, {-1, 3, -3}};
        case CABornSubprocess::UUbar_TTbar:
            return {{1, 2, -2}, {1, 4, -4}, {-1, 2, -2}, {-1, 4, -4}};
        case CABornSubprocess::BBbar_TTbar:
            return {{1, 5, -5}, {-1, 5, -5}};
    }
    throw std::invalid_argument("ca_born_pdf_combinations: unsupported subprocess");
}

CAResult compute_ca_integrand(const CAInput& input) {
    validate_input(input);

    const CSChannel channel = ca_channel(input.subprocess);
    const double shat = shat_from_momenta(input.momenta);
    const double costheta = costheta_from_momenta(input.momenta, input.mt);
    const double gs = std::sqrt(4.0 * kPi * input.alpha_s);
    const double born = born_me2(channel, shat, input.mt, input.muR, gs, costheta);

    const auto groups = build_groups(input.subprocess);
    const auto splitting = dipole_splitting(groups);
    const SplittingValues split = splitting_values(splitting, input.x1, input.x2, input.z1, input.z2);
    const MassiveValues massive = massive_values(input, groups, splitting);
    const std::array<FourVector, 5> legs{{{}, input.momenta[0], input.momenta[1], input.momenta[2], input.momenta[3]}};

    CAResult result;
    result.born_me2 = born;
    result.shat = shat;
    result.costheta = costheta;

    const double alpha_s_2pi = input.alpha_s / (2.0 * kPi);
    const std::array<double, 3> iT2_ap{{0.0,
        channel == CSChannel::GG ? 1.0 / kCA : 1.0 / kCF,
        channel == CSChannel::GG ? 1.0 / kCA : 1.0 / kCF}};
    const std::array<double, 3> gamma_a_T2_ap{{0.0,
        channel == CSChannel::GG ? cs_gamma_g() / kCA : cs_gamma_q() / kCF,
        channel == CSChannel::GG ? cs_gamma_g() / kCA : cs_gamma_q() / kCF}};

    for (const auto& group : groups) {
        CACollinearContribution contribution;
        contribution.name = group.name;
        contribution.type = group.type;
        contribution.emitter = group.emitter;
        contribution.in_collinear = group.in_collinear;
        contribution.spectators = spectators_for_emitter(group.emitter);
        contribution.me2_cf.reserve(contribution.spectators.size());
        for (const int spectator : contribution.spectators) {
            contribution.me2_cf.push_back(spectator == 0
                ? born
                : color_me2(channel, shat, input.mt, input.muR, gs, costheta, group.emitter, spectator));
        }

        for (std::size_t j = 0; j < contribution.spectators.size(); ++j) {
            const int spectator = contribution.spectators[j];
            const double me2 = contribution.me2_cf[j];
            const int type = group.type;
            const int emitter = group.emitter;
            const double z = z_for_emitter(input, emitter);
            const double x = x_for_emitter(input, emitter);

            if (spectator == 0) {
                contribution.sum_K[0] += (split.Kbar[type][emitter] + split.Kbar_plus[type][emitter]) * me2;
                contribution.sum_K[1] += (-split.Kbar_plus[type][emitter]) * z * me2;
                contribution.sum_K[2] += (split.Kbar_delta[type] - split.intKbar_plus[type][emitter]) * me2;
            }
            else if (spectator < 3) {
                contribution.sum_K[0] += -iT2_ap[emitter] * (split.Kt[type][emitter] + split.Kt_plus[type][emitter]) * me2;
                contribution.sum_K[1] += -iT2_ap[emitter] * (-split.Kt_plus[type][emitter]) * z * me2;
                contribution.sum_K[2] += -iT2_ap[emitter] * (split.Kt_delta[type] - split.intKt_plus[type][emitter]) * me2;
            }
            else {
                const double m2 = input.mt * input.mt;
                const std::array<int, 2> pair{{std::min(emitter, spectator), std::max(emitter, spectator)}};
                const double s = massive.s[pair[0]][pair[1]];
                const double sx = massive.sx[pair[0]][pair[1]];
                contribution.sum_K[0] += -(
                    massive.Kit[type][emitter][spectator]
                    + massive.Kit_plus_x[type][emitter][spectator]
                    + 2.0 / (1.0 - z) * massive.Kit_plus_outside_x[type][emitter][spectator]) * me2;
                contribution.sum_K[1] += -(
                    -massive.Kit_plus_1[type][emitter][spectator]
                    - 2.0 / (1.0 - z) * massive.Kit_plus_outside_1[type][emitter][spectator]) * z * me2;
                contribution.sum_K[2] += -(
                    massive.Kit_delta[type][emitter][spectator]
                    - massive.intKit_plus[type][emitter][spectator]
                    - 2.0 * std::log(1.0 - x) * massive.Kit_plus_outside_1[type][emitter][spectator]) * me2;

                contribution.sum_K[0] += -iT2_ap[emitter] * split.P_reg[type][emitter]
                    * std::log(((1.0 - z) * sx) / ((1.0 - z) * sx + m2)) * me2;
                if (group.in_collinear[0] == 1) {
                    contribution.sum_K[2] += -gamma_a_T2_ap[emitter]
                        * (std::log((s - 2.0 * input.mt * std::sqrt(s + m2) + 2.0 * m2) / s)
                           + 2.0 * input.mt / (std::sqrt(s + m2) + input.mt)) * me2;
                }
            }

            if (spectator != 0) {
                const std::array<int, 2> pair{{std::min(emitter, spectator), std::max(emitter, spectator)}};
                const double ln_papi = std::log(matrix_sij(legs[pair[0]], legs[pair[1]]));
                const double logscale2_fact_papi = 2.0 * std::log(input.muF) - ln_papi;
                contribution.sum_P[0] += iT2_ap[emitter] * (split.P[type][emitter] + split.P_plus[type][emitter]) * logscale2_fact_papi * me2;
                contribution.sum_P[1] += iT2_ap[emitter] * (-split.P_plus[type][emitter]) * z * logscale2_fact_papi * me2;
                contribution.sum_P[2] += iT2_ap[emitter] * (split.P_delta[type] - split.intP_plus[type][emitter]) * logscale2_fact_papi * me2;
            }
        }

        for (int slot = 0; slot < 3; ++slot) {
            contribution.kp[slot] = alpha_s_2pi * (contribution.sum_K[slot] + contribution.sum_P[slot]);
        }
        contribution.pdf_factor_tsv = pdf_factor_raw(group, input.x1, input.x2, input.z1, input.z2, input.muF);
        contribution.pdf_factor_cv = pdf_factor_cv(group, input.x1, input.x2, input.z1, input.z2, input.muF);

        const auto tsv_slots = assemble_integrand_slots(input, group, contribution.kp, contribution.pdf_factor_tsv);
        const auto cv_slots = assemble_integrand_slots(input, group, contribution.kp, contribution.pdf_factor_cv);
        contribution.integrand_tsv = tsv_slots[2];
        contribution.integrand_cv = cv_slots[2];
        result.sum_integrand_tsv += contribution.integrand_tsv;
        result.sum_integrand_cv += contribution.integrand_cv;
        result.collinear.push_back(contribution);
    }

    return result;
}

}  // namespace ttbar_qtsub