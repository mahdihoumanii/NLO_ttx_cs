#include "color_correlation.hpp"
#include "born.hpp"
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

#include "common.hpp"

namespace ttbar_qtsub {
namespace {

using Complex = std::complex<double>;

constexpr Complex kI(0.0, 1.0);
constexpr Complex kPiC(kPi, 0.0);

double li_series_real(double x, int power) {
    if (x == 0.0) {
        return 0.0;
    }
double     sum = 0.0;
    double xk = x;
    for (int k = 1; k <= 500000; ++k) {
        const double denom = (power == 2) ? static_cast<double>(k) * static_cast<double>(k)
                                           : static_cast<double>(k) * static_cast<double>(k) * static_cast<double>(k);
        const double term = xk / denom;
        sum += term;
        if (std::abs(term) < 1.0e-15) {
            break;
        }
xk         *= x;
    }
return     sum;
}

Complex li2(const Complex& z) {
    const double x = std::real(z);
    constexpr double zeta2 = (kPi * kPi) / 6.0;

    if (std::abs(x - 1.0) < 1.0e-14) {
        return Complex(zeta2, 0.0);
    }
if     (x > 1.0) {
        const double li2_1mx = std::real(li2(Complex(1.0 - x, 0.0)));
        const double lx = std::log(x);
        const double lxm1 = std::log(x - 1.0);
        return Complex(zeta2 - lx * lxm1 - li2_1mx, -kPi * lx);
    }
if     (x < 0.0) {
        const double y = x / (x - 1.0);
        const double ly = std::log(1.0 - x);
        return -li2(Complex(y, 0.0)) - Complex(0.5 * ly * ly, 0.0);
    }
if     (x > 0.5) {
        const double l1mx = std::log(1.0 - x);
        return Complex(zeta2, 0.0) - Complex(std::log(x) * l1mx, 0.0) - li2(Complex(1.0 - x, 0.0));
    }

return     Complex(li_series_real(x, 2), 0.0);
}

Complex li3(const Complex& z) {
    const double x = std::real(z);
    constexpr double zeta3 = 1.2020569031595942854;

    if (std::abs(x - 1.0) < 1.0e-14) {
        return Complex(zeta3, 0.0);
    }
if     (x > 1.0) {
        const Complex lminusx = std::log(Complex(-x, 0.0));
        return li3(Complex(1.0 / x, 0.0)) - (lminusx * lminusx * lminusx) / 6.0 - (kPi * kPi / 6.0) * lminusx;
    }
if     (x < -1.0) {
        const Complex lminusx = std::log(Complex(-x, 0.0));
        return li3(Complex(1.0 / x, 0.0)) - (lminusx * lminusx * lminusx) / 6.0 - (kPi * kPi / 6.0) * lminusx;
    }

return     Complex(li_series_real(x, 3), 0.0);
}

}  // namespace
double color_correlation_qq_1_2(double shat, double mt_in, double mu_in, double gs_in,double costheta) {
    const Complex s(shat, 0.0);
    const Complex mt(mt_in, 0.0);
    const Complex mu(mu_in, 0.0);
    const Complex gs(gs_in, 0.0);
	const Complex c(costheta, 0.0);
    if (shat <= 4.0 * mt_in * mt_in) {
        return 0.0;
    }
    const Complex beta = std::sqrt(s * (s - 4.0 * std::pow(mt, 2.0))) / s;
	const Complex expr = 0.03703703703703704*(2. - 1.*std::pow(beta,2) + std::pow(beta,2)*std::pow(c,2))*std::pow(gs,4);

    return std::real(expr);

	}

double color_correlation_qq_1_3(double shat, double mt_in, double mu_in, double gs_in,double costheta) {
    const Complex s(shat, 0.0);
    const Complex mt(mt_in, 0.0);
    const Complex mu(mu_in, 0.0);
    const Complex gs(gs_in, 0.0);
	const Complex c(costheta, 0.0);
    if (shat <= 4.0 * mt_in * mt_in) {
        return 0.0;
    }
    const Complex beta = std::sqrt(s * (s - 4.0 * std::pow(mt, 2.0))) / s;
	const Complex expr = -0.2592592592592593*(2. - 1.*std::pow(beta,2) + std::pow(beta,2)*std::pow(c,2))*std::pow(gs,4);

    return std::real(expr);

	}

double color_correlation_qq_1_4(double shat, double mt_in, double mu_in, double gs_in,double costheta) {
    const Complex s(shat, 0.0);
    const Complex mt(mt_in, 0.0);
    const Complex mu(mu_in, 0.0);
    const Complex gs(gs_in, 0.0);
	const Complex c(costheta, 0.0);
    if (shat <= 4.0 * mt_in * mt_in) {
        return 0.0;
    }
    const Complex beta = std::sqrt(s * (s - 4.0 * std::pow(mt, 2.0))) / s;
	const Complex expr = -0.07407407407407407*(2. - 1.*std::pow(beta,2) + std::pow(beta,2)*std::pow(c,2))*std::pow(gs,4);

    return std::real(expr);

	}

double color_correlation_qq_2_3(double shat, double mt_in, double mu_in, double gs_in,double costheta) {
    const Complex s(shat, 0.0);
    const Complex mt(mt_in, 0.0);
    const Complex mu(mu_in, 0.0);
    const Complex gs(gs_in, 0.0);
	const Complex c(costheta, 0.0);
    if (shat <= 4.0 * mt_in * mt_in) {
        return 0.0;
    }
    const Complex beta = std::sqrt(s * (s - 4.0 * std::pow(mt, 2.0))) / s;
	const Complex expr = -0.07407407407407407*(2. - 1.*std::pow(beta,2) + std::pow(beta,2)*std::pow(c,2))*std::pow(gs,4);

    return std::real(expr);

	}

double color_correlation_qq_2_4(double shat, double mt_in, double mu_in, double gs_in,double costheta) {
    const Complex s(shat, 0.0);
    const Complex mt(mt_in, 0.0);
    const Complex mu(mu_in, 0.0);
    const Complex gs(gs_in, 0.0);
	const Complex c(costheta, 0.0);
    if (shat <= 4.0 * mt_in * mt_in) {
        return 0.0;
    }
    const Complex beta = std::sqrt(s * (s - 4.0 * std::pow(mt, 2.0))) / s;
	const Complex expr = -0.2592592592592593*(2. - 1.*std::pow(beta,2) + std::pow(beta,2)*std::pow(c,2))*std::pow(gs,4);

    return std::real(expr);

	}

double color_correlation_qq_3_4(double shat, double mt_in, double mu_in, double gs_in,double costheta) {
    const Complex s(shat, 0.0);
    const Complex mt(mt_in, 0.0);
    const Complex mu(mu_in, 0.0);
    const Complex gs(gs_in, 0.0);
	const Complex c(costheta, 0.0);
    if (shat <= 4.0 * mt_in * mt_in) {
        return 0.0;
    }
    const Complex beta = std::sqrt(s * (s - 4.0 * std::pow(mt, 2.0))) / s;
	const Complex expr = 0.03703703703703704*(2. - 1.*std::pow(beta,2) + std::pow(beta,2)*std::pow(c,2))*std::pow(gs,4);

    return std::real(expr);

	}

double color_correlation_gg_1_2(double shat, double mt_in, double mu_in, double gs_in,double costheta) {
    const Complex s(shat, 0.0);
    const Complex mt(mt_in, 0.0);
    const Complex mu(mu_in, 0.0);
    const Complex gs(gs_in, 0.0);
	const Complex c(costheta, 0.0);
    if (shat <= 4.0 * mt_in * mt_in) {
        return 0.0;
    }
    const Complex beta = std::sqrt(s * (s - 4.0 * std::pow(mt, 2.0))) / s;
	const Complex expr = (0.28125*(1. + std::pow(beta,2)*std::pow(c,2))*(-1. - 2.*std::pow(beta,2) + 2.*std::pow(beta,4) + 2.*std::pow(beta,2)*std::pow(c,2) - 2.*std::pow(beta,4)*std::pow(c,2) + std::pow(beta,4)*std::pow(c,4))*std::pow(gs,4))/(std::pow(-1. + beta*c,2)*std::pow(1. + beta*c,2));

    return std::real(expr);

	}

double color_correlation_gg_1_3(double shat, double mt_in, double mu_in, double gs_in,double costheta) {
    const Complex s(shat, 0.0);
    const Complex mt(mt_in, 0.0);
    const Complex mu(mu_in, 0.0);
    const Complex gs(gs_in, 0.0);
	const Complex c(costheta, 0.0);
    if (shat <= 4.0 * mt_in * mt_in) {
        return 0.0;
    }
    const Complex beta = std::sqrt(s * (s - 4.0 * std::pow(mt, 2.0))) / s;
	const Complex expr = (0.015625*(1. + 3.*beta*c)*(5. + 3.*beta*c)*(-1. - 2.*std::pow(beta,2) + 2.*std::pow(beta,4) + 2.*std::pow(beta,2)*std::pow(c,2) - 2.*std::pow(beta,4)*std::pow(c,2) + std::pow(beta,4)*std::pow(c,4))*std::pow(gs,4))/(std::pow(-1. + beta*c,2)*std::pow(1. + beta*c,2));

    return std::real(expr);

	}

double color_correlation_gg_1_4(double shat, double mt_in, double mu_in, double gs_in,double costheta) {
    const Complex s(shat, 0.0);
    const Complex mt(mt_in, 0.0);
    const Complex mu(mu_in, 0.0);
    const Complex gs(gs_in, 0.0);
	const Complex c(costheta, 0.0);
    if (shat <= 4.0 * mt_in * mt_in) {
        return 0.0;
    }
    const Complex beta = std::sqrt(s * (s - 4.0 * std::pow(mt, 2.0))) / s;
	const Complex expr = (0.015625*(-5. + 3.*beta*c)*(-1. + 3.*beta*c)*(-1. - 2.*std::pow(beta,2) + 2.*std::pow(beta,4) + 2.*std::pow(beta,2)*std::pow(c,2) - 2.*std::pow(beta,4)*std::pow(c,2) + std::pow(beta,4)*std::pow(c,4))*std::pow(gs,4))/(std::pow(-1. + beta*c,2)*std::pow(1. + beta*c,2));

    return std::real(expr);

	}

double color_correlation_gg_2_3(double shat, double mt_in, double mu_in, double gs_in,double costheta) {
    const Complex s(shat, 0.0);
    const Complex mt(mt_in, 0.0);
    const Complex mu(mu_in, 0.0);
    const Complex gs(gs_in, 0.0);
	const Complex c(costheta, 0.0);
    if (shat <= 4.0 * mt_in * mt_in) {
        return 0.0;
    }
    const Complex beta = std::sqrt(s * (s - 4.0 * std::pow(mt, 2.0))) / s;
	const Complex expr = (0.015625*(-5. + 3.*beta*c)*(-1. + 3.*beta*c)*(-1. - 2.*std::pow(beta,2) + 2.*std::pow(beta,4) + 2.*std::pow(beta,2)*std::pow(c,2) - 2.*std::pow(beta,4)*std::pow(c,2) + std::pow(beta,4)*std::pow(c,4))*std::pow(gs,4))/(std::pow(-1. + beta*c,2)*std::pow(1. + beta*c,2));

    return std::real(expr);

	}

double color_correlation_gg_2_4(double shat, double mt_in, double mu_in, double gs_in,double costheta) {
    const Complex s(shat, 0.0);
    const Complex mt(mt_in, 0.0);
    const Complex mu(mu_in, 0.0);
    const Complex gs(gs_in, 0.0);
	const Complex c(costheta, 0.0);
    if (shat <= 4.0 * mt_in * mt_in) {
        return 0.0;
    }
    const Complex beta = std::sqrt(s * (s - 4.0 * std::pow(mt, 2.0))) / s;
	const Complex expr = (0.015625*(1. + 3.*beta*c)*(5. + 3.*beta*c)*(-1. - 2.*std::pow(beta,2) + 2.*std::pow(beta,4) + 2.*std::pow(beta,2)*std::pow(c,2) - 2.*std::pow(beta,4)*std::pow(c,2) + std::pow(beta,4)*std::pow(c,4))*std::pow(gs,4))/(std::pow(-1. + beta*c,2)*std::pow(1. + beta*c,2));

    return std::real(expr);

	}

double color_correlation_gg_3_4(double shat, double mt_in, double mu_in, double gs_in,double costheta) {
    const Complex s(shat, 0.0);
    const Complex mt(mt_in, 0.0);
    const Complex mu(mu_in, 0.0);
    const Complex gs(gs_in, 0.0);
	const Complex c(costheta, 0.0);
    if (shat <= 4.0 * mt_in * mt_in) {
        return 0.0;
    }
    const Complex beta = std::sqrt(s * (s - 4.0 * std::pow(mt, 2.0))) / s;
	const Complex expr = (-0.003472222222222222*(-11. + 9.*std::pow(beta,2)*std::pow(c,2))*(-1. - 2.*std::pow(beta,2) + 2.*std::pow(beta,4) + 2.*std::pow(beta,2)*std::pow(c,2) - 2.*std::pow(beta,4)*std::pow(c,2) + std::pow(beta,4)*std::pow(c,4))*std::pow(gs,4))/(std::pow(-1. + beta*c,2)*std::pow(1. + beta*c,2));

    return std::real(expr);

	}

// QQ dispatcher: color_correlation_qq(shat, mt, mu, gs, c, i, j)
// Diagonal (i,i): CF * born_qq   [CF = 4/3]
// Off-diagonal: symmetric, delegates to specific function
// Valid particle indices: 1..4
double color_correlation_qq(double shat, double mt, double mu, double gs, double c, int i, int j) {
    constexpr double CF = 4.0 / 3.0;
    if (i == j) {
        return CF * born_qq(shat, mt, mu, gs, c);
    }
    if (i > j) { int tmp = i; i = j; j = tmp; }  // symmetrise
    switch (i * 10 + j) {
        case 12: return color_correlation_qq_1_2(shat, mt, mu, gs, c);
        case 13: return color_correlation_qq_1_3(shat, mt, mu, gs, c);
        case 14: return color_correlation_qq_1_4(shat, mt, mu, gs, c);
        case 23: return color_correlation_qq_2_3(shat, mt, mu, gs, c);
        case 24: return color_correlation_qq_2_4(shat, mt, mu, gs, c);
        case 34: return color_correlation_qq_3_4(shat, mt, mu, gs, c);
        default: throw std::invalid_argument("color_correlation_qq: invalid (i,j) pair");
    }
}

// GG dispatcher: color_correlation_gg(shat, mt, mu, gs, c, i, j)
// Diagonal (i,i): emitter Casimir times born_gg
//   initial gluons 1,2 -> CA * born_gg
//   final heavy quarks 3,4 -> CF * born_gg
// Off-diagonal: symmetric, delegates to specific function
// Valid particle indices: 1..4
double color_correlation_gg(double shat, double mt, double mu, double gs, double c, int i, int j) {
    constexpr double CA = 3.0;
    constexpr double CF = 4.0 / 3.0;
    if (i == j) {
        if (i == 1 || i == 2) {
            return CA * born_gg(shat, mt, mu, gs, c);
        }
        if (i == 3 || i == 4) {
            return CF * born_gg(shat, mt, mu, gs, c);
        }
        throw std::invalid_argument("color_correlation_gg: invalid diagonal index");
    }
    if (i > j) { int tmp = i; i = j; j = tmp; }  // symmetrise
    switch (i * 10 + j) {
        case 12: return color_correlation_gg_1_2(shat, mt, mu, gs, c);
        case 13: return color_correlation_gg_1_3(shat, mt, mu, gs, c);
        case 14: return color_correlation_gg_1_4(shat, mt, mu, gs, c);
        case 23: return color_correlation_gg_2_3(shat, mt, mu, gs, c);
        case 24: return color_correlation_gg_2_4(shat, mt, mu, gs, c);
        case 34: return color_correlation_gg_3_4(shat, mt, mu, gs, c);
        default: throw std::invalid_argument("color_correlation_gg: invalid (i,j) pair");
    }
}

}
