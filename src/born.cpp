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
double born_qq(double shat, double mt_in, double mu_in, double gs_in,double costheta) {
    const Complex s(shat, 0.0);
    const Complex mt(mt_in, 0.0);
    const Complex mu(mu_in, 0.0);
    const Complex gs(gs_in, 0.0);
	const Complex c(costheta, 0.0);
    if (shat <= 4.0 * mt_in * mt_in) {
        return 0.0;
    }
    const Complex beta = std::sqrt(s * (s - 4.0 * std::pow(mt, 2.0))) / s;
	const Complex expr = 0.2222222222222222*(2. - 1.*std::pow(beta,2) + std::pow(beta,2)*std::pow(c,2))*std::pow(gs,4);

    return std::real(expr);

	}

double born_gg(double shat, double mt_in, double mu_in, double gs_in, double costheta) {
    const Complex s(shat, 0.0);
    const Complex mt(mt_in, 0.0);
    const Complex mu(mu_in, 0.0);
    const Complex gs(gs_in, 0.0);
	const Complex c(costheta, 0.0);
    if (shat <= 4.0 * mt_in * mt_in) {
        return 0.0;
    }
    const Complex beta = std::sqrt(s * (s - 4.0 * std::pow(mt, 2.0))) / s;
	const Complex expr = (-0.02083333333333333*(7. + 9.*std::pow(beta,2)*std::pow(c,2))*(-1. - 2.*std::pow(beta,2) + 2.*std::pow(beta,4) + 2.*std::pow(beta,2)*std::pow(c,2) - 2.*std::pow(beta,4)*std::pow(c,2) + std::pow(beta,4)*std::pow(c,4))*std::pow(gs,4))/(std::pow(-1. + beta*c,2)*std::pow(1. + beta*c,2));

    return std::real(expr);

	}

}