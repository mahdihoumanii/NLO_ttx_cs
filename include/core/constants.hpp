#pragma once
// =============================================================================
// core/constants.hpp
// Central numeric constants for the ttbar_qtsub project.
//
// All other headers and source files should include THIS header (directly or
// via common.hpp) instead of defining their own pi, tolerances, or color
// factors.  The old include/common.hpp now re-exports this header for backward
// compatibility.
// =============================================================================

#include <cmath>

namespace ttbar_qtsub {

// -----------------------------------------------------------------------------
// Mathematical constants
// -----------------------------------------------------------------------------

/// pi, to full double precision.
constexpr double kPi   = 3.141592653589793238462643383279502884;
constexpr double kPiSq = kPi * kPi;   // pi^2, appears in Li2/Li3 expressions.

// -----------------------------------------------------------------------------
// Physical unit conversion
// -----------------------------------------------------------------------------

/// 1 GeV^{-2} in picobarns  (hbar^2 c^2).
constexpr double kPbPerGeV2 = 0.389379338e9;

// -----------------------------------------------------------------------------
// SU(3) color group constants
// -----------------------------------------------------------------------------

constexpr int    kNc = 3;              ///< Number of colors.
constexpr double kCF = 4.0 / 3.0;     ///< Quark Casimir  C_F = (N_c^2-1)/(2 N_c).
constexpr double kCA = 3.0;           ///< Gluon Casimir  C_A = N_c.
constexpr double kTR = 0.5;           ///< Fundamental normalization  T_R = 1/2.

// -----------------------------------------------------------------------------
// Default physics parameters  (used as RunConfig defaults; NOT hardwired in
// any amplitude expression —  those always receive explicit arguments)
// -----------------------------------------------------------------------------

constexpr double kMtDefault  = 173.0; ///< Default top pole mass  [GeV].
constexpr int    kNfDefault  = 5;     ///< Default number of light quark flavors.
constexpr double kSqrtSDefault = 13000.0; ///< Default collider sqrt(S)  [GeV].

// -----------------------------------------------------------------------------
// Numerical tolerances
// -----------------------------------------------------------------------------

/// A shat value is considered "at or below threshold" when
///    shat  <=  4 mt^2 * (1 + kEpsilonKin).
/// Keeping this at exactly 0 means the check is shat > 4 mt^2 strictly.
constexpr double kEpsilonKin      = 0.0;

/// Absolute guard on  |costheta| <= 1.
constexpr double kEpsilonCosTheta = 1.0e-12;

// -----------------------------------------------------------------------------
// Convenience helpers
// -----------------------------------------------------------------------------

/// Square of x.  Avoids std::pow(x, 2) in hot paths.
inline double sqr(double x) noexcept { return x * x; }

}  // namespace ttbar_qtsub
