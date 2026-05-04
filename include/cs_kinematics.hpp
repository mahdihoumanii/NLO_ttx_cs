#pragma once

#include <array>

namespace ttbar_qtsub {

struct FourVector {
    double E;
    double px;
    double py;
    double pz;
};

double dot(const FourVector& a, const FourVector& b) noexcept;

// Pairwise invariant convention used by the subtraction kernels: s_ij = 2 p_i . p_j.
double matrix_sij(const FourVector& pi, const FourVector& pj) noexcept;

// Project convention:
//   same in/out status   -> (p_i + p_j)^2
//   mixed in/out status  -> (p_i - p_j)^2
double project_sij(const FourVector& pi,
                   const FourVector& pj,
                   bool i_incoming,
                   bool j_incoming) noexcept;

// Pair tables use 1-based indexing to mirror the phase-space leg labels.
// Entry [0][*] and [*][0] are intentionally unused.
std::array<std::array<double, 5>, 5> build_matrix_sij_2to2(
    const FourVector& p1,
    const FourVector& p2,
    const FourVector& p3,
    const FourVector& p4) noexcept;

// This helper assumes the standard real-channel ordering used in the project:
//   legs 1,2 incoming; legs 3,4,5 outgoing.
// Entry [0][*] and [*][0] are intentionally unused.
std::array<std::array<double, 6>, 6> build_project_sij_2to3(
    const FourVector& p1,
    const FourVector& p2,
    const FourVector& p3,
    const FourVector& p4,
    const FourVector& p5) noexcept;

}  // namespace ttbar_qtsub