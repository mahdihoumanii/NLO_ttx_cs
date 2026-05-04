#include "cs_kinematics.hpp"

#include <array>

namespace ttbar_qtsub {
namespace {

FourVector add(const FourVector& a, const FourVector& b) noexcept {
    return {a.E + b.E, a.px + b.px, a.py + b.py, a.pz + b.pz};
}

FourVector subtract(const FourVector& a, const FourVector& b) noexcept {
    return {a.E - b.E, a.px - b.px, a.py - b.py, a.pz - b.pz};
}

double mass_sq(const FourVector& p) noexcept {
    return dot(p, p);
}

}  // namespace

double dot(const FourVector& a, const FourVector& b) noexcept {
    return a.E * b.E - a.px * b.px - a.py * b.py - a.pz * b.pz;
}

double matrix_sij(const FourVector& pi, const FourVector& pj) noexcept {
    return 2.0 * dot(pi, pj);
}

double project_sij(const FourVector& pi,
                   const FourVector& pj,
                   bool i_incoming,
                   bool j_incoming) noexcept {
    // The sign conversion is explicit here by construction.
    // Same side of the scattering process  -> combine momenta.
    // Opposite sides                        -> take their difference.
    if (i_incoming == j_incoming) {
        return mass_sq(add(pi, pj));
    }
    return mass_sq(subtract(pi, pj));
}

std::array<std::array<double, 5>, 5> build_matrix_sij_2to2(
    const FourVector& p1,
    const FourVector& p2,
    const FourVector& p3,
    const FourVector& p4) noexcept {
    std::array<std::array<double, 5>, 5> sij{};
    const std::array<FourVector, 5> legs{{{}, p1, p2, p3, p4}};

    for (int i = 1; i <= 4; ++i) {
        for (int j = i + 1; j <= 4; ++j) {
            const double value = matrix_sij(legs[i], legs[j]);
            sij[i][j] = value;
            sij[j][i] = value;
        }
    }

    return sij;
}

std::array<std::array<double, 6>, 6> build_project_sij_2to3(
    const FourVector& p1,
    const FourVector& p2,
    const FourVector& p3,
    const FourVector& p4,
    const FourVector& p5) noexcept {
    std::array<std::array<double, 6>, 6> sij{};
    const std::array<FourVector, 6> legs{{{}, p1, p2, p3, p4, p5}};
    const std::array<bool, 6> incoming{{false, true, true, false, false, false}};

    for (int i = 1; i <= 5; ++i) {
        for (int j = i + 1; j <= 5; ++j) {
            const double value = project_sij(legs[i], legs[j], incoming[i], incoming[j]);
            sij[i][j] = value;
            sij[j][i] = value;
        }
    }

    return sij;
}

}  // namespace ttbar_qtsub