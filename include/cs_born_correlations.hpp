#pragma once

#include <array>
#include <string>

#include "cs_kinematics.hpp"

namespace ttbar_qtsub {

enum class CSChannel {
    QQbar,
    GG,
};

struct OrderedColorCorrelation {
    int emitter;
    int spectator;
    double Bij;
    double matrix_sij;
    bool emitter_massive;
    bool spectator_massive;
    std::string emitter_type;
    std::string spectator_type;
};

// Build the 12 ordered emitter-spectator pairs for 2 -> 2 Born kinematics.
// The validated physics is delegated to the existing color_correlation.cpp
// dispatchers; this layer only exposes the results in a stable ordered form.
std::array<OrderedColorCorrelation, 12> build_ordered_born_correlations(
    CSChannel channel,
    double shat,
    double mt,
    double mu,
    double gs,
    double costheta,
    const FourVector& p1,
    const FourVector& p2,
    const FourVector& p3,
    const FourVector& p4);

}  // namespace ttbar_qtsub