#pragma once

#include <vector>

#include "cs_kinematics.hpp"
#include "cs_born_correlations.hpp"

namespace ttbar_qtsub {

struct VAKernelResult {
    double V;
    double X;
    double I;

    double V_plus_I;
    double VA_total;
};

// Assemble the Born-kinematics VA kernel from the already validated finite
// virtual, standalone CT shift, and finite CDST I-operator pieces.
// The input momenta are expected in the standard project order p1, p2, p3, p4.
VAKernelResult compute_va_kernel(
    CSChannel channel,
    const std::vector<FourVector>& momenta,
    double alpha_s,
    double muR,
    double mt);

}  // namespace ttbar_qtsub