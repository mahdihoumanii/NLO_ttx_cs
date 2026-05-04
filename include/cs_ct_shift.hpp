#pragma once

#include <vector>

#include "cs_born_correlations.hpp"
#include "cs_kinematics.hpp"

namespace ttbar_qtsub {

struct CTShiftResult {
    double total;
    double decoupling;
};

CTShiftResult compute_ct_shift_decoupling_approx(
    CSChannel channel,
    double Born,
    double alpha_s,
    double muR,
    double mt);

CTShiftResult compute_ct_shift_local_exact(
    CSChannel channel,
    const std::vector<FourVector>& momenta,
    double alpha_s,
    double muR,
    double mt);

// Backward-compatible alias for the standalone approximate decoupling shift.
CTShiftResult compute_ct_shift(
    CSChannel channel,
    double Born,
    double alpha_s,
    double muR,
    double mt);

}  // namespace ttbar_qtsub