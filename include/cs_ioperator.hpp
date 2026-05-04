#pragma once

#include <vector>

#include "cs_born_correlations.hpp"

namespace ttbar_qtsub {

struct CDSTPairCoefficients {
    int emitter;
    int spectator;
    double Bij;
    double matrix_sij;

    double VS0;
    double VS1;
    double VS2;
    double VNS;
    double finite0;
    double finite1;

    double log_mu2_over_sij;
    double log2_half;

    double bracket_finite;
    double contribution_finite;
};

struct CSIOperatorResult {
    double finite;
    std::vector<CDSTPairCoefficients> pairs;
};

// Assemble the finite CDST I-operator contribution for the validated Born
// channels using the ordered B_ij interface. This function keeps every pair's
// intermediate coefficients exposed in the returned diagnostics.
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
    const FourVector& p4);

}  // namespace ttbar_qtsub