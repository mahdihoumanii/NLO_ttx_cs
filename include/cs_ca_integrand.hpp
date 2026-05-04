#pragma once

#include <array>
#include <string>
#include <vector>

#include "cs_born_correlations.hpp"
#include "cs_kinematics.hpp"

namespace ttbar_qtsub {

enum class CABornSubprocess {
    GG_TTbar,
    DDbar_TTbar,
    UUbar_TTbar,
    BBbar_TTbar,
};

struct CAPDFCombination {
    int direction = 1;
    int flav1 = 0;
    int flav2 = 0;
};

struct CAInput {
    CABornSubprocess subprocess = CABornSubprocess::UUbar_TTbar;
    std::vector<FourVector> momenta;
    double x1 = 0.0;
    double x2 = 0.0;
    double z1 = 1.0;
    double z2 = 1.0;
    double alpha_s = 0.0;
    double muF = 0.0;
    double muR = 0.0;
    double mt = 0.0;
    double g_z1 = 1.0;
    double g_z2 = 1.0;
};

struct CACollinearContribution {
    std::string name;
    int type = 0;
    int emitter = 0;
    std::array<int, 3> in_collinear{};
    std::vector<int> spectators;
    std::vector<double> me2_cf;
    std::array<double, 3> sum_K{};
    std::array<double, 3> sum_P{};
    std::array<double, 3> kp{};
    std::array<double, 3> pdf_factor_tsv{};
    std::array<double, 3> pdf_factor_cv{};
    double integrand_tsv = 0.0;
    double integrand_cv = 0.0;
};

struct CAResult {
    std::vector<CACollinearContribution> collinear;
    double sum_integrand_tsv = 0.0;
    double sum_integrand_cv = 0.0;
    double born_me2 = 0.0;
    double shat = 0.0;
    double costheta = 0.0;
};

CSChannel ca_channel(CABornSubprocess subprocess);
const char* ca_subprocess_name(CABornSubprocess subprocess);
std::vector<CAPDFCombination> ca_born_pdf_combinations(CABornSubprocess subprocess);

CAResult compute_ca_integrand(const CAInput& input);

}  // namespace ttbar_qtsub