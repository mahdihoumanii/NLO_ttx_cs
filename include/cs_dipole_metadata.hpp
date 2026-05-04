// Generated from Dipoles.py — do not edit manually.
#pragma once

#include <array>
#include <string>
#include <vector>

namespace ttbar_qtsub {

enum class RealChannel {
    QQbar_TTbarG,
    GG_TTbarG,
    QG_TTbarQ,
    QbarG_TTbarQbar,
    GQ_TTbarQ,
    GQbar_TTbarQbar,
};

enum class DipoleClass {
    FF,
    FI,
    IF,
    II,
};

enum class DenominatorKind {
    Pair,
    Beam,
};

struct CSDipoleMetadata {
    DipoleClass dipole_class;
    std::string symbol;
    std::array<int, 2> emitter_legs;
    int spectator_leg;
    std::string kernel;
    std::string born_channel;
    int corr_left;
    std::array<int, 2> corr_right;
    DenominatorKind denominator_kind;
    std::array<int, 2> denominator_legs;
};

const std::vector<RealChannel>& all_real_channels();
std::vector<CSDipoleMetadata> get_dipoles(RealChannel channel);
std::string to_string(RealChannel channel);
std::string to_string(DipoleClass dipole_class);
std::string to_string(DenominatorKind denominator_kind);

}  // namespace ttbar_qtsub
