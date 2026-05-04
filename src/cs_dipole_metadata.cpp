// Generated from Dipoles.py — do not edit manually.
#include "cs_dipole_metadata.hpp"

#include <stdexcept>

namespace ttbar_qtsub {
namespace {

const std::vector<RealChannel> kAllRealChannels{
    RealChannel::QQbar_TTbarG,
    RealChannel::GG_TTbarG,
    RealChannel::QG_TTbarQ,
    RealChannel::QbarG_TTbarQbar,
    RealChannel::GQ_TTbarQ,
    RealChannel::GQbar_TTbarQbar,
};

const std::vector<CSDipoleMetadata> kQQbar_TTbarGDipoles{
    {DipoleClass::FF, "FF[3,5,4]", std::array<int, 2>{3, 5}, 4, "Qg", "qqbttb", 4, std::array<int, 2>{3, 5}, DenominatorKind::Pair, std::array<int, 2>{3, 5}},
    {DipoleClass::FF, "FF[4,5,3]", std::array<int, 2>{4, 5}, 3, "Qg", "qqbttb", 3, std::array<int, 2>{4, 5}, DenominatorKind::Pair, std::array<int, 2>{4, 5}},
    {DipoleClass::FI, "FI[3,5,1]", std::array<int, 2>{3, 5}, 1, "Qg", "qqbttb", 1, std::array<int, 2>{3, 5}, DenominatorKind::Pair, std::array<int, 2>{3, 5}},
    {DipoleClass::FI, "FI[3,5,2]", std::array<int, 2>{3, 5}, 2, "Qg", "qqbttb", 2, std::array<int, 2>{3, 5}, DenominatorKind::Pair, std::array<int, 2>{3, 5}},
    {DipoleClass::FI, "FI[4,5,1]", std::array<int, 2>{4, 5}, 1, "Qg", "qqbttb", 1, std::array<int, 2>{4, 5}, DenominatorKind::Pair, std::array<int, 2>{4, 5}},
    {DipoleClass::FI, "FI[4,5,2]", std::array<int, 2>{4, 5}, 2, "Qg", "qqbttb", 2, std::array<int, 2>{4, 5}, DenominatorKind::Pair, std::array<int, 2>{4, 5}},
    {DipoleClass::IF, "IF[1,5,3]", std::array<int, 2>{1, 5}, 3, "qg", "qqbttb", 3, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::IF, "IF[1,5,4]", std::array<int, 2>{1, 5}, 4, "qg", "qqbttb", 4, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::IF, "IF[2,5,3]", std::array<int, 2>{2, 5}, 3, "qbg", "qqbttb", 3, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}},
    {DipoleClass::IF, "IF[2,5,4]", std::array<int, 2>{2, 5}, 4, "qbg", "qqbttb", 4, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}},
    {DipoleClass::II, "II[1,5,2]", std::array<int, 2>{1, 5}, 2, "qg", "qqbttb", 2, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::II, "II[2,5,1]", std::array<int, 2>{2, 5}, 1, "qbg", "qqbttb", 1, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}}
};

const std::vector<CSDipoleMetadata> kGG_TTbarGDipoles{
    {DipoleClass::FF, "FF[3,5,4]", std::array<int, 2>{3, 5}, 4, "Qg", "ggttb", 4, std::array<int, 2>{3, 5}, DenominatorKind::Pair, std::array<int, 2>{3, 5}},
    {DipoleClass::FF, "FF[4,5,3]", std::array<int, 2>{4, 5}, 3, "Qg", "ggttb", 3, std::array<int, 2>{4, 5}, DenominatorKind::Pair, std::array<int, 2>{4, 5}},
    {DipoleClass::FI, "FI[3,5,1]", std::array<int, 2>{3, 5}, 1, "Qg", "ggttb", 1, std::array<int, 2>{3, 5}, DenominatorKind::Pair, std::array<int, 2>{3, 5}},
    {DipoleClass::FI, "FI[3,5,2]", std::array<int, 2>{3, 5}, 2, "Qg", "ggttb", 2, std::array<int, 2>{3, 5}, DenominatorKind::Pair, std::array<int, 2>{3, 5}},
    {DipoleClass::FI, "FI[4,5,1]", std::array<int, 2>{4, 5}, 1, "Qg", "ggttb", 1, std::array<int, 2>{4, 5}, DenominatorKind::Pair, std::array<int, 2>{4, 5}},
    {DipoleClass::FI, "FI[4,5,2]", std::array<int, 2>{4, 5}, 2, "Qg", "ggttb", 2, std::array<int, 2>{4, 5}, DenominatorKind::Pair, std::array<int, 2>{4, 5}},
    {DipoleClass::IF, "IF[1,5,3]", std::array<int, 2>{1, 5}, 3, "gg", "ggttb", 3, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::IF, "IF[1,5,4]", std::array<int, 2>{1, 5}, 4, "gg", "ggttb", 4, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::IF, "IF[2,5,3]", std::array<int, 2>{2, 5}, 3, "gg", "ggttb", 3, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}},
    {DipoleClass::IF, "IF[2,5,4]", std::array<int, 2>{2, 5}, 4, "gg", "ggttb", 4, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}},
    {DipoleClass::II, "II[1,5,2]", std::array<int, 2>{1, 5}, 2, "gg", "ggttb", 2, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::II, "II[2,5,1]", std::array<int, 2>{2, 5}, 1, "gg", "ggttb", 1, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}}
};

const std::vector<CSDipoleMetadata> kQG_TTbarQDipoles{
    {DipoleClass::IF, "IF[1,5,3]", std::array<int, 2>{1, 5}, 3, "QQ", "ggttb", 4, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::IF, "IF[1,5,4]", std::array<int, 2>{1, 5}, 4, "QQ", "ggttb", 3, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::IF, "IF[2,5,3]", std::array<int, 2>{2, 5}, 3, "gq", "qqbttb", 4, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}},
    {DipoleClass::IF, "IF[2,5,4]", std::array<int, 2>{2, 5}, 4, "gq", "qqbttb", 3, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}},
    {DipoleClass::II, "II[1,5,2]", std::array<int, 2>{1, 5}, 2, "QQ", "ggttb", 2, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::II, "II[2,5,1]", std::array<int, 2>{2, 5}, 1, "gq", "qqbttb", 1, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}}
};

const std::vector<CSDipoleMetadata> kQbarG_TTbarQbarDipoles{
    {DipoleClass::IF, "IF[1,5,3]", std::array<int, 2>{1, 5}, 3, "QBarQBar", "ggttb", 4, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::IF, "IF[1,5,4]", std::array<int, 2>{1, 5}, 4, "QBarQBar", "ggttb", 3, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::IF, "IF[2,5,3]", std::array<int, 2>{2, 5}, 3, "gqb", "qqbttb", 3, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}},
    {DipoleClass::IF, "IF[2,5,4]", std::array<int, 2>{2, 5}, 4, "gqb", "qqbttb", 4, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}},
    {DipoleClass::II, "II[1,5,2]", std::array<int, 2>{1, 5}, 2, "QBarQBar", "ggttb", 2, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::II, "II[2,5,1]", std::array<int, 2>{2, 5}, 1, "gqb", "qqbttb", 1, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}}
};

const std::vector<CSDipoleMetadata> kGQ_TTbarQDipoles{
    {DipoleClass::IF, "IF[1,5,3]", std::array<int, 2>{1, 5}, 3, "gq", "qqbttb", 4, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::IF, "IF[1,5,4]", std::array<int, 2>{1, 5}, 4, "gq", "qqbttb", 3, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::IF, "IF[2,5,3]", std::array<int, 2>{2, 5}, 3, "QQ", "ggttb", 4, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}},
    {DipoleClass::IF, "IF[2,5,4]", std::array<int, 2>{2, 5}, 4, "QQ", "ggttb", 3, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}},
    {DipoleClass::II, "II[1,5,2]", std::array<int, 2>{1, 5}, 2, "gq", "qqbttb", 2, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::II, "II[2,5,1]", std::array<int, 2>{2, 5}, 1, "QQ", "ggttb", 1, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}}
};

const std::vector<CSDipoleMetadata> kGQbar_TTbarQbarDipoles{
    {DipoleClass::IF, "IF[1,5,3]", std::array<int, 2>{1, 5}, 3, "gqb", "qqbttb", 3, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::IF, "IF[1,5,4]", std::array<int, 2>{1, 5}, 4, "gqb", "qqbttb", 4, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::IF, "IF[2,5,3]", std::array<int, 2>{2, 5}, 3, "QBarQBar", "ggttb", 4, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}},
    {DipoleClass::IF, "IF[2,5,4]", std::array<int, 2>{2, 5}, 4, "QBarQBar", "ggttb", 3, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}},
    {DipoleClass::II, "II[1,5,2]", std::array<int, 2>{1, 5}, 2, "gqb", "qqbttb", 2, std::array<int, 2>{1, 5}, DenominatorKind::Beam, std::array<int, 2>{1, 5}},
    {DipoleClass::II, "II[2,5,1]", std::array<int, 2>{2, 5}, 1, "QBarQBar", "ggttb", 1, std::array<int, 2>{2, 5}, DenominatorKind::Beam, std::array<int, 2>{2, 5}}
};

}  // namespace

const std::vector<RealChannel>& all_real_channels() {
    return kAllRealChannels;
}

std::vector<CSDipoleMetadata> get_dipoles(RealChannel channel) {
    switch (channel) {
        case RealChannel::QQbar_TTbarG: return kQQbar_TTbarGDipoles;
        case RealChannel::GG_TTbarG: return kGG_TTbarGDipoles;
        case RealChannel::QG_TTbarQ: return kQG_TTbarQDipoles;
        case RealChannel::QbarG_TTbarQbar: return kQbarG_TTbarQbarDipoles;
        case RealChannel::GQ_TTbarQ: return kGQ_TTbarQDipoles;
        case RealChannel::GQbar_TTbarQbar: return kGQbar_TTbarQbarDipoles;
    }

    throw std::invalid_argument("Unsupported RealChannel");
}

std::string to_string(RealChannel channel) {
    switch (channel) {
        case RealChannel::QQbar_TTbarG: return "qqbar -> ttbar g";
        case RealChannel::GG_TTbarG: return "gg -> ttbar g";
        case RealChannel::QG_TTbarQ: return "qg -> ttbar q";
        case RealChannel::QbarG_TTbarQbar: return "qbar g -> ttbar qbar";
        case RealChannel::GQ_TTbarQ: return "gq -> ttbar q";
        case RealChannel::GQbar_TTbarQbar: return "gqbar -> ttbar qbar";
    }

    throw std::invalid_argument("Unsupported RealChannel");
}

std::string to_string(DipoleClass dipole_class) {
    switch (dipole_class) {
        case DipoleClass::FF: return "FF";
        case DipoleClass::FI: return "FI";
        case DipoleClass::IF: return "IF";
        case DipoleClass::II: return "II";
    }

    throw std::invalid_argument("Unsupported DipoleClass");
}

std::string to_string(DenominatorKind denominator_kind) {
    switch (denominator_kind) {
        case DenominatorKind::Pair: return "pair";
        case DenominatorKind::Beam: return "beam";
    }

    throw std::invalid_argument("Unsupported DenominatorKind");
}

}  // namespace ttbar_qtsub
