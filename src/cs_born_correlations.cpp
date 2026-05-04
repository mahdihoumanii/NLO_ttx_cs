#include "cs_born_correlations.hpp"

#include <array>
#include <stdexcept>

#include "color_correlation.hpp"

namespace ttbar_qtsub {
namespace {

struct LegMetadata {
    bool massive;
    const char* type;
};

LegMetadata metadata_for_leg(CSChannel channel, int leg) {
    switch (channel) {
        case CSChannel::QQbar:
            switch (leg) {
                case 1:
                case 2:
                    return {false, "quark"};
                case 3:
                case 4:
                    return {true, "massive quark"};
                default:
                    break;
            }
            break;
        case CSChannel::GG:
            switch (leg) {
                case 1:
                case 2:
                    return {false, "gluon"};
                case 3:
                case 4:
                    return {true, "massive quark"};
                default:
                    break;
            }
            break;
    }

    throw std::invalid_argument("metadata_for_leg: leg must be in 1..4");
}

double bij_for_pair(CSChannel channel,
                    double shat,
                    double mt,
                    double mu,
                    double gs,
                    double costheta,
                    int emitter,
                    int spectator) {
    switch (channel) {
        case CSChannel::QQbar:
            return color_correlation_qq(shat, mt, mu, gs, costheta, emitter, spectator);
        case CSChannel::GG:
            return color_correlation_gg(shat, mt, mu, gs, costheta, emitter, spectator);
    }

    throw std::invalid_argument("bij_for_pair: unsupported channel");
}

}  // namespace

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
    const FourVector& p4) {
    const std::array<std::array<int, 2>, 12> ordered_pairs{{
        {{1, 2}}, {{1, 3}}, {{1, 4}},
        {{2, 1}}, {{2, 3}}, {{2, 4}},
        {{3, 1}}, {{3, 2}}, {{3, 4}},
        {{4, 1}}, {{4, 2}}, {{4, 3}},
    }};

    const std::array<FourVector, 5> legs{{{}, p1, p2, p3, p4}};
    std::array<OrderedColorCorrelation, 12> correlations{};

    for (std::size_t index = 0; index < ordered_pairs.size(); ++index) {
        const int emitter = ordered_pairs[index][0];
        const int spectator = ordered_pairs[index][1];
        const LegMetadata emitter_meta = metadata_for_leg(channel, emitter);
        const LegMetadata spectator_meta = metadata_for_leg(channel, spectator);

        correlations[index] = OrderedColorCorrelation{
            emitter,
            spectator,
            bij_for_pair(channel, shat, mt, mu, gs, costheta, emitter, spectator),
            matrix_sij(legs[emitter], legs[spectator]),
            emitter_meta.massive,
            spectator_meta.massive,
            emitter_meta.type,
            spectator_meta.type,
        };
    }

    return correlations;
}

}  // namespace ttbar_qtsub