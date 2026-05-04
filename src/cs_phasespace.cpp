#include "cs_phasespace.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace cs_phasespace {
namespace {

double parse_double(const std::string& text, const std::string& label) {
    try {
        std::size_t parsed = 0;
        const double value = std::stod(text, &parsed);
        if (parsed != text.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return value;
    } catch (const std::exception&) {
        throw std::runtime_error("Failed to parse " + label + ": " + text);
    }
}

int parse_int(const std::string& text, const std::string& label) {
    try {
        std::size_t parsed = 0;
        const int value = std::stoi(text, &parsed);
        if (parsed != text.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return value;
    } catch (const std::exception&) {
        throw std::runtime_error("Failed to parse " + label + ": " + text);
    }
}

double safe_sqrt(double value) {
    return std::sqrt(std::max(0.0, value));
}

int bit_for_leg(int leg) {
    return 1 << (leg - 1);
}

double mass_for_bit(int bit, double top_mass) {
    if (bit == 4 || bit == 8) {
        return top_mass;
    }
    return 0.0;
}

double sqrt_min_for_mask(int mask, double top_mass) {
    double result = 0.0;
    for (int bit = 1; bit <= 16; bit <<= 1) {
        if ((mask & bit) != 0) {
            result += mass_for_bit(bit, top_mass);
        }
    }
    return result;
}

std::vector<int> parse_int_list(const std::string& text) {
    std::vector<int> values;
    for (const std::string& token : whitespace_tokens(text)) {
        values.push_back(parse_int(token, "phase-space integer"));
    }
    return values;
}

ChannelDefinition parse_channel_content(const std::string& content) {
    std::vector<std::string> groups = split(content, '-');
    while (groups.size() < 4) {
        groups.emplace_back();
    }
    ChannelDefinition channel;
    channel.p = parse_int_list(groups[0]);
    channel.f = parse_int_list(groups[1]);
    channel.t = parse_int_list(groups[2]);
    channel.d = parse_int_list(groups[3]);
    return channel;
}

FourVector decay_momentum(const FourVector& parent,
                          double parent_s,
                          double child1_s,
                          double child2_s,
                          double r1,
                          double r2,
                          bool rotate_system) {
    const double sqrt_parent_s = std::sqrt(parent_s);
    const double sqrt_lambda = std::sqrt(std::max(0.0, lambda(parent_s, child1_s, child2_s)));
    FourVector k10{(parent_s + child1_s - child2_s) / (2.0 * sqrt_parent_s), 0.0, 0.0, sqrt_lambda / (2.0 * sqrt_parent_s)};
    const double costheta = 2.0 * r1 - 1.0;
    const double phi = 2.0 * kPi * r2;
    if (rotate_system) {
        return k10.rotation(phi, costheta).rotation(parent).boost(parent.pinv(), parent_s);
    }
    return k10.rotation(phi, costheta).boost(parent.pinv(), parent_s);
}

}  // namespace

double FourVector::r2() const { return x * x + y * y + z * z; }

double FourVector::r() const { return std::sqrt(r2()); }

double FourVector::m2() const { return e * e - r2(); }

double FourVector::m() const { return std::sqrt(std::max(0.0, m2())); }

double FourVector::sp3(const FourVector& other) const { return x * other.x + y * other.y + z * other.z; }

FourVector FourVector::pinv() const { return {e, -x, -y, -z}; }

FourVector FourVector::operator+(const FourVector& other) const { return {e + other.e, x + other.x, y + other.y, z + other.z}; }

FourVector FourVector::operator-(const FourVector& other) const { return {e - other.e, x - other.x, y - other.y, z - other.z}; }

FourVector FourVector::operator-() const { return {-e, -x, -y, -z}; }

FourVector FourVector::operator*(double factor) const { return {factor * e, factor * x, factor * y, factor * z}; }

FourVector FourVector::operator/(double factor) const { return {e / factor, x / factor, y / factor, z / factor}; }

FourVector FourVector::rotation(double phi, double costheta) const {
    const double sintheta = std::sqrt(std::max(0.0, 1.0 - costheta * costheta));
    const double sinphi = std::sin(phi);
    const double cosphi = std::cos(phi);
    return {
        e,
        cosphi * costheta * x - sinphi * y + cosphi * sintheta * z,
        sinphi * costheta * x + cosphi * y + sinphi * sintheta * z,
        -sintheta * x + costheta * z,
    };
}

FourVector FourVector::rotation_beam(double phi, double costheta) const {
    const double sintheta = std::sqrt(std::max(0.0, 1.0 - costheta * costheta));
    return {e, std::cos(phi) * sintheta * z, std::sin(phi) * sintheta * z, costheta * z};
}

FourVector FourVector::rotation(const FourVector& axis) const {
    double phi = 0.0;
    double costheta = 0.0;
    if (axis.r() > 0.0) {
        costheta = axis.z / axis.r();
        if (axis.x > 0.0) {
            phi = std::atan(axis.y / axis.x);
        } else if (axis.x < 0.0) {
            phi = std::atan(axis.y / axis.x) + kPi;
        } else if (axis.y > 0.0) {
            phi = kPi / 2.0;
        } else if (axis.y < 0.0) {
            phi = 3.0 * kPi / 2.0;
        }
    }
    return rotation(phi, costheta);
}

FourVector FourVector::boost(const FourVector& frame, double frame_m2) const {
    const double frame_m = std::sqrt(frame_m2);
    const double spatial_dot = sp3(frame);
    const double common = spatial_dot / (frame_m + frame.e) - e;
    return {
        (e * frame.e - spatial_dot) / frame_m,
        x + frame.x / frame_m * common,
        y + frame.y / frame_m * common,
        z + frame.z / frame_m * common,
    };
}

FourVector FourVector::boost(const FourVector& frame) const {
    return boost(frame, frame.m2());
}

FourVector operator*(double factor, const FourVector& vector) {
    return vector * factor;
}

std::string trim(const std::string& value) {
    const std::size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const std::size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::vector<std::string> split(const std::string& text, char delimiter) {
    std::vector<std::string> pieces;
    std::string piece;
    std::stringstream stream(text);
    while (std::getline(stream, piece, delimiter)) {
        pieces.push_back(piece);
    }
    if (!text.empty() && text.back() == delimiter) {
        pieces.emplace_back();
    }
    return pieces;
}

std::vector<std::string> whitespace_tokens(const std::string& text) {
    std::vector<std::string> tokens;
    std::stringstream stream(text);
    std::string token;
    while (stream >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

PhaseSpaceBlock parse_phase_space_block(const std::string& path, const std::string& identifier) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Failed to open phase-space data: " + path);
    }

    PhaseSpaceBlock block;
    bool active = false;
    std::string line;
    while (std::getline(input, line)) {
        if (trim(line) == identifier) {
            active = true;
            continue;
        }
        if (active && trim(line).empty()) {
            break;
        }
        if (!active || line.find('|') == std::string::npos) {
            continue;
        }

        std::vector<std::string> pieces = split(line, '|');
        if (pieces.empty()) {
            continue;
        }
        const std::vector<std::string> head = whitespace_tokens(pieces[0]);
        if (head.size() < 2) {
            continue;
        }
        const std::string name = head[0];
        const int count = parse_int(head[1], name + " count");
        if (static_cast<int>(pieces.size()) < count + 1) {
            throw std::runtime_error("Malformed phase-space line for " + name);
        }

        if (name == "smin") {
            block.smin.clear();
            for (int i = 0; i < count; ++i) {
                block.smin.push_back(parse_int_list(pieces[i + 1]));
            }
        } else if (name == "smax") {
            block.smax.clear();
            for (int i = 0; i < count; ++i) {
                const std::vector<int> values = parse_int_list(pieces[i + 1]);
                if (values.size() < 2) {
                    throw std::runtime_error("Malformed smax entry");
                }
                SmaxDefinition smax;
                smax.process = values[0];
                smax.out.assign(values.begin() + 1, values.end());
                block.smax.push_back(smax);
            }
        } else if (name == "p") {
            block.p.clear();
            for (int i = 0; i < count; ++i) {
                const std::vector<int> values = parse_int_list(pieces[i + 1]);
                if (values.size() != 4) {
                    throw std::runtime_error("Malformed propagator entry");
                }
                block.p.push_back({values[0], values[1], values[2], values[3]});
            }
        } else if (name == "f") {
            block.f.clear();
            for (int i = 0; i < count; ++i) {
                const std::vector<int> values = parse_int_list(pieces[i + 1]);
                if (values.size() != 3) {
                    throw std::runtime_error("Malformed timelike-invariant entry");
                }
                block.f.push_back({values[0], values[1], values[2]});
            }
        } else if (name == "t") {
            block.t.clear();
            for (int i = 0; i < count; ++i) {
                const std::vector<int> values = parse_int_list(pieces[i + 1]);
                if (values.size() != 6) {
                    throw std::runtime_error("Malformed t-channel entry");
                }
                block.t.push_back({values[0], values[1], values[2], values[3], values[4], values[5]});
            }
        } else if (name == "d") {
            block.d.clear();
            for (int i = 0; i < count; ++i) {
                const std::vector<int> values = parse_int_list(pieces[i + 1]);
                if (values.size() != 3) {
                    throw std::runtime_error("Malformed decay entry");
                }
                block.d.push_back({values[0], values[1], values[2]});
            }
        } else if (name == "channel") {
            block.channels.clear();
            for (int i = 0; i < count; ++i) {
                block.channels.push_back(parse_channel_content(pieces[i + 1]));
            }
        }
    }

    if (block.channels.empty()) {
        throw std::runtime_error("Did not find channel table for " + identifier);
    }
    return block;
}

double lambda(double x, double y, double z) {
    if (x == 0.0) {
        return std::pow(y - z, 2);
    }
    if (y == 0.0) {
        return std::pow(x - z, 2);
    }
    if (z == 0.0) {
        return std::pow(y - x, 2);
    }
    return std::pow(x - y - z, 2) - 4.0 * y * z;
}

double h_propto_pot(double r, double smin, double smax, double exponent, double cut_technical) {
    if (exponent != 1.0) {
        return std::pow(r * std::pow(smax + cut_technical, 1.0 - exponent)
                            + (1.0 - r) * std::pow(smin + cut_technical, 1.0 - exponent),
                        1.0 / (1.0 - exponent))
             - cut_technical;
    }
    return std::exp(r * std::log(smax + cut_technical) + (1.0 - r) * std::log(smin + cut_technical)) - cut_technical;
}

double g_propto_pot(double s, double smin, double smax, double exponent, double cut_technical) {
    if (exponent != 1.0) {
        return (1.0 - exponent) / ((std::pow(smax + cut_technical, 1.0 - exponent) - std::pow(smin + cut_technical, 1.0 - exponent)) * std::pow(s + cut_technical, exponent));
    }
    return 1.0 / ((std::log(smax + cut_technical) - std::log(smin + cut_technical)) * (s + cut_technical));
}

double h_propto_pot_mod(double r, double smin, double smax, double exponent, double cut_technical) {
    if (exponent != 1.0) {
        return std::pow(r * std::pow(smax + cut_technical, 1.0 - exponent)
                            + (1.0 - r) * std::pow(smin + cut_technical, 1.0 - exponent),
                        1.0 / (1.0 - exponent))
             - cut_technical;
    }
    return smin + (smax - smin) * (std::exp(r * std::log(1.0 + cut_technical) + (1.0 - r) * std::log(cut_technical)) - cut_technical);
}

double g_propto_pot_mod(double s, double smin, double smax, double exponent, double cut_technical) {
    if (exponent != 1.0) {
        return (1.0 - exponent) / ((std::pow(smax + cut_technical, 1.0 - exponent) - std::pow(smin + cut_technical, 1.0 - exponent)) * std::pow(s + cut_technical, exponent));
    }
    return 1.0 / ((std::log(1.0 + cut_technical) - std::log(cut_technical)) * ((s - smin) + cut_technical * (smax - smin)));
}

double c_phi(double r) {
    return 2.0 * kPi * r;
}

double smin_value(const PhaseSpaceBlock& phase_space, int index, double top_mass) {
    const std::vector<int>& entry = phase_space.smin.at(index);
    int mask = 0;
    for (int value : entry) {
        mask += value;
    }
    return std::pow(sqrt_min_for_mask(mask, top_mass), 2);
}

double smax_value(const PhaseSpaceBlock& phase_space, int index, double sqrt_shat, double top_mass) {
    const SmaxDefinition& entry = phase_space.smax.at(index);
    if (entry.process != 0) {
        throw std::runtime_error("Only process-0 smax is implemented in cs_phasespace");
    }
    int mask = 0;
    for (int value : entry.out) {
        mask += value;
    }
    return std::pow(sqrt_shat - sqrt_min_for_mask(mask, top_mass), 2);
}

double mapped_vanishing_width_s(double r, double smin, double smax, double mass2_mapping, double exponent) {
    if (exponent == 1.0) {
        const double low = std::log(smin - mass2_mapping);
        const double high = std::log(smax - mass2_mapping);
        return std::exp(r * high + (1.0 - r) * low) + mass2_mapping;
    }
    const double low = std::pow(smin - mass2_mapping, 1.0 - exponent);
    const double high = std::pow(smax - mass2_mapping, 1.0 - exponent);
    return std::pow(r * high + (1.0 - r) * low, 1.0 / (1.0 - exponent)) + mass2_mapping;
}

double mapped_vanishing_width_g(double s, double smin, double smax, double mass2_mapping, double exponent) {
    if (exponent == 1.0) {
        const double low = std::log(smin - mass2_mapping);
        const double high = std::log(smax - mass2_mapping);
        return 1.0 / ((high - low) * (s - mass2_mapping));
    }
    const double low = std::pow(smin - mass2_mapping, 1.0 - exponent);
    const double high = std::pow(smax - mass2_mapping, 1.0 - exponent);
    return (1.0 - exponent) / ((high - low) * std::pow(s - mass2_mapping, exponent));
}

double mapped_tchannel_s(double r, double tmin, double tmax, double mass2_mapping, double exponent) {
    if (exponent == 1.0) {
        const double low = std::log((-tmin) - mass2_mapping);
        const double high = std::log((-tmax) - mass2_mapping);
        return -std::exp(r * high + (1.0 - r) * low) - mass2_mapping;
    }
    const double low = std::pow((-tmin) - mass2_mapping, 1.0 - exponent);
    const double high = std::pow((-tmax) - mass2_mapping, 1.0 - exponent);
    return -std::pow(r * high + (1.0 - r) * low, 1.0 / (1.0 - exponent)) - mass2_mapping;
}

double mapped_tchannel_g(double s_in1_out1, double tmin, double tmax, double mass2_mapping, double exponent) {
    if (exponent == 1.0) {
        return -1.0 / ((std::log((-tmax) - mass2_mapping) - std::log((-tmin) - mass2_mapping)) * ((-s_in1_out1) - mass2_mapping));
    }
    const double low = std::pow((-tmin) - mass2_mapping, 1.0 - exponent);
    const double high = std::pow((-tmax) - mass2_mapping, 1.0 - exponent);
    return -(1.0 - exponent) / ((high - low) * std::pow((-s_in1_out1) - mass2_mapping, exponent));
}

double timelike_invariant_density(const PhaseSpaceBlock& phase_space, int index, double sqrt_shat, double top_mass) {
    const TimelikeInvariantDefinition& invariant = phase_space.f.at(index);
    const double f_smin = smin_value(phase_space, invariant.smin, top_mass);
    const double f_smax = smax_value(phase_space, invariant.smax, sqrt_shat, top_mass);
    return 1.0 / (f_smax - f_smin);
}

double decay_density(double parent_s, double child1_s, double child2_s, double gis1, double gis2) {
    const double sqrt_lambda = std::sqrt(std::max(0.0, lambda(parent_s, child1_s, child2_s)));
    return 0.5 * (4.0 * parent_s) / (kPi * sqrt_lambda) * gis1 * gis2;
}

double tchannel_mass2_mapping(int mass_id, double top_mass, double mass0) {
    const double mass2 = (mass_id == 6) ? top_mass * top_mass : 0.0;
    return (mass2 != 0.0) ? -mass2 : mass0;
}

double apply_propagator(const PhaseSpaceBlock& phase_space,
                        const std::map<int, TaggedBlock>& p_blocks,
                        int index,
                        double sqrt_shat,
                        double top_mass,
                        double top_mass2,
                        double mass0,
                        double nuxs,
                        std::map<int, double>& local_s) {
    const PropagatorDefinition& propagator = phase_space.p.at(index);
    const TaggedBlock& block = p_blocks.at(index);
    const double p_smin = block.has_smin ? block.smin : smin_value(phase_space, propagator.smin, top_mass);
    const double p_smax = block.has_smax ? block.smax : smax_value(phase_space, propagator.smax, sqrt_shat, top_mass);
    const double mass2_mapping = ((propagator.mass_id == 6) ? top_mass2 : 0.0) + mass0;
    local_s[propagator.out] = block.has_s ? block.s : mapped_vanishing_width_s(block.r1, p_smin, p_smax, mass2_mapping, nuxs);
    return mapped_vanishing_width_g(local_s[propagator.out], p_smin, p_smax, mass2_mapping, nuxs) * p_blocks.at(index).gis1;
}

double apply_timelike_invariant(const PhaseSpaceBlock& phase_space,
                                const std::map<int, TaggedBlock>& f_blocks,
                                int index,
                                double sqrt_shat,
                                double top_mass,
                                std::map<int, double>& local_s) {
    const TimelikeInvariantDefinition& invariant = phase_space.f.at(index);
    const TaggedBlock& block = f_blocks.at(index);
    const double f_smin = block.has_smin ? block.smin : smin_value(phase_space, invariant.smin, top_mass);
    const double f_smax = block.has_smax ? block.smax : smax_value(phase_space, invariant.smax, sqrt_shat, top_mass);
    local_s[invariant.out] = block.has_s ? block.s : block.r1 * f_smax + (1.0 - block.r1) * f_smin;
    return block.gis1 / (f_smax - f_smin);
}

double apply_tchannel(const PhaseSpaceBlock& phase_space,
                      const std::map<int, TaggedBlock>& t_blocks,
                      int index,
                      double top_mass,
                      double mass0,
                      double nuxt,
                      std::map<int, double>& local_s,
                      std::map<int, FourVector>& local_p) {
    const TChannelDefinition& tchannel = phase_space.t.at(index);
    const double s_out = local_s.at(tchannel.out);
    const double s_out1 = local_s.at(tchannel.out1);
    const double s_out2 = local_s.at(tchannel.out2);
    const double s_in1 = local_s.at(tchannel.in1);
    const double s_in2 = local_s.at(tchannel.in2);

    const double term1 = (s_out + s_out1 - s_out2) * (s_out + s_in1 - s_in2);
    const double term2in = std::sqrt(std::max(0.0, lambda(s_out, s_in1, s_in2)));
    const double term2out = std::sqrt(std::max(0.0, lambda(s_out, s_out1, s_out2)));
    const double term2 = term2in * term2out;
    const double term3 = 2.0 * std::sqrt(s_out);
    double tmin = s_out1 + s_in1 - (term1 + term2) / (2.0 * s_out);
    double tmax = s_out1 + s_in1 - (term1 - term2) / (2.0 * s_out);
    if (tmax > 0.0) {
        tmax = 0.0;
    }

    const double mass2_mapping = tchannel_mass2_mapping(tchannel.mass_id, top_mass, mass0);
    const int in1_out1 = tchannel.in1 + tchannel.out1;
    local_s[in1_out1] = mapped_tchannel_s(t_blocks.at(index).r1, tmin, tmax, mass2_mapping, nuxt);

    const FourVector p10{(s_out + s_in1 - s_in2) / term3, 0.0, 0.0, term2in / term3};
    const FourVector p20{(s_out - s_in1 + s_in2) / term3, 0.0, 0.0, -term2in / term3};
    FourVector k10{(s_out + s_out1 - s_out2) / term3, 0.0, 0.0, term2out / term3};
    const FourVector p1b = local_p.at(tchannel.in1).boost(local_p.at(tchannel.out), s_out);
    const double phi = 2.0 * kPi * t_blocks.at(index).r2;
    const double costheta = ((2.0 * s_out * (local_s.at(in1_out1) - s_out1 - s_in1)) + term1) / term2;
    k10 = k10.rotation_beam(phi, costheta);
    const FourVector k20 = p10 + p20 - k10;

    local_p[tchannel.out1] = k10.rotation(p1b).boost(local_p.at(tchannel.out).pinv(), s_out);
    local_p[tchannel.out2] = k20.rotation(p1b).boost(local_p.at(tchannel.out).pinv(), s_out);
    local_p[in1_out1] = local_p.at(tchannel.in1) - local_p.at(tchannel.out1);

    return mapped_tchannel_g(local_s.at(in1_out1), tmin, tmax, mass2_mapping, nuxt) * (2.0 * term2in) / kPi * t_blocks.at(index).gis1 * t_blocks.at(index).gis2;
}

double apply_decay(const PhaseSpaceBlock& phase_space,
                   const std::map<int, TaggedBlock>& d_blocks,
                   int index,
                   bool switch_decay_system,
                   int unrotated_all_final_bit,
                   std::map<int, double>& local_s,
                   std::map<int, FourVector>& local_p) {
    const DecayDefinition& decay = phase_space.d.at(index);
    const TaggedBlock& block = d_blocks.at(index);
    const bool rotate_system = switch_decay_system && decay.out != unrotated_all_final_bit;
    local_p[decay.out1] = decay_momentum(local_p.at(decay.out), local_s.at(decay.out), local_s.at(decay.out1), local_s.at(decay.out2), block.r1, block.r2, rotate_system);
    local_p[decay.out2] = local_p.at(decay.out) - local_p.at(decay.out1);
    const double density_s = block.has_s ? block.s : local_s.at(decay.out);
    const double density_s1 = block.has_s1 ? block.s1 : local_s.at(decay.out1);
    const double density_s2 = block.has_s2 ? block.s2 : local_s.at(decay.out2);
    const double density_sqrt_lambda = block.has_sqrt_lambda ? block.sqrt_lambda : std::sqrt(std::max(0.0, lambda(density_s, density_s1, density_s2)));
    const double density_gcostheta = block.has_gcostheta ? block.gcostheta : 0.5;
    return density_gcostheta * (4.0 * density_s) / (kPi * density_sqrt_lambda) * block.gis1 * block.gis2;
}

void fill_born_invariants(BranchState& state) {
    state.p[0] = state.p.at(1) + state.p.at(2);
    state.s[0] = state.p.at(0).m2();
    state.p[12] = state.p.at(4) + state.p.at(8);
    state.s[12] = state.p.at(12).m2();
    state.p[5] = state.p.at(1) - state.p.at(4);
    state.p[9] = state.p.at(1) - state.p.at(8);
    state.p[6] = state.p.at(2) - state.p.at(4);
    state.p[10] = state.p.at(2) - state.p.at(8);
    state.s[5] = state.p.at(5).m2();
    state.s[9] = state.p.at(9).m2();
    state.s[6] = state.p.at(6).m2();
    state.s[10] = state.p.at(10).m2();
}

}  // namespace cs_phasespace