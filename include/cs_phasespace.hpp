#pragma once

#include <map>
#include <string>
#include <vector>

namespace cs_phasespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct FourVector {
    double e = 0.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    double r2() const;
    double r() const;
    double m2() const;
    double m() const;
    double sp3(const FourVector& other) const;
    FourVector pinv() const;

    FourVector operator+(const FourVector& other) const;
    FourVector operator-(const FourVector& other) const;
    FourVector operator-() const;
    FourVector operator*(double factor) const;
    FourVector operator/(double factor) const;

    FourVector rotation(double phi, double costheta) const;
    FourVector rotation_beam(double phi, double costheta) const;
    FourVector rotation(const FourVector& axis) const;
    FourVector boost(const FourVector& frame, double frame_m2) const;
    FourVector boost(const FourVector& frame) const;
};

FourVector operator*(double factor, const FourVector& vector);

struct ChannelDefinition {
    std::vector<int> p;
    std::vector<int> f;
    std::vector<int> t;
    std::vector<int> d;
};

struct PropagatorDefinition {
    int mass_id = 0;
    int out = 0;
    int smin = 0;
    int smax = 0;
};

struct TimelikeInvariantDefinition {
    int out = 0;
    int smin = 0;
    int smax = 0;
};

struct SmaxDefinition {
    int process = 0;
    std::vector<int> out;
};

struct DecayDefinition {
    int out = 0;
    int out1 = 0;
    int out2 = 0;
};

struct TChannelDefinition {
    int mass_id = 0;
    int out = 0;
    int in1 = 0;
    int in2 = 0;
    int out1 = 0;
    int out2 = 0;
};

struct PhaseSpaceBlock {
    std::vector<std::vector<int>> smin;
    std::vector<SmaxDefinition> smax;
    std::vector<PropagatorDefinition> p;
    std::vector<TimelikeInvariantDefinition> f;
    std::vector<TChannelDefinition> t;
    std::vector<DecayDefinition> d;
    std::vector<ChannelDefinition> channels;
};

struct TaggedBlock {
    double r1 = 0.0;
    double r2 = 0.0;
    double gis1 = 1.0;
    double gis2 = 1.0;
    double s = 0.0;
    double s1 = 0.0;
    double s2 = 0.0;
    double smin = 0.0;
    double smax = 0.0;
    double sqrt_lambda = 0.0;
    double gcostheta = 0.5;
    bool has_s = false;
    bool has_s1 = false;
    bool has_s2 = false;
    bool has_smin = false;
    bool has_smax = false;
    bool has_sqrt_lambda = false;
    bool has_gcostheta = false;
};

struct BranchState {
    std::map<int, double> s;
    std::map<int, FourVector> p;
};

std::string trim(const std::string& value);
std::vector<std::string> split(const std::string& text, char delimiter);
std::vector<std::string> whitespace_tokens(const std::string& text);

PhaseSpaceBlock parse_phase_space_block(const std::string& path, const std::string& identifier);

double lambda(double x, double y, double z);
double h_propto_pot(double r, double smin, double smax, double exponent, double cut_technical);
double g_propto_pot(double s, double smin, double smax, double exponent, double cut_technical);
double h_propto_pot_mod(double r, double smin, double smax, double exponent, double cut_technical);
double g_propto_pot_mod(double s, double smin, double smax, double exponent, double cut_technical);
double c_phi(double r);

double smin_value(const PhaseSpaceBlock& phase_space, int index, double top_mass);
double smax_value(const PhaseSpaceBlock& phase_space, int index, double sqrt_shat, double top_mass);

double mapped_vanishing_width_s(double r, double smin, double smax, double mass2_mapping, double exponent);
double mapped_vanishing_width_g(double s, double smin, double smax, double mass2_mapping, double exponent);
double mapped_tchannel_s(double r, double tmin, double tmax, double mass2_mapping, double exponent);
double mapped_tchannel_g(double s_in1_out1, double tmin, double tmax, double mass2_mapping, double exponent);
double timelike_invariant_density(const PhaseSpaceBlock& phase_space, int index, double sqrt_shat, double top_mass);
double decay_density(double parent_s, double child1_s, double child2_s, double gis1, double gis2);
double tchannel_mass2_mapping(int mass_id, double top_mass, double mass0);

double apply_propagator(const PhaseSpaceBlock& phase_space,
                        const std::map<int, TaggedBlock>& p_blocks,
                        int index,
                        double sqrt_shat,
                        double top_mass,
                        double top_mass2,
                        double mass0,
                        double nuxs,
                        std::map<int, double>& local_s);

double apply_timelike_invariant(const PhaseSpaceBlock& phase_space,
                                const std::map<int, TaggedBlock>& f_blocks,
                                int index,
                                double sqrt_shat,
                                double top_mass,
                                std::map<int, double>& local_s);

double apply_tchannel(const PhaseSpaceBlock& phase_space,
                      const std::map<int, TaggedBlock>& t_blocks,
                      int index,
                      double top_mass,
                      double mass0,
                      double nuxt,
                      std::map<int, double>& local_s,
                      std::map<int, FourVector>& local_p);

double apply_decay(const PhaseSpaceBlock& phase_space,
                   const std::map<int, TaggedBlock>& d_blocks,
                   int index,
                   bool switch_decay_system,
                   int unrotated_all_final_bit,
                   std::map<int, double>& local_s,
                   std::map<int, FourVector>& local_p);

void fill_born_invariants(BranchState& state);

}  // namespace cs_phasespace