#include "cs_runner.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <system_error>

#include <LHAPDF/LHAPDF.h>

#include "born.hpp"
#include "common.hpp"
#include "cs_ca_integrand.hpp"
#include "cs_phasespace.hpp"
#include "cs_plot.hpp"
#include "cs_va_integrand.hpp"

namespace cs_ppttb {
namespace {

constexpr double kHcf = 1.98812436727894843e+07 * 8.0 * ttbar_qtsub::kPi * ttbar_qtsub::kPi * ttbar_qtsub::kPi;
constexpr double kTauExponent = 1.0;
constexpr double kMass0 = -1.0e-4;
constexpr double kNuxs = 1.1;
constexpr double kNuxt = 0.9;
constexpr double kFbPerGeV2 = ttbar_qtsub::kPbPerGeV2 * 1.0e3;

struct PDFCombination {
    int flav1 = 0;
    int flav2 = 0;
};

struct BornSubprocessSpec {
    std::string name;
    bool gluon_channel = false;
    std::vector<PDFCombination> pdf;
};

struct CASubprocessSpec {
    std::string name;
    ttbar_qtsub::CABornSubprocess subprocess = ttbar_qtsub::CABornSubprocess::UUbar_TTbar;
    std::string phase_space_id;
};

struct VASubprocessSpec {
    std::string name;
    ttbar_qtsub::CSChannel channel = ttbar_qtsub::CSChannel::QQbar;
    std::string phase_space_id;
    std::vector<PDFCombination> pdf;
};

struct Summary {
    long generated = 0;
    long accepted = 0;
    double sum_weights = 0.0;
    double sum_weight_squares = 0.0;

    void add(double weight, bool accepted_event) {
        ++generated;
        if (accepted_event) {
            ++accepted;
        }
        sum_weights += weight;
        sum_weight_squares += weight * weight;
    }

    double mean() const {
        return generated > 0 ? sum_weights / static_cast<double>(generated) : 0.0;
    }

    double variance() const {
        if (generated <= 1) {
            return 0.0;
        }
        const double n = static_cast<double>(generated);
        return std::max(0.0, (sum_weight_squares - sum_weights * sum_weights / n) / (n - 1.0));
    }

    double standard_error() const {
        return generated > 0 ? std::sqrt(variance() / static_cast<double>(generated)) : 0.0;
    }
};

struct PhaseSettings {
    cs_phasespace::PhaseSpaceBlock phase_space;
    std::vector<double> alpha;
};

struct GeneratedCAEvent {
    cs_phasespace::BranchState state;
    int selected_channel = -1;
    double tau = 0.0;
    double x1 = 0.0;
    double x2 = 0.0;
    double z1 = 1.0;
    double z2 = 1.0;
    double g_tau = 0.0;
    double g_x1x2 = 0.0;
    double g_mc = 0.0;
    double g_tot = 0.0;
    double g_z1 = 1.0;
    double g_z2 = 1.0;
    double ps_factor = 0.0;
};

struct GeneratedVAEvent {
    cs_phasespace::BranchState state;
    int selected_channel = -1;
    double tau = 0.0;
    double x1 = 0.0;
    double x2 = 0.0;
    double g_tau = 0.0;
    double g_x1x2 = 0.0;
    double g_mc = 0.0;
    double g_tot = 0.0;
    double ps_factor = 0.0;
};

struct ComponentComputation {
    double value = 0.0;
    double error = 0.0;
    std::vector<HistogramSeries> histograms;
};

struct ObservableDefinition {
    std::string key;
    std::string stem;
    std::string x_label;
    std::string y_label;
    std::string lo_title;
    std::string nlo_title;
};

void configure_lhapdf_environment() {
    LHAPDF::setVerbosity(0);
    if (std::getenv("LHAPDF_DATA_PATH") != nullptr) {
        return;
    }
    if (const char* conda_prefix = std::getenv("CONDA_PREFIX")) {
        const std::filesystem::path candidate = std::filesystem::path(conda_prefix) / "share" / "LHAPDF";
        if (std::filesystem::exists(candidate)) {
            setenv("LHAPDF_DATA_PATH", candidate.c_str(), 0);
            return;
        }
    }
    const std::filesystem::path fallback = "/opt/homebrew/Caskroom/miniconda/base/envs/lhapdf/share/LHAPDF";
    if (std::filesystem::exists(fallback)) {
        setenv("LHAPDF_DATA_PATH", fallback.c_str(), 0);
    }
}

std::unique_ptr<LHAPDF::PDF> make_pdf(const RunConfig& config) {
    configure_lhapdf_environment();
    return std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(config.pdf_set, config.pdf_member));
}

double uniform(std::mt19937_64& rng) {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double value = distribution(rng);
    if (value <= 0.0) {
        value = std::numeric_limits<double>::min();
    }
    if (value >= 1.0) {
        value = std::nextafter(1.0, 0.0);
    }
    return value;
}

std::string lowercase_copy(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return text;
}

const ObservableDefinition& observable_definition(const std::string& requested_observable) {
    static const ObservableDefinition mtt{"mtt", "mtt", "m_tt [GeV]", "d sigma / d m_tt [fb / GeV]", "LO m_tt spectrum", "NLO m_tt spectrum"};
    static const ObservableDefinition top_pt{"top_pt", "top_pt", "pT(t) [GeV]", "d sigma / d pT(t) [fb / GeV]", "LO top pT spectrum", "NLO top pT spectrum"};
    static const ObservableDefinition antitop_pt{"tb_pt", "tb_pt", "pT(tbar) [GeV]", "d sigma / d pT(tbar) [fb / GeV]", "LO antitop pT spectrum", "NLO antitop pT spectrum"};
    static const ObservableDefinition ytt{"y_tt", "ytt", "y_tt", "d sigma / d y_tt [fb]", "LO y_tt spectrum", "NLO y_tt spectrum"};

    const std::string key = lowercase_copy(requested_observable);
    if (key == "mtt") {
        return mtt;
    }
    if (key == "top_pt" || key == "t_pt" || key == "tpt") {
        return top_pt;
    }
    if (key == "tb_pt" || key == "antitop_pt" || key == "tbar_pt" || key == "tbpt") {
        return antitop_pt;
    }
    if (key == "ytt" || key == "y_tt") {
        return ytt;
    }
    throw std::runtime_error("Unsupported observable: " + requested_observable);
}

const ObservableDefinition& observable_definition(const RunConfig& config) {
    return observable_definition(config.observable);
}

double transverse_momentum(const cs_phasespace::FourVector& momentum) {
    return std::sqrt(momentum.x * momentum.x + momentum.y * momentum.y);
}

double rapidity_from_energy_pz(double energy, double pz) {
    const double plus = energy + pz;
    const double minus = energy - pz;
    if (!(plus > 0.0) || !(minus > 0.0)) {
        return 0.0;
    }
    return 0.5 * std::log(plus / minus);
}

double born_top_momentum_magnitude(double shat, double mt_gev) {
    return std::sqrt(std::max(0.0, 0.25 * shat - mt_gev * mt_gev));
}

std::vector<BornSubprocessSpec> born_subprocesses() {
    return {
        {"gg_tt~", true, {{0, 0}}},
        {"dd~_tt~", false, {{1, -1}, {-1, 1}, {3, -3}, {-3, 3}}},
        {"uu~_tt~", false, {{2, -2}, {-2, 2}, {4, -4}, {-4, 4}}},
        {"bb~_tt~", false, {{5, -5}, {-5, 5}}},
    };
}

std::vector<CASubprocessSpec> ca_subprocesses() {
    return {
        {"gg_tt~", ttbar_qtsub::CABornSubprocess::GG_TTbar, "gg_ttx 2 0 0"},
        {"dd~_tt~", ttbar_qtsub::CABornSubprocess::DDbar_TTbar, "ddx_ttx 2 0 0"},
        {"uu~_tt~", ttbar_qtsub::CABornSubprocess::UUbar_TTbar, "uux_ttx 2 0 0"},
        {"bb~_tt~", ttbar_qtsub::CABornSubprocess::BBbar_TTbar, "bbx_ttx 2 0 0"},
    };
}

std::vector<VASubprocessSpec> va_subprocesses() {
    return {
        {"gg_tt~", ttbar_qtsub::CSChannel::GG, "gg_ttx 2 0 0", {{0, 0}}},
        {"dd~_tt~", ttbar_qtsub::CSChannel::QQbar, "ddx_ttx 2 0 0", {{1, -1}, {-1, 1}, {3, -3}, {-3, 3}}},
        {"uu~_tt~", ttbar_qtsub::CSChannel::QQbar, "uux_ttx 2 0 0", {{2, -2}, {-2, 2}, {4, -4}, {-4, 4}}},
        {"bb~_tt~", ttbar_qtsub::CSChannel::QQbar, "bbx_ttx 2 0 0", {{5, -5}, {-5, 5}}},
    };
}

double hadronic_s(const RunConfig& config) {
    return config.sqrt_s_gev * config.sqrt_s_gev;
}

double xfx(LHAPDF::PDF& pdf, int flavor, double x, double mu_f) {
    return pdf.xfxQ(flavor, x, mu_f);
}

double pdf_factor(const std::vector<PDFCombination>& combinations,
                  LHAPDF::PDF& pdf,
                  double x1,
                  double x2,
                  double mu_f) {
    const double tau = x1 * x2;
    double result = 0.0;
    for (const PDFCombination& combination : combinations) {
        result += xfx(pdf, combination.flav1, x1, mu_f) * xfx(pdf, combination.flav2, x2, mu_f);
    }
    return result / tau;
}

cs_phasespace::TaggedBlock random_block(std::mt19937_64& rng, bool two_randoms) {
    cs_phasespace::TaggedBlock block;
    block.r1 = uniform(rng);
    block.r2 = two_randoms ? uniform(rng) : 0.0;
    return block;
}

std::map<int, cs_phasespace::TaggedBlock> random_blocks_for_indices(const std::vector<int>& indices,
                                                                    bool two_randoms,
                                                                    std::mt19937_64& rng) {
    std::map<int, cs_phasespace::TaggedBlock> blocks;
    for (int index : indices) {
        blocks[index] = random_block(rng, two_randoms);
    }
    return blocks;
}

int choose_channel(const std::vector<double>& alpha, std::mt19937_64& rng) {
    const double draw = uniform(rng);
    double cumulative = 0.0;
    for (std::size_t index = 0; index < alpha.size(); ++index) {
        cumulative += alpha[index];
        if (draw < cumulative || index + 1 == alpha.size()) {
            return static_cast<int>(index);
        }
    }
    throw std::runtime_error("Failed to choose a phase-space channel");
}

cs_phasespace::BranchState make_born_state(const cs_phasespace::PhaseSpaceBlock& phase_space,
                                           int channel_index,
                                           double shat,
                                           const RunConfig& config,
                                           std::mt19937_64& rng) {
    const double sqrt_shat = std::sqrt(shat);
    cs_phasespace::BranchState state;
    state.s[0] = shat;
    state.s[1] = 0.0;
    state.s[2] = 0.0;
    state.s[4] = config.mt_gev * config.mt_gev;
    state.s[8] = config.mt_gev * config.mt_gev;
    state.s[12] = shat;
    state.p[1] = {sqrt_shat / 2.0, 0.0, 0.0, sqrt_shat / 2.0};
    state.p[2] = {sqrt_shat / 2.0, 0.0, 0.0, -sqrt_shat / 2.0};
    state.p[0] = state.p[1] + state.p[2];
    state.p[12] = state.p[0];

    const cs_phasespace::ChannelDefinition& channel = phase_space.channels.at(static_cast<std::size_t>(channel_index));
    const auto p_blocks = random_blocks_for_indices(channel.p, false, rng);
    const auto f_blocks = random_blocks_for_indices(channel.f, false, rng);
    const auto t_blocks = random_blocks_for_indices(channel.t, true, rng);
    const auto d_blocks = random_blocks_for_indices(channel.d, true, rng);

    for (int index : channel.p) {
        cs_phasespace::apply_propagator(phase_space, p_blocks, index, sqrt_shat, config.mt_gev, config.mt_gev * config.mt_gev, kMass0, kNuxs, state.s);
    }
    for (int index : channel.f) {
        cs_phasespace::apply_timelike_invariant(phase_space, f_blocks, index, sqrt_shat, config.mt_gev, state.s);
    }
    for (int index : channel.t) {
        cs_phasespace::apply_tchannel(phase_space, t_blocks, index, config.mt_gev, kMass0, kNuxt, state.s, state.p);
    }
    for (int index : channel.d) {
        cs_phasespace::apply_decay(phase_space, d_blocks, index, true, 12, state.s, state.p);
    }
    cs_phasespace::fill_born_invariants(state);
    return state;
}

double born_channel_density(const cs_phasespace::PhaseSpaceBlock& phase_space,
                            const cs_phasespace::BranchState& state,
                            const cs_phasespace::ChannelDefinition& channel,
                            const RunConfig& config) {
    const double sqrt_shat = std::sqrt(state.s.at(0));
    double product = 1.0;
    for (int index : channel.p) {
        const cs_phasespace::PropagatorDefinition& propagator = phase_space.p.at(index);
        const double p_smin = cs_phasespace::smin_value(phase_space, propagator.smin, config.mt_gev);
        const double p_smax = cs_phasespace::smax_value(phase_space, propagator.smax, sqrt_shat, config.mt_gev);
        const double mass2_mapping = ((propagator.mass_id == 6) ? config.mt_gev * config.mt_gev : 0.0) + kMass0;
        product *= cs_phasespace::mapped_vanishing_width_g(state.s.at(propagator.out), p_smin, p_smax, mass2_mapping, kNuxs);
    }
    for (int index : channel.f) {
        product *= cs_phasespace::timelike_invariant_density(phase_space, index, sqrt_shat, config.mt_gev);
    }
    for (int index : channel.t) {
        const cs_phasespace::TChannelDefinition& tchannel = phase_space.t.at(index);
        const double s_out = state.s.at(tchannel.out);
        const double s_out1 = state.s.at(tchannel.out1);
        const double s_out2 = state.s.at(tchannel.out2);
        const double s_in1 = state.s.at(tchannel.in1);
        const double s_in2 = state.s.at(tchannel.in2);
        const double term1 = (s_out + s_out1 - s_out2) * (s_out + s_in1 - s_in2);
        const double term2in = std::sqrt(std::max(0.0, cs_phasespace::lambda(s_out, s_in1, s_in2)));
        const double term2out = std::sqrt(std::max(0.0, cs_phasespace::lambda(s_out, s_out1, s_out2)));
        const double term2 = term2in * term2out;
        if (!(term2 > 0.0)) {
            return 0.0;
        }
        double tmin = s_out1 + s_in1 - (term1 + term2) / (2.0 * s_out);
        double tmax = s_out1 + s_in1 - (term1 - term2) / (2.0 * s_out);
        if (tmax > 0.0) {
            tmax = 0.0;
        }
        const double mass2_mapping = cs_phasespace::tchannel_mass2_mapping(tchannel.mass_id, config.mt_gev, kMass0);
        const double t_value = (state.p.at(tchannel.in1) - state.p.at(tchannel.out1)).m2();
        product *= cs_phasespace::mapped_tchannel_g(t_value, tmin, tmax, mass2_mapping, kNuxt) * (2.0 * term2in) / ttbar_qtsub::kPi;
    }
    for (int index : channel.d) {
        const cs_phasespace::DecayDefinition& decay = phase_space.d.at(index);
        product *= cs_phasespace::decay_density(state.s.at(decay.out), state.s.at(decay.out1), state.s.at(decay.out2), 1.0, 1.0);
    }
    return std::isfinite(product) && product > 0.0 ? product : 0.0;
}

PhaseSettings make_settings(const RunConfig& config, const std::string& phase_space_id) {
    PhaseSettings settings;
    settings.phase_space = cs_phasespace::parse_phase_space_block(config.born_phase_space.string(), phase_space_id);
    settings.alpha.assign(settings.phase_space.channels.size(), 1.0 / static_cast<double>(settings.phase_space.channels.size()));
    return settings;
}

std::vector<ttbar_qtsub::FourVector> born_momenta(const cs_phasespace::BranchState& state) {
    return {
        {state.p.at(1).e, state.p.at(1).x, state.p.at(1).y, state.p.at(1).z},
        {state.p.at(2).e, state.p.at(2).x, state.p.at(2).y, state.p.at(2).z},
        {state.p.at(4).e, state.p.at(4).x, state.p.at(4).y, state.p.at(4).z},
        {state.p.at(8).e, state.p.at(8).x, state.p.at(8).y, state.p.at(8).z},
    };
}

GeneratedCAEvent generate_ca_event(const RunConfig& config,
                                   const PhaseSettings& settings,
                                   std::mt19937_64& rng) {
    GeneratedCAEvent event;
    const double tau0 = 4.0 * config.mt_gev * config.mt_gev / hadronic_s(config);
    event.tau = cs_phasespace::h_propto_pot(uniform(rng), tau0, 1.0, kTauExponent, 0.0);
    event.g_tau = cs_phasespace::g_propto_pot(event.tau, tau0, 1.0, kTauExponent, 0.0);
    event.x1 = cs_phasespace::h_propto_pot(uniform(rng), event.tau, 1.0, 1.0, 0.0);
    event.x2 = event.tau / event.x1;
    const double min_x1 = std::max(event.tau, 1.0e-7);
    event.g_x1x2 = 1.0 / (-std::log(min_x1));
    event.selected_channel = choose_channel(settings.alpha, rng);
    event.state = make_born_state(settings.phase_space, event.selected_channel, event.tau * hadronic_s(config), config, rng);
    for (std::size_t channel_index = 0; channel_index < settings.phase_space.channels.size(); ++channel_index) {
        event.g_mc += settings.alpha[channel_index] * born_channel_density(settings.phase_space, event.state, settings.phase_space.channels[channel_index], config);
    }
    event.z1 = std::exp(uniform(rng) * std::log(event.x1));
    event.z2 = std::exp(uniform(rng) * std::log(event.x2));
    event.g_z1 = 1.0 / (-std::log(event.x1));
    event.g_z2 = 1.0 / (-std::log(event.x2));
    event.g_tot = event.g_mc * event.g_tau * event.g_x1x2;
    event.ps_factor = kHcf / (event.g_tot * event.state.s.at(0));
    return event;
}

GeneratedVAEvent generate_va_event(const RunConfig& config,
                                   const PhaseSettings& settings,
                                   std::mt19937_64& rng) {
    GeneratedVAEvent event;
    const double tau0 = 4.0 * config.mt_gev * config.mt_gev / hadronic_s(config);
    event.tau = cs_phasespace::h_propto_pot(uniform(rng), tau0, 1.0, kTauExponent, 0.0);
    event.g_tau = cs_phasespace::g_propto_pot(event.tau, tau0, 1.0, kTauExponent, 0.0);
    event.x1 = cs_phasespace::h_propto_pot(uniform(rng), event.tau, 1.0, 1.0, 0.0);
    event.x2 = event.tau / event.x1;
    const double min_x1 = std::max(event.tau, 1.0e-7);
    event.g_x1x2 = 1.0 / (-std::log(min_x1));
    event.selected_channel = choose_channel(settings.alpha, rng);
    event.state = make_born_state(settings.phase_space, event.selected_channel, event.tau * hadronic_s(config), config, rng);
    for (std::size_t channel_index = 0; channel_index < settings.phase_space.channels.size(); ++channel_index) {
        event.g_mc += settings.alpha[channel_index] * born_channel_density(settings.phase_space, event.state, settings.phase_space.channels[channel_index], config);
    }
    event.g_tot = event.g_mc * event.g_tau * event.g_x1x2;
    event.ps_factor = kHcf / (event.g_tot * event.state.s.at(0));
    return event;
}

double mtt_from_state(const cs_phasespace::BranchState& state) {
    return std::sqrt(std::max(0.0, state.s.at(12)));
}

double hadronic_pair_rapidity(double x1, double x2) {
    return rapidity_from_energy_pz(x1 + x2, x1 - x2);
}

double observable_from_lo_event(const ObservableDefinition& observable,
                                const RunConfig& config,
                                double shat,
                                double x1,
                                double x2,
                                double costheta) {
    if (observable.key == "mtt") {
        return std::sqrt(shat);
    }
    if (observable.key == "top_pt" || observable.key == "tb_pt") {
        const double sin_theta = std::sqrt(std::max(0.0, 1.0 - costheta * costheta));
        return born_top_momentum_magnitude(shat, config.mt_gev) * sin_theta;
    }
    if (observable.key == "y_tt") {
        return hadronic_pair_rapidity(x1, x2);
    }
    throw std::runtime_error("Unsupported observable in LO generation: " + observable.key);
}

double observable_from_state(const ObservableDefinition& observable,
                             const RunConfig& config,
                             const cs_phasespace::BranchState& state,
                             double x1 = std::numeric_limits<double>::quiet_NaN(),
                             double x2 = std::numeric_limits<double>::quiet_NaN()) {
    if (observable.key == "mtt") {
        return mtt_from_state(state);
    }
    if (observable.key == "top_pt") {
        return transverse_momentum(state.p.at(4));
    }
    if (observable.key == "tb_pt") {
        return transverse_momentum(state.p.at(8));
    }
    if (observable.key == "y_tt") {
        const auto& top = state.p.at(4);
        const auto& antitop = state.p.at(8);
        double rapidity = rapidity_from_energy_pz(top.e + antitop.e, top.z + antitop.z);
        if (std::isfinite(x1) && std::isfinite(x2) && x1 > 0.0 && x2 > 0.0) {
            rapidity += hadronic_pair_rapidity(x1, x2);
        }
        return rapidity;
    }
    throw std::runtime_error("Unsupported observable from state: " + observable.key);
}

std::map<std::string, std::string> read_tabbed_entries(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Failed to open artifact " + path.string());
    }
    std::map<std::string, std::string> entries;
    std::string line;
    while (std::getline(input, line)) {
        if (line == "Per-subprocess") {
            break;
        }
        std::stringstream stream(line);
        std::string candidate;
        if (std::getline(stream, candidate, '\t')) {
            std::string value;
            if (std::getline(stream, value)) {
                entries.emplace(candidate, value);
            }
        }
    }
    return entries;
}

double read_tabbed_value(const std::map<std::string, std::string>& entries,
                         const std::filesystem::path& path,
                         const std::string& key) {
    const auto found = entries.find(key);
    if (found != entries.end()) {
        return std::stod(found->second);
    }
    throw std::runtime_error("Key '" + key + "' not found in artifact " + path.string());
}

std::optional<double> optional_tabbed_value(const std::map<std::string, std::string>& entries,
                                           const std::string& key) {
    const auto found = entries.find(key);
    if (found == entries.end()) {
        return std::nullopt;
    }
    return std::stod(found->second);
}

bool nearly_equal(double left, double right) {
    const double scale = std::max({1.0, std::abs(left), std::abs(right)});
    return std::abs(left - right) <= 1.0e-8 * scale;
}

void require_matching_artifact_value(const std::map<std::string, std::string>& entries,
                                     const std::filesystem::path& path,
                                     const std::string& key,
                                     double expected) {
    const std::optional<double> found = optional_tabbed_value(entries, key);
    if (!found.has_value()) {
        throw std::runtime_error("RA artifact " + path.string() + " is missing required metadata key '" + key + "'");
    }
    if (!nearly_equal(*found, expected)) {
        std::ostringstream message;
        message << "RA artifact " << path.string() << " has " << key << "=" << *found
                << " but the run config requests " << expected;
        throw std::runtime_error(message.str());
    }
}

void validate_ra_artifact_metadata(const RunConfig& config,
                                   const std::map<std::string, std::string>& entries,
                                   const std::filesystem::path& path) {
    const bool has_modern_metadata = entries.count("sqrt_s_gev") > 0
        || entries.count("mt_gev") > 0
        || entries.count("muF_gev") > 0
        || entries.count("muR_gev") > 0;

    if (!has_modern_metadata) {
        if (!nearly_equal(config.sqrt_s_gev, 8000.0)) {
            std::ostringstream message;
            message << "RA artifact " << path.string()
                    << " has no sqrt_s_gev metadata, so it is treated as a legacy 8 TeV artifact. "
                    << "Regenerate the RA artifact with --sqrt-s-gev " << config.sqrt_s_gev
                    << " or point ra_artifact to a matching metadata-rich file.";
            throw std::runtime_error(message.str());
        }
        return;
    }

    require_matching_artifact_value(entries, path, "sqrt_s_gev", config.sqrt_s_gev);
    require_matching_artifact_value(entries, path, "mt_gev", config.mt_gev);
    require_matching_artifact_value(entries, path, "muF_gev", config.muF_gev);
    require_matching_artifact_value(entries, path, "muR_gev", config.muR_gev);

    const auto pdf_set = entries.find("pdf_set");
    if (pdf_set != entries.end() && pdf_set->second != config.pdf_set) {
        throw std::runtime_error("RA artifact " + path.string() + " uses pdf_set=" + pdf_set->second
            + " but the run config requests " + config.pdf_set);
    }
    const std::optional<double> pdf_member = optional_tabbed_value(entries, "pdf_member");
    if (pdf_member.has_value() && static_cast<int>(*pdf_member) != config.pdf_member) {
        std::ostringstream message;
        message << "RA artifact " << path.string() << " uses pdf_member=" << static_cast<int>(*pdf_member)
                << " but the run config requests " << config.pdf_member;
        throw std::runtime_error(message.str());
    }
}

HistogramSeries load_histogram_csv(const std::filesystem::path& path,
                                   const std::string& name,
                                   const std::string& x_label,
                                   const std::string& y_label) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Failed to open histogram CSV " + path.string());
    }
    HistogramSeries series;
    series.name = name;
    series.x_label = x_label;
    series.y_label = y_label;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        if (line.find("bin_low_gev") != std::string::npos) {
            continue;
        }
        std::stringstream stream(line);
        std::string low_text;
        std::string high_text;
        std::string value_text;
        std::string error_text;
        if (!std::getline(stream, low_text, ',')) {
            continue;
        }
        if (!std::getline(stream, high_text, ',')) {
            throw std::runtime_error("Malformed histogram CSV row in " + path.string());
        }
        if (!std::getline(stream, value_text, ',')) {
            throw std::runtime_error("Malformed histogram CSV row in " + path.string());
        }
        if (!std::getline(stream, error_text, ',')) {
            throw std::runtime_error("Malformed histogram CSV row in " + path.string());
        }
        series.bins.push_back({std::stod(low_text), std::stod(high_text), std::stod(value_text), std::stod(error_text)});
    }
    return series;
}

void write_histogram_csv(const HistogramSeries& series, const std::filesystem::path& path) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("Failed to write histogram CSV " + path.string());
    }
    output << std::scientific << std::setprecision(17);
    output << "bin_low_gev,bin_high_gev,value,error\n";
    for (const HistogramBin& bin : series.bins) {
        output << bin.low << ',' << bin.high << ',' << bin.value << ',' << bin.error << '\n';
    }
}

void write_text_file(const std::filesystem::path& path, const std::string& text) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("Failed to write file " + path.string());
    }
    output << text;
}

void remove_if_exists(const std::filesystem::path& path) {
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
}

void clear_distribution_outputs(const ObservableDefinition& observable,
                                const std::filesystem::path& csv_dir,
                                const std::filesystem::path& png_dir,
                                const std::filesystem::path& pdf_dir) {
    for (const char* suffix : {"_lo.csv", "_ra.csv", "_ca.csv", "_va.csv", "_nlo.csv", "_nlo_unavailable.txt"}) {
        remove_if_exists(csv_dir / (observable.stem + suffix));
    }
    for (const char* suffix : {"_lo.png", "_nlo.png"}) {
        remove_if_exists(png_dir / (observable.stem + suffix));
    }
    for (const char* suffix : {"_lo.pdf", "_nlo.pdf"}) {
        remove_if_exists(pdf_dir / (observable.stem + suffix));
    }
}

std::vector<HistogramAccumulator> make_histogram_accumulators(const RunConfig& config) {
    std::vector<HistogramAccumulator> accumulators;
    accumulators.reserve(config.distributions.size());
    for (const DistributionConfig& distribution : config.distributions) {
        accumulators.emplace_back(distribution.bins, distribution.min_gev, distribution.max_gev);
    }
    return accumulators;
}

std::vector<HistogramSeries> finalize_histograms(const std::vector<HistogramAccumulator>& accumulators,
                                                 const std::vector<DistributionConfig>& distributions,
                                                 const std::string& suffix) {
    std::vector<HistogramSeries> histograms;
    histograms.reserve(accumulators.size());
    for (std::size_t index = 0; index < accumulators.size(); ++index) {
        const ObservableDefinition& observable = observable_definition(distributions[index].observable);
        histograms.push_back(accumulators[index].finalize(
            observable.stem + suffix,
            observable.x_label,
            observable.y_label));
    }
    return histograms;
}

ComponentComputation compute_lo(const RunConfig& config, LHAPDF::PDF& pdf) {
    const double tau0 = 4.0 * config.mt_gev * config.mt_gev / hadronic_s(config);
    const double tau_log_range = std::log(1.0 / tau0);
    const double alpha_s = pdf.alphasQ(config.muR_gev);
    const double gs = std::sqrt(4.0 * ttbar_qtsub::kPi * alpha_s);
    const auto subprocesses = born_subprocesses();

    std::mt19937_64 rng(config.seed);
    Summary summary;
    std::vector<HistogramAccumulator> histograms = make_histogram_accumulators(config);

    for (long event = 0; event < config.lo_samples; ++event) {
        const double tau = tau0 * std::exp(uniform(rng) * tau_log_range);
        const double g_tau = 1.0 / (tau * tau_log_range);
        const double x1_log_range = std::log(1.0 / tau);
        const double x1 = tau * std::exp(uniform(rng) * x1_log_range);
        const double x2 = tau / x1;
        const double g_x1 = 1.0 / (x1 * x1_log_range);
        const double costheta = -1.0 + 2.0 * uniform(rng);
        const double g_costheta = 0.5;
        const double shat = tau * hadronic_s(config);
        const double beta = std::sqrt(std::max(0.0, 1.0 - 4.0 * config.mt_gev * config.mt_gev / shat));
        const double phase_flux = beta / (32.0 * ttbar_qtsub::kPi * shat);

        double dsigma_dcostheta = 0.0;
        for (const BornSubprocessSpec& subprocess : subprocesses) {
            const double born = subprocess.gluon_channel
                ? ttbar_qtsub::born_gg(shat, config.mt_gev, config.muR_gev, gs, costheta)
                : ttbar_qtsub::born_qq(shat, config.mt_gev, config.muR_gev, gs, costheta);
            dsigma_dcostheta += pdf_factor(subprocess.pdf, pdf, x1, x2, config.muF_gev) * phase_flux * born * kFbPerGeV2 / x1;
        }

        const double weight = dsigma_dcostheta / (g_tau * g_x1 * g_costheta);
        summary.add(weight, true);
        for (std::size_t distribution_index = 0; distribution_index < config.distributions.size(); ++distribution_index) {
            const ObservableDefinition& observable = observable_definition(config.distributions[distribution_index].observable);
            histograms[distribution_index].observe(
                observable_from_lo_event(observable, config, shat, x1, x2, costheta),
                weight);
        }
    }

    return {summary.mean(), summary.standard_error(), finalize_histograms(histograms, config.distributions, "_lo")};
}

ComponentComputation compute_ca(const RunConfig& config, LHAPDF::PDF& pdf) {
    (void)pdf;
    const double alpha_s = pdf.alphasQ(config.muR_gev);
    std::vector<std::vector<HistogramSeries>> subprocess_histograms(config.distributions.size());
    double combined_value = 0.0;
    double combined_error2 = 0.0;
    const auto specs = ca_subprocesses();

    for (std::size_t spec_index = 0; spec_index < specs.size(); ++spec_index) {
        const CASubprocessSpec& spec = specs[spec_index];
        const PhaseSettings settings = make_settings(config, spec.phase_space_id);
        std::mt19937_64 rng(config.seed + 1000ULL * static_cast<unsigned long long>(spec_index + 1));
        Summary summary;
        std::vector<HistogramAccumulator> histograms = make_histogram_accumulators(config);

        for (long generated = 0; generated < config.ca_samples; ++generated) {
            GeneratedCAEvent event = generate_ca_event(config, settings, rng);
            bool accepted = event.g_tot > 0.0 && std::isfinite(event.ps_factor);
            double weight = 0.0;
            if (accepted) {
                try {
                    const ttbar_qtsub::CAInput input{
                        spec.subprocess,
                        born_momenta(event.state),
                        event.x1,
                        event.x2,
                        event.z1,
                        event.z2,
                        alpha_s,
                        config.muF_gev,
                        config.muR_gev,
                        config.mt_gev,
                        event.g_z1,
                        event.g_z2,
                    };
                    const auto result = ttbar_qtsub::compute_ca_integrand(input);
                    weight = event.ps_factor * result.sum_integrand_cv;
                } catch (const std::exception&) {
                    accepted = false;
                    weight = 0.0;
                }
            }
            if (!std::isfinite(weight)) {
                accepted = false;
                weight = 0.0;
            }
            summary.add(weight, accepted);
            for (std::size_t distribution_index = 0; distribution_index < config.distributions.size(); ++distribution_index) {
                const ObservableDefinition& observable = observable_definition(config.distributions[distribution_index].observable);
                histograms[distribution_index].observe(
                    observable_from_state(observable, config, event.state, event.x1, event.x2),
                    weight);
            }
        }

        combined_value += summary.mean();
        combined_error2 += summary.standard_error() * summary.standard_error();
        std::vector<HistogramSeries> finalized = finalize_histograms(histograms, config.distributions, "_ca_" + spec.name);
        for (std::size_t distribution_index = 0; distribution_index < finalized.size(); ++distribution_index) {
            subprocess_histograms[distribution_index].push_back(finalized[distribution_index]);
        }
    }

    std::vector<HistogramSeries> combined_histograms;
    combined_histograms.reserve(config.distributions.size());
    for (std::size_t distribution_index = 0; distribution_index < config.distributions.size(); ++distribution_index) {
        const ObservableDefinition& observable = observable_definition(config.distributions[distribution_index].observable);
        combined_histograms.push_back(combine_histogram_series(
            observable.stem + "_ca",
            observable.x_label,
            observable.y_label,
            subprocess_histograms[distribution_index]));
    }
    return {combined_value, std::sqrt(combined_error2), combined_histograms};
}

ComponentComputation compute_va(const RunConfig& config, LHAPDF::PDF& pdf) {
    const double alpha_s = pdf.alphasQ(config.muR_gev);
    std::vector<std::vector<HistogramSeries>> subprocess_histograms(config.distributions.size());
    double combined_value = 0.0;
    double combined_error2 = 0.0;
    const auto specs = va_subprocesses();

    for (std::size_t spec_index = 0; spec_index < specs.size(); ++spec_index) {
        const VASubprocessSpec& spec = specs[spec_index];
        const PhaseSettings settings = make_settings(config, spec.phase_space_id);
        std::mt19937_64 rng(config.seed + 2000ULL * static_cast<unsigned long long>(spec_index + 1));
        Summary summary;
        std::vector<HistogramAccumulator> histograms = make_histogram_accumulators(config);

        for (long generated = 0; generated < config.va_samples; ++generated) {
            GeneratedVAEvent event = generate_va_event(config, settings, rng);
            bool accepted = event.g_tot > 0.0 && std::isfinite(event.ps_factor);
            double weight = 0.0;
            if (accepted) {
                try {
                    const ttbar_qtsub::VAIntegrandFactors factors{
                        event.ps_factor,
                        pdf_factor(spec.pdf, pdf, event.x1, event.x2, config.muF_gev),
                        1.0,
                        1.0,
                    };
                    const auto result = ttbar_qtsub::compute_va_integrand(
                        spec.channel,
                        born_momenta(event.state),
                        alpha_s,
                        config.muR_gev,
                        config.mt_gev,
                        factors);
                    weight = result.integrand;
                } catch (const std::exception&) {
                    accepted = false;
                    weight = 0.0;
                }
            }
            if (!std::isfinite(weight)) {
                accepted = false;
                weight = 0.0;
            }
            summary.add(weight, accepted);
            for (std::size_t distribution_index = 0; distribution_index < config.distributions.size(); ++distribution_index) {
                const ObservableDefinition& observable = observable_definition(config.distributions[distribution_index].observable);
                histograms[distribution_index].observe(
                    observable_from_state(observable, config, event.state, event.x1, event.x2),
                    weight);
            }
        }

        combined_value += summary.mean();
        combined_error2 += summary.standard_error() * summary.standard_error();
        std::vector<HistogramSeries> finalized = finalize_histograms(histograms, config.distributions, "_va_" + spec.name);
        for (std::size_t distribution_index = 0; distribution_index < finalized.size(); ++distribution_index) {
            subprocess_histograms[distribution_index].push_back(finalized[distribution_index]);
        }
    }

    std::vector<HistogramSeries> combined_histograms;
    combined_histograms.reserve(config.distributions.size());
    for (std::size_t distribution_index = 0; distribution_index < config.distributions.size(); ++distribution_index) {
        const ObservableDefinition& observable = observable_definition(config.distributions[distribution_index].observable);
        combined_histograms.push_back(combine_histogram_series(
            observable.stem + "_va",
            observable.x_label,
            observable.y_label,
            subprocess_histograms[distribution_index]));
    }
    return {combined_value, std::sqrt(combined_error2), combined_histograms};
}

ComponentResult load_ra_component(const RunConfig& config) {
    if (config.ra_mode != "artifact") {
        throw std::runtime_error("Only ra_mode=artifact is supported in the clean package at present");
    }
    const auto artifact_entries = read_tabbed_entries(config.ra_artifact);
    validate_ra_artifact_metadata(config, artifact_entries, config.ra_artifact);
    ComponentResult result;
    result.name = "RA";
    result.value = read_tabbed_value(artifact_entries, config.ra_artifact, "integral_estimate");
    result.error = read_tabbed_value(artifact_entries, config.ra_artifact, "standard_error");
    result.source = config.ra_artifact.string();
    result.artifact_only = true;
    return result;
}

double quadrature_error(const std::vector<ComponentResult>& components) {
    double error2 = 0.0;
    for (const ComponentResult& component : components) {
        error2 += component.error * component.error;
    }
    return std::sqrt(error2);
}

std::string format_value(double value) {
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(17) << value;
    return stream.str();
}

}  // namespace

RunReport run_package(const RunConfig& config) {
    auto pdf = make_pdf(config);

    const auto lo = compute_lo(config, *pdf);
    const auto ca = compute_ca(config, *pdf);
    const auto va = compute_va(config, *pdf);
    auto ra = load_ra_component(config);

    std::vector<DistributionResult> distributions;
    distributions.reserve(config.distributions.size());
    bool any_ra_histogram = false;
    for (std::size_t distribution_index = 0; distribution_index < config.distributions.size(); ++distribution_index) {
        const DistributionConfig& distribution_config = config.distributions[distribution_index];
        const ObservableDefinition& observable = observable_definition(distribution_config.observable);
        DistributionResult distribution;
        distribution.config = distribution_config;
        distribution.lo_histogram = lo.histograms[distribution_index];
        distribution.ca_histogram = ca.histograms[distribution_index];
        distribution.va_histogram = va.histograms[distribution_index];

        if (distribution_config.write_nlo
            && !distribution_config.ra_histogram_csv.empty()
            && std::filesystem::exists(distribution_config.ra_histogram_csv)) {
            distribution.ra_histogram = load_histogram_csv(
                distribution_config.ra_histogram_csv,
                observable.stem + "_ra",
                observable.x_label,
                observable.y_label);
            distribution.has_ra_histogram = true;
            any_ra_histogram = true;
        }

        if (distribution.has_ra_histogram) {
            distribution.has_nlo_histogram = true;
            distribution.nlo_histogram = combine_histogram_series(
                observable.stem + "_nlo",
                observable.x_label,
                observable.y_label,
                {distribution.lo_histogram, distribution.ra_histogram, distribution.ca_histogram, distribution.va_histogram});
        }
        distributions.push_back(distribution);
    }
    ra.has_histogram = any_ra_histogram;

    RunReport report;
    report.config = config;
    report.components = {
        {"LO", lo.value, lo.error, "local Monte Carlo (analytic Born)", false, !lo.histograms.empty(), lo.histograms.empty() ? HistogramSeries{} : lo.histograms.front()},
        ra,
        {"CA", ca.value, ca.error, "local Monte Carlo (analytic collinear counterterm)", false, !ca.histograms.empty(), ca.histograms.empty() ? HistogramSeries{} : ca.histograms.front()},
        {"VA", va.value, va.error, "local Monte Carlo (analytic virtual+integrated dipoles)", false, !va.histograms.empty(), va.histograms.empty() ? HistogramSeries{} : va.histograms.front()},
    };
    report.distributions = distributions;

    report.nlo_total = std::accumulate(report.components.begin(), report.components.end(), 0.0,
        [](double value, const ComponentResult& component) { return value + component.value; });
    report.nlo_error = quadrature_error(report.components);
    if (!distributions.empty()) {
        report.lo_histogram = distributions.front().lo_histogram;
        report.has_nlo_histogram = distributions.front().has_nlo_histogram;
        report.nlo_histogram = distributions.front().nlo_histogram;
    }

    for (const DistributionResult& distribution : distributions) {
        if (distribution.config.write_nlo && !distribution.has_nlo_histogram) {
            report.notes.push_back(
                distribution.config.observable
                + " NLO distribution was not assembled because no RA histogram was configured or found.");
        }
    }
    report.notes.push_back("RA is intentionally packaged as a validated artifact until a fully local spin-correlated real-emission kernel is implemented.");
    return report;
}

void write_run_outputs(const RunReport& report) {
    const std::filesystem::path output_dir = report.config.output_dir;
    const std::filesystem::path csv_dir = output_dir / "csv";
    const std::filesystem::path png_dir = output_dir / "plots" / "png";
    const std::filesystem::path pdf_dir = output_dir / "plots" / "pdf";

    std::filesystem::create_directories(csv_dir);
    std::filesystem::create_directories(png_dir);
    std::filesystem::create_directories(pdf_dir);

    for (const DistributionResult& distribution : report.distributions) {
        clear_distribution_outputs(observable_definition(distribution.config.observable), csv_dir, png_dir, pdf_dir);
    }

    {
        std::ofstream output(output_dir / "summary.txt");
        if (!output) {
            throw std::runtime_error("Failed to open summary output");
        }
        output << std::scientific << std::setprecision(17);
        output << "package	cs_ppttb\n";
        output << "config	" << report.config.config_path.string() << '\n';
        output << "sqrt_s_gev	" << report.config.sqrt_s_gev << '\n';
        output << "mt_gev	" << report.config.mt_gev << '\n';
        output << "muF_gev	" << report.config.muF_gev << '\n';
        output << "muR_gev	" << report.config.muR_gev << '\n';
        output << "pdf_set	" << report.config.pdf_set << '\n';
        output << "pdf_member	" << report.config.pdf_member << '\n';
        output << "observable\t" << report.config.observable << '\n';
        output << "lo_samples	" << report.config.lo_samples << '\n';
        output << "ca_samples	" << report.config.ca_samples << '\n';
        output << "va_samples	" << report.config.va_samples << '\n';
        output << "seed	" << report.config.seed << '\n';
        output << "nlo_total_fb	" << report.nlo_total << '\n';
        output << "nlo_total_pb	" << report.nlo_total * 1.0e-3 << '\n';
        output << "nlo_standard_error_fb	" << report.nlo_error << '\n';
        output << "nlo_standard_error_pb	" << report.nlo_error * 1.0e-3 << '\n';
        output << "ra_mode	" << report.config.ra_mode << '\n';
        output << "distributions	" << report.distributions.size() << '\n';
        for (const DistributionResult& distribution : report.distributions) {
            output << "distribution	" << distribution.config.observable
                   << "	bins=" << distribution.config.bins
                   << "	min=" << distribution.config.min_gev
                   << "	max=" << distribution.config.max_gev
                   << "	lo=" << (distribution.config.write_lo ? 1 : 0)
                   << "	nlo=" << (distribution.config.write_nlo ? 1 : 0)
                   << "	has_ra_histogram=" << (distribution.has_ra_histogram ? 1 : 0)
                   << '\n';
        }
        output << '\n';
        output << "Component breakdown\n";
        output << "component	value_fb	standard_error_fb	source	artifact_only	has_histogram\n";
        for (const ComponentResult& component : report.components) {
            output << component.name << '\t'
                   << component.value << '\t'
                   << component.error << '\t'
                   << component.source << '\t'
                   << (component.artifact_only ? 1 : 0) << '\t'
                   << (component.has_histogram ? 1 : 0) << '\n';
        }
        if (!report.notes.empty()) {
            output << '\n';
            output << "Notes\n";
            for (const std::string& note : report.notes) {
                output << "- " << note << '\n';
            }
        }
    }

    {
        std::ofstream output(csv_dir / "component_totals.tsv");
        if (!output) {
            throw std::runtime_error("Failed to write component totals TSV");
        }
        output << std::scientific << std::setprecision(17);
        output << "component\tvalue_fb\tstandard_error_fb\tsource\tartifact_only\thas_histogram\n";
        for (const ComponentResult& component : report.components) {
            output << component.name << '\t'
                   << component.value << '\t'
                   << component.error << '\t'
                   << component.source << '\t'
                   << (component.artifact_only ? 1 : 0) << '\t'
                   << (component.has_histogram ? 1 : 0) << '\n';
        }
        output << "NLO\t" << report.nlo_total << '\t' << report.nlo_error << "\tassembled\t0\t" << (report.has_nlo_histogram ? 1 : 0) << '\n';
    }

    if (report.config.write_csv) {
        for (const DistributionResult& distribution : report.distributions) {
            const ObservableDefinition& observable = observable_definition(distribution.config.observable);
            if (distribution.config.write_lo) {
                write_histogram_csv(distribution.lo_histogram, csv_dir / (observable.stem + "_lo.csv"));
            }
            if (distribution.config.write_nlo) {
                write_histogram_csv(distribution.ca_histogram, csv_dir / (observable.stem + "_ca.csv"));
                write_histogram_csv(distribution.va_histogram, csv_dir / (observable.stem + "_va.csv"));
                if (distribution.has_ra_histogram) {
                    write_histogram_csv(distribution.ra_histogram, csv_dir / (observable.stem + "_ra.csv"));
                }
                if (distribution.has_nlo_histogram) {
                    write_histogram_csv(distribution.nlo_histogram, csv_dir / (observable.stem + "_nlo.csv"));
                } else {
                    std::ostringstream note;
                    note << distribution.config.observable
                         << " NLO distribution was not assembled because no RA histogram was configured or found.\n";
                    write_text_file(csv_dir / (observable.stem + "_nlo_unavailable.txt"), note.str());
                }
            }
        }
    }

    if (report.config.write_plots) {
        for (const DistributionResult& distribution : report.distributions) {
            const ObservableDefinition& observable = observable_definition(distribution.config.observable);
            HistogramSeries lo_reference;
            const HistogramSeries* lo_overlay = nullptr;
            if (!distribution.config.reference_lo_csv.empty() && std::filesystem::exists(distribution.config.reference_lo_csv)) {
                lo_reference = load_histogram_csv(distribution.config.reference_lo_csv, "reference_lo", observable.x_label, observable.y_label);
                lo_overlay = &lo_reference;
            }

            if (distribution.config.write_lo) {
                write_histogram_plot_png(distribution.lo_histogram, png_dir / (observable.stem + "_lo.png"), observable.lo_title, lo_overlay);
                write_histogram_plot_pdf(distribution.lo_histogram, pdf_dir / (observable.stem + "_lo.pdf"), observable.lo_title, lo_overlay);
            }

            if (distribution.config.write_nlo && distribution.has_nlo_histogram) {
                HistogramSeries nlo_reference;
                const HistogramSeries* nlo_overlay = nullptr;
                if (!distribution.config.reference_nlo_csv.empty() && std::filesystem::exists(distribution.config.reference_nlo_csv)) {
                    nlo_reference = load_histogram_csv(distribution.config.reference_nlo_csv, "reference_nlo", observable.x_label, observable.y_label);
                    nlo_overlay = &nlo_reference;
                }
                write_histogram_plot_png(distribution.nlo_histogram, png_dir / (observable.stem + "_nlo.png"), observable.nlo_title, nlo_overlay);
                write_histogram_plot_pdf(distribution.nlo_histogram, pdf_dir / (observable.stem + "_nlo.pdf"), observable.nlo_title, nlo_overlay);
            }
        }
    }
}

}  // namespace cs_ppttb