#include "cs_config.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace cs_ppttb {
namespace {

std::string trim(const std::string& value) {
    const std::size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const std::size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::string lowercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

bool parse_bool(const std::string& value, const std::string& key) {
    const std::string normalized = lowercase(trim(value));
    if (normalized == "1" || normalized == "true" || normalized == "yes" || normalized == "on") {
        return true;
    }
    if (normalized == "0" || normalized == "false" || normalized == "no" || normalized == "off") {
        return false;
    }
    throw std::runtime_error("Failed to parse boolean for " + key + ": " + value);
}

long parse_long(const std::string& value, const std::string& key) {
    try {
        std::size_t parsed = 0;
        const long result = std::stol(value, &parsed);
        if (parsed != value.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return result;
    } catch (const std::exception&) {
        throw std::runtime_error("Failed to parse integer for " + key + ": " + value);
    }
}

int parse_int(const std::string& value, const std::string& key) {
    const long result = parse_long(value, key);
    return static_cast<int>(result);
}

unsigned long long parse_seed(const std::string& value, const std::string& key) {
    try {
        std::size_t parsed = 0;
        const unsigned long long result = std::stoull(value, &parsed);
        if (parsed != value.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return result;
    } catch (const std::exception&) {
        throw std::runtime_error("Failed to parse integer for " + key + ": " + value);
    }
}

double parse_double(const std::string& value, const std::string& key) {
    try {
        std::size_t parsed = 0;
        const double result = std::stod(value, &parsed);
        if (parsed != value.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return result;
    } catch (const std::exception&) {
        throw std::runtime_error("Failed to parse floating-point value for " + key + ": " + value);
    }
}

std::filesystem::path discover_package_root(const std::filesystem::path& config_path) {
    std::filesystem::path probe = std::filesystem::absolute(config_path).parent_path();
    while (!probe.empty()) {
        if (std::filesystem::exists(probe / "include") && std::filesystem::exists(probe / "src")) {
            return probe;
        }
        if (probe == probe.root_path()) {
            break;
        }
        probe = probe.parent_path();
    }
    throw std::runtime_error("Could not discover cs_ppttb package root from config path " + config_path.string());
}

std::filesystem::path resolve_path(const std::filesystem::path& root, const std::string& value) {
    if (trim(value).empty()) {
        return {};
    }
    const std::filesystem::path path_value(value);
    if (path_value.is_absolute()) {
        return path_value;
    }
    return root / path_value;
}

DistributionConfig default_distribution_config(const std::string& observable) {
    DistributionConfig distribution;
    distribution.observable = observable;
    distribution.write_lo = false;
    distribution.write_nlo = false;
    if (observable == "top_pt" || observable == "tb_pt") {
        distribution.bins = 100;
        distribution.min_gev = 0.0;
        distribution.max_gev = 500.0;
    } else if (observable == "ytt") {
        distribution.bins = 80;
        distribution.min_gev = -4.0;
        distribution.max_gev = 4.0;
    } else {
        distribution.observable = "mtt";
        distribution.bins = 40;
    }
    return distribution;
}

bool parse_distribution_key(const std::string& key, std::string& observable, std::string& field) {
    const std::string prefixes[] = {"top_pt", "tb_pt", "mtt", "ytt", "y_tt"};
    for (const std::string& prefix : prefixes) {
        const std::string prefix_with_separator = prefix + "_";
        if (key.rfind(prefix_with_separator, 0) == 0) {
            observable = (prefix == "y_tt") ? "ytt" : prefix;
            field = key.substr(prefix_with_separator.size());
            return true;
        }
    }
    return false;
}

DistributionConfig single_legacy_distribution(const RunConfig& config) {
    DistributionConfig distribution;
    distribution.observable = config.observable;
    distribution.bins = config.mtt_bins;
    distribution.min_gev = config.mtt_min_gev;
    distribution.max_gev = config.mtt_max_gev;
    distribution.write_lo = config.write_lo_distribution;
    distribution.write_nlo = config.write_nlo_distribution;
    distribution.ra_histogram_csv = config.ra_histogram_csv;
    distribution.reference_lo_csv = config.reference_lo_csv;
    distribution.reference_nlo_csv = config.reference_nlo_csv;
    return distribution;
}

}  // namespace

RunConfig load_config(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Failed to open config file " + path.string());
    }

    std::map<std::string, std::string> entries;
    std::string line;
    int line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) {
            line = line.substr(0, comment);
        }
        line = trim(line);
        if (line.empty()) {
            continue;
        }
        const std::size_t equals = line.find('=');
        if (equals == std::string::npos) {
            throw std::runtime_error("Malformed config line " + std::to_string(line_number) + ": " + line);
        }
        const std::string key = trim(line.substr(0, equals));
        const std::string value = trim(line.substr(equals + 1));
        if (key.empty()) {
            throw std::runtime_error("Malformed config line " + std::to_string(line_number) + ": empty key");
        }
        entries[key] = value;
    }

    RunConfig config;
    config.config_path = std::filesystem::absolute(path);
    config.package_root = discover_package_root(config.config_path);

    std::map<std::string, DistributionConfig> distribution_options{
        {"mtt", default_distribution_config("mtt")},
        {"top_pt", default_distribution_config("top_pt")},
        {"tb_pt", default_distribution_config("tb_pt")},
        {"ytt", default_distribution_config("ytt")},
    };
    bool has_distribution_options = false;

    for (const auto& entry : entries) {
        const std::string& key = entry.first;
        const std::string& value = entry.second;
        if (key == "sqrt_s_gev") {
            config.sqrt_s_gev = parse_double(value, key);
        } else if (key == "mt_gev") {
            config.mt_gev = parse_double(value, key);
        } else if (key == "muF_gev") {
            config.muF_gev = parse_double(value, key);
        } else if (key == "muR_gev") {
            config.muR_gev = parse_double(value, key);
        } else if (key == "pdf_set") {
            config.pdf_set = value;
        } else if (key == "pdf_member") {
            config.pdf_member = parse_int(value, key);
        } else if (key == "ra_mode") {
            config.ra_mode = value;
        } else if (key == "observable") {
            config.observable = value;
        } else if (key == "lo_samples") {
            config.lo_samples = parse_long(value, key);
        } else if (key == "ca_samples") {
            config.ca_samples = parse_long(value, key);
        } else if (key == "va_samples") {
            config.va_samples = parse_long(value, key);
        } else if (key == "seed") {
            config.seed = parse_seed(value, key);
        } else if (key == "mtt_bins") {
            config.mtt_bins = parse_int(value, key);
            distribution_options["mtt"].bins = config.mtt_bins;
        } else if (key == "mtt_min_gev") {
            config.mtt_min_gev = parse_double(value, key);
            distribution_options["mtt"].min_gev = config.mtt_min_gev;
        } else if (key == "mtt_max_gev") {
            config.mtt_max_gev = parse_double(value, key);
            distribution_options["mtt"].max_gev = config.mtt_max_gev;
        } else if (key == "write_csv") {
            config.write_csv = parse_bool(value, key);
        } else if (key == "write_plots") {
            config.write_plots = parse_bool(value, key);
        } else if (key == "write_lo_distribution" || key == "mtt_lo") {
            config.write_lo_distribution = parse_bool(value, key);
            distribution_options["mtt"].write_lo = config.write_lo_distribution;
            has_distribution_options = true;
        } else if (key == "write_nlo_distribution" || key == "mtt_nlo") {
            config.write_nlo_distribution = parse_bool(value, key);
            distribution_options["mtt"].write_nlo = config.write_nlo_distribution;
            has_distribution_options = true;
        } else if (key == "output_dir") {
            config.output_dir = resolve_path(config.package_root, value);
        } else if (key == "born_phase_space") {
            config.born_phase_space = resolve_path(config.package_root, value);
        } else if (key == "ra_artifact") {
            config.ra_artifact = resolve_path(config.package_root, value);
        } else if (key == "ra_histogram_csv") {
            config.ra_histogram_csv = resolve_path(config.package_root, value);
            distribution_options["mtt"].ra_histogram_csv = config.ra_histogram_csv;
        } else if (key == "reference_lo_csv") {
            config.reference_lo_csv = resolve_path(config.package_root, value);
            distribution_options["mtt"].reference_lo_csv = config.reference_lo_csv;
        } else if (key == "reference_nlo_csv") {
            config.reference_nlo_csv = resolve_path(config.package_root, value);
            distribution_options["mtt"].reference_nlo_csv = config.reference_nlo_csv;
        } else {
            std::string distribution_observable;
            std::string distribution_field;
            if (!parse_distribution_key(key, distribution_observable, distribution_field)) {
                throw std::runtime_error("Unknown config key: " + key);
            }
            DistributionConfig& distribution = distribution_options[distribution_observable];
            has_distribution_options = true;
            if (distribution_field == "lo") {
                distribution.write_lo = parse_bool(value, key);
            } else if (distribution_field == "nlo") {
                distribution.write_nlo = parse_bool(value, key);
            } else if (distribution_field == "bins") {
                distribution.bins = parse_int(value, key);
            } else if (distribution_field == "min_gev" || distribution_field == "min") {
                distribution.min_gev = parse_double(value, key);
            } else if (distribution_field == "max_gev" || distribution_field == "max") {
                distribution.max_gev = parse_double(value, key);
            } else if (distribution_field == "ra_histogram_csv") {
                distribution.ra_histogram_csv = resolve_path(config.package_root, value);
            } else if (distribution_field == "reference_lo_csv") {
                distribution.reference_lo_csv = resolve_path(config.package_root, value);
            } else if (distribution_field == "reference_nlo_csv") {
                distribution.reference_nlo_csv = resolve_path(config.package_root, value);
            } else {
                throw std::runtime_error("Unknown distribution config key: " + key);
            }
        }
    }

    if (config.output_dir.empty()) {
        config.output_dir = config.package_root / "results" / "8tev_run";
    }
    if (config.born_phase_space.empty()) {
        config.born_phase_space = config.package_root / "data" / "references" / "ppttx20.phasespace.born.dat";
    }
    if (config.ra_artifact.empty()) {
        config.ra_artifact = config.package_root / "validation" / "local_ra_integral_all_channels_100k_refresh.txt";
    }

    if (config.sqrt_s_gev <= 0.0) {
        throw std::runtime_error("sqrt_s_gev must be positive");
    }
    if (config.mt_gev <= 0.0 || config.muF_gev <= 0.0 || config.muR_gev <= 0.0) {
        throw std::runtime_error("Mass and scale parameters must be positive");
    }
    if (config.lo_samples <= 0 || config.ca_samples <= 0 || config.va_samples <= 0) {
        throw std::runtime_error("All sample counts must be positive");
    }
    if (config.observable.empty()) {
        throw std::runtime_error("observable must not be empty");
    }
    if (config.mtt_bins <= 0) {
        throw std::runtime_error("mtt_bins must be positive");
    }
    if (std::isnan(config.mtt_min_gev)) {
        config.mtt_min_gev = 2.0 * config.mt_gev;
    }
    if (std::isnan(config.mtt_max_gev)) {
        config.mtt_max_gev = config.sqrt_s_gev;
    }
    if (config.mtt_min_gev >= config.mtt_max_gev) {
        throw std::runtime_error("mtt histogram range is invalid");
    }

    if (!has_distribution_options) {
        config.distributions.push_back(single_legacy_distribution(config));
    } else {
        const std::string ordered_observables[] = {"mtt", "top_pt", "tb_pt", "ytt"};
        for (const std::string& observable_key : ordered_observables) {
            DistributionConfig distribution = distribution_options[observable_key];
            if (std::isnan(distribution.min_gev)) {
                distribution.min_gev = 2.0 * config.mt_gev;
            }
            if (std::isnan(distribution.max_gev)) {
                distribution.max_gev = config.sqrt_s_gev;
            }
            if (distribution.bins <= 0) {
                throw std::runtime_error(distribution.observable + "_bins must be positive");
            }
            if (distribution.min_gev >= distribution.max_gev) {
                throw std::runtime_error(distribution.observable + " histogram range is invalid");
            }
            config.distributions.push_back(distribution);
        }
    }

    return config;
}

}  // namespace cs_ppttb