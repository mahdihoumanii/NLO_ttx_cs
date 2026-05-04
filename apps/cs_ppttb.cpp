#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>

#include "cs_config.hpp"
#include "cs_runner.hpp"

namespace {

long parse_long(const std::string& value, const std::string& key) {
    try {
        std::size_t parsed = 0;
        const long result = std::stol(value, &parsed);
        if (parsed != value.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return result;
    } catch (const std::exception&) {
        throw std::runtime_error("Failed to parse " + key + ": " + value);
    }
}

unsigned long long parse_seed(const std::string& value) {
    try {
        std::size_t parsed = 0;
        const unsigned long long result = std::stoull(value, &parsed);
        if (parsed != value.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return result;
    } catch (const std::exception&) {
        throw std::runtime_error("Failed to parse seed: " + value);
    }
}

void print_help() {
    std::cout << "Usage: cs_ppttb [--config PATH] [--lo-samples N] [--ca-samples N] [--va-samples N] [--seed N] [--output-dir PATH]\n";
}

void print_run_summary(const cs_ppttb::RunReport& report) {
    constexpr double pb_per_fb = 1.0e-3;

    std::cout << std::endl;
    std::cout << "************************************************************" << std::endl;
    std::cout << "** cs_ppttb NLO cross section summary                     **" << std::endl;
    std::cout << "************************************************************" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "** sigma_NLO = " << report.nlo_total * pb_per_fb
              << " +/- " << report.nlo_error * pb_per_fb << " pb **" << std::endl;
    std::cout << "************************************************************" << std::endl;
    std::cout << std::endl;

    std::cout << "Component breakdown [pb]" << std::endl;
    for (const cs_ppttb::ComponentResult& component : report.components) {
        std::cout << "  " << std::setw(2) << component.name << " = "
                  << std::setw(12) << component.value * pb_per_fb
                  << " +/- " << std::setw(10) << component.error * pb_per_fb
                  << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Distribution output" << std::endl;
    for (const cs_ppttb::DistributionResult& distribution : report.distributions) {
        std::cout << "  " << distribution.config.observable
                  << "  LO " << (distribution.config.write_lo ? "enabled" : "disabled")
                  << "  NLO " << (distribution.config.write_nlo ? "enabled" : "disabled");
        if (distribution.config.write_nlo) {
            std::cout << "  RA histogram " << (distribution.has_ra_histogram ? "found" : "missing");
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Output directory: " << report.config.output_dir << std::endl;
    for (const std::string& note : report.notes) {
        std::cout << "note: " << note << std::endl;
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::string config_path = "config/default_8tev.cfg";
        std::string output_override;
        long lo_override = -1;
        long ca_override = -1;
        long va_override = -1;
        unsigned long long seed_override = 0;
        bool has_seed_override = false;

        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            const auto require_value = [&](const std::string& flag) -> std::string {
                if (index + 1 >= argc) {
                    throw std::runtime_error(flag + " requires a value");
                }
                return argv[++index];
            };
            const auto after_equals = [&](const std::string& prefix) {
                return argument.substr(prefix.size());
            };

            if (argument == "--config") {
                config_path = require_value(argument);
            } else if (argument.rfind("--config=", 0) == 0) {
                config_path = after_equals("--config=");
            } else if (argument == "--output-dir") {
                output_override = require_value(argument);
            } else if (argument.rfind("--output-dir=", 0) == 0) {
                output_override = after_equals("--output-dir=");
            } else if (argument == "--lo-samples") {
                lo_override = parse_long(require_value(argument), "lo-samples");
            } else if (argument.rfind("--lo-samples=", 0) == 0) {
                lo_override = parse_long(after_equals("--lo-samples="), "lo-samples");
            } else if (argument == "--ca-samples") {
                ca_override = parse_long(require_value(argument), "ca-samples");
            } else if (argument.rfind("--ca-samples=", 0) == 0) {
                ca_override = parse_long(after_equals("--ca-samples="), "ca-samples");
            } else if (argument == "--va-samples") {
                va_override = parse_long(require_value(argument), "va-samples");
            } else if (argument.rfind("--va-samples=", 0) == 0) {
                va_override = parse_long(after_equals("--va-samples="), "va-samples");
            } else if (argument == "--seed") {
                seed_override = parse_seed(require_value(argument));
                has_seed_override = true;
            } else if (argument.rfind("--seed=", 0) == 0) {
                seed_override = parse_seed(after_equals("--seed="));
                has_seed_override = true;
            } else if (argument == "--help" || argument == "-h") {
                print_help();
                return 0;
            } else {
                throw std::runtime_error("Unknown argument: " + argument);
            }
        }

        cs_ppttb::RunConfig config = cs_ppttb::load_config(config_path);
        if (!output_override.empty()) {
            config.output_dir = output_override;
        }
        if (lo_override > 0) {
            config.lo_samples = lo_override;
        }
        if (ca_override > 0) {
            config.ca_samples = ca_override;
        }
        if (va_override > 0) {
            config.va_samples = va_override;
        }
        if (has_seed_override) {
            config.seed = seed_override;
        }

        const cs_ppttb::RunReport report = cs_ppttb::run_package(config);
        cs_ppttb::write_run_outputs(report);
        print_run_summary(report);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "cs_ppttb: " << error.what() << '\n';
        return 1;
    }
}