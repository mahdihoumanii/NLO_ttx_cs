#pragma once

#include <filesystem>
#include <limits>
#include <string>
#include <vector>

namespace cs_ppttb {

struct DistributionConfig {
    std::string observable;
    int bins = 40;
    double min_gev = std::numeric_limits<double>::quiet_NaN();
    double max_gev = std::numeric_limits<double>::quiet_NaN();
    bool write_lo = true;
    bool write_nlo = true;
    std::filesystem::path ra_histogram_csv;
    std::filesystem::path reference_lo_csv;
    std::filesystem::path reference_nlo_csv;
};

struct RunConfig {
    std::filesystem::path config_path;
    std::filesystem::path package_root;
    std::filesystem::path output_dir;
    std::filesystem::path born_phase_space;
    std::filesystem::path ra_artifact;
    std::filesystem::path ra_histogram_csv;
    std::filesystem::path reference_lo_csv;
    std::filesystem::path reference_nlo_csv;

    double sqrt_s_gev = 8000.0;
    double mt_gev = 173.3;
    double muF_gev = 173.2;
    double muR_gev = 173.2;

    std::string pdf_set = "NNPDF31_nlo_as_0118";
    int pdf_member = 0;
    std::string ra_mode = "artifact";
    std::string observable = "mtt";

    long lo_samples = 1000000;
    long ca_samples = 1000000;
    long va_samples = 1000000;
    unsigned long long seed = 8675309ULL;

    int mtt_bins = 40;
    double mtt_min_gev = std::numeric_limits<double>::quiet_NaN();
    double mtt_max_gev = std::numeric_limits<double>::quiet_NaN();

    bool write_csv = true;
    bool write_plots = true;
    bool write_lo_distribution = true;
    bool write_nlo_distribution = true;
    std::vector<DistributionConfig> distributions;
};

RunConfig load_config(const std::filesystem::path& path);

}  // namespace cs_ppttb