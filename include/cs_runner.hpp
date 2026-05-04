#pragma once

#include <string>
#include <vector>

#include "cs_config.hpp"
#include "cs_histogram.hpp"

namespace cs_ppttb {

struct ComponentResult {
    std::string name;
    double value = 0.0;
    double error = 0.0;
    std::string source;
    bool artifact_only = false;
    bool has_histogram = false;
    HistogramSeries histogram;
};

struct DistributionResult {
    DistributionConfig config;
    HistogramSeries lo_histogram;
    HistogramSeries ra_histogram;
    HistogramSeries ca_histogram;
    HistogramSeries va_histogram;
    HistogramSeries nlo_histogram;
    bool has_ra_histogram = false;
    bool has_nlo_histogram = false;
};

struct RunReport {
    RunConfig config;
    std::vector<ComponentResult> components;
    std::vector<DistributionResult> distributions;
    double nlo_total = 0.0;
    double nlo_error = 0.0;
    bool has_nlo_histogram = false;
    HistogramSeries lo_histogram;
    HistogramSeries nlo_histogram;
    std::vector<std::string> notes;
};

RunReport run_package(const RunConfig& config);
void write_run_outputs(const RunReport& report);

}  // namespace cs_ppttb