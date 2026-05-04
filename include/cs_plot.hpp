#pragma once

#include <filesystem>
#include <string>

#include "cs_histogram.hpp"

namespace cs_ppttb {

void write_histogram_plot_png(const HistogramSeries& primary,
                              const std::filesystem::path& path,
                              const std::string& title,
                              const HistogramSeries* overlay = nullptr);

void write_histogram_plot_pdf(const HistogramSeries& primary,
                              const std::filesystem::path& path,
                              const std::string& title,
                              const HistogramSeries* overlay = nullptr);

}  // namespace cs_ppttb