#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace cs_ppttb {

struct HistogramBin {
    double low = 0.0;
    double high = 0.0;
    double value = 0.0;
    double error = 0.0;
};

struct HistogramSeries {
    std::string name;
    std::string x_label;
    std::string y_label;
    std::vector<HistogramBin> bins;
};

class HistogramAccumulator {
public:
    HistogramAccumulator(int bin_count, double min_x, double max_x)
        : bin_count_(bin_count),
          min_x_(min_x),
          max_x_(max_x),
          bin_width_((max_x - min_x) / static_cast<double>(bin_count)),
          sums_(static_cast<std::size_t>(bin_count), 0.0),
          sums2_(static_cast<std::size_t>(bin_count), 0.0) {
        if (bin_count_ <= 0 || !(max_x_ > min_x_)) {
            throw std::invalid_argument("Invalid histogram definition");
        }
    }

    void observe(double x, double weight) {
        ++samples_;
        if (!std::isfinite(x) || !std::isfinite(weight)) {
            return;
        }
        if (x < min_x_ || x >= max_x_) {
            return;
        }
        int index = static_cast<int>((x - min_x_) / bin_width_);
        if (index < 0) {
            return;
        }
        if (index >= bin_count_) {
            index = bin_count_ - 1;
        }
        const double density_weight = weight / bin_width_;
        sums_[static_cast<std::size_t>(index)] += density_weight;
        sums2_[static_cast<std::size_t>(index)] += density_weight * density_weight;
    }

    HistogramSeries finalize(const std::string& name,
                             const std::string& x_label,
                             const std::string& y_label) const {
        HistogramSeries series;
        series.name = name;
        series.x_label = x_label;
        series.y_label = y_label;
        series.bins.reserve(static_cast<std::size_t>(bin_count_));

        const double n = static_cast<double>(samples_);
        for (int index = 0; index < bin_count_; ++index) {
            const double low = min_x_ + static_cast<double>(index) * bin_width_;
            const double high = low + bin_width_;
            const double sum = sums_[static_cast<std::size_t>(index)];
            const double sum2 = sums2_[static_cast<std::size_t>(index)];
            const double mean = samples_ > 0 ? sum / n : 0.0;
            double variance = 0.0;
            if (samples_ > 1) {
                variance = std::max(0.0, (sum2 - sum * sum / n) / (n - 1.0));
            }
            const double error = samples_ > 0 ? std::sqrt(variance / n) : 0.0;
            series.bins.push_back({low, high, mean, error});
        }
        return series;
    }

private:
    int bin_count_ = 0;
    long samples_ = 0;
    double min_x_ = 0.0;
    double max_x_ = 0.0;
    double bin_width_ = 0.0;
    std::vector<double> sums_;
    std::vector<double> sums2_;
};

inline HistogramSeries combine_histogram_series(const std::string& name,
                                                const std::string& x_label,
                                                const std::string& y_label,
                                                const std::vector<HistogramSeries>& parts) {
    if (parts.empty()) {
        return {name, x_label, y_label, {}};
    }
    HistogramSeries combined = parts.front();
    combined.name = name;
    combined.x_label = x_label;
    combined.y_label = y_label;
    for (std::size_t part_index = 1; part_index < parts.size(); ++part_index) {
        const HistogramSeries& part = parts[part_index];
        if (part.bins.size() != combined.bins.size()) {
            throw std::runtime_error("Histogram bin counts do not match");
        }
        for (std::size_t bin_index = 0; bin_index < combined.bins.size(); ++bin_index) {
            HistogramBin& target = combined.bins[bin_index];
            const HistogramBin& source = part.bins[bin_index];
            if (std::abs(target.low - source.low) > 1.0e-12 || std::abs(target.high - source.high) > 1.0e-12) {
                throw std::runtime_error("Histogram bin edges do not match");
            }
            target.value += source.value;
            target.error = std::sqrt(target.error * target.error + source.error * source.error);
        }
    }
    return combined;
}

}  // namespace cs_ppttb