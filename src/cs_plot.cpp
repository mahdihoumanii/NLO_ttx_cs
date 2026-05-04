#include "cs_plot.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cs_ppttb {
namespace {

struct Color {
    std::uint8_t r;
    std::uint8_t g;
    std::uint8_t b;
};

constexpr int kWidth = 1200;
constexpr int kHeight = 800;
constexpr int kLeftMargin = 90;
constexpr int kRightMargin = 40;
constexpr int kTopMargin = 50;
constexpr int kBottomMargin = 70;

constexpr Color kWhite{255, 255, 255};
constexpr Color kGrid{220, 226, 232};
constexpr Color kAxis{60, 72, 88};
constexpr Color kPrimary{27, 73, 101};
constexpr Color kOverlay{220, 94, 30};

struct Image {
    int width = kWidth;
    int height = kHeight;
    std::vector<std::uint8_t> pixels = std::vector<std::uint8_t>(static_cast<std::size_t>(kWidth * kHeight * 3), 255);

    void set_pixel(int x, int y, Color color) {
        if (x < 0 || x >= width || y < 0 || y >= height) {
            return;
        }
        const std::size_t index = static_cast<std::size_t>((y * width + x) * 3);
        pixels[index] = color.r;
        pixels[index + 1] = color.g;
        pixels[index + 2] = color.b;
    }

    void fill_rect(int x0, int y0, int x1, int y1, Color color) {
        if (x0 > x1) {
            std::swap(x0, x1);
        }
        if (y0 > y1) {
            std::swap(y0, y1);
        }
        for (int y = y0; y <= y1; ++y) {
            for (int x = x0; x <= x1; ++x) {
                set_pixel(x, y, color);
            }
        }
    }

    void draw_line(int x0, int y0, int x1, int y1, Color color) {
        const int dx = std::abs(x1 - x0);
        const int sx = x0 < x1 ? 1 : -1;
        const int dy = -std::abs(y1 - y0);
        const int sy = y0 < y1 ? 1 : -1;
        int error = dx + dy;
        while (true) {
            set_pixel(x0, y0, color);
            if (x0 == x1 && y0 == y1) {
                break;
            }
            const int e2 = 2 * error;
            if (e2 >= dy) {
                error += dy;
                x0 += sx;
            }
            if (e2 <= dx) {
                error += dx;
                y0 += sy;
            }
        }
    }
};

void append_uint32_be(std::vector<std::uint8_t>& out, std::uint32_t value) {
    out.push_back(static_cast<std::uint8_t>((value >> 24) & 0xffu));
    out.push_back(static_cast<std::uint8_t>((value >> 16) & 0xffu));
    out.push_back(static_cast<std::uint8_t>((value >> 8) & 0xffu));
    out.push_back(static_cast<std::uint8_t>(value & 0xffu));
}

std::uint32_t crc32(const std::vector<std::uint8_t>& bytes) {
    std::uint32_t crc = 0xffffffffu;
    for (std::uint8_t byte : bytes) {
        crc ^= byte;
        for (int bit = 0; bit < 8; ++bit) {
            const std::uint32_t mask = static_cast<std::uint32_t>(-(static_cast<int>(crc & 1u)));
            crc = (crc >> 1) ^ (0xedb88320u & mask);
        }
    }
    return ~crc;
}

std::uint32_t adler32(const std::vector<std::uint8_t>& bytes) {
    constexpr std::uint32_t kMod = 65521u;
    std::uint32_t a = 1u;
    std::uint32_t b = 0u;
    for (std::uint8_t byte : bytes) {
        a = (a + byte) % kMod;
        b = (b + a) % kMod;
    }
    return (b << 16) | a;
}

void append_chunk(std::vector<std::uint8_t>& png,
                  const std::array<char, 4>& type,
                  const std::vector<std::uint8_t>& data) {
    append_uint32_be(png, static_cast<std::uint32_t>(data.size()));
    std::vector<std::uint8_t> crc_input;
    crc_input.reserve(4 + data.size());
    for (char ch : type) {
        const auto byte = static_cast<std::uint8_t>(ch);
        png.push_back(byte);
        crc_input.push_back(byte);
    }
    png.insert(png.end(), data.begin(), data.end());
    crc_input.insert(crc_input.end(), data.begin(), data.end());
    append_uint32_be(png, crc32(crc_input));
}

std::vector<std::uint8_t> build_png_bytes(const Image& image) {
    std::vector<std::uint8_t> raw;
    raw.reserve(static_cast<std::size_t>((image.width * 3 + 1) * image.height));
    for (int y = 0; y < image.height; ++y) {
        raw.push_back(0u);
        const std::size_t offset = static_cast<std::size_t>(y * image.width * 3);
        raw.insert(raw.end(), image.pixels.begin() + static_cast<std::ptrdiff_t>(offset), image.pixels.begin() + static_cast<std::ptrdiff_t>(offset + image.width * 3));
    }

    std::vector<std::uint8_t> zlib;
    zlib.push_back(0x78u);
    zlib.push_back(0x01u);
    std::size_t cursor = 0;
    while (cursor < raw.size()) {
        const std::size_t remaining = raw.size() - cursor;
        const std::size_t block_size = std::min<std::size_t>(remaining, 65535u);
        const bool final_block = (cursor + block_size == raw.size());
        zlib.push_back(final_block ? 0x01u : 0x00u);
        const std::uint16_t len = static_cast<std::uint16_t>(block_size);
        const std::uint16_t nlen = static_cast<std::uint16_t>(~len);
        zlib.push_back(static_cast<std::uint8_t>(len & 0xffu));
        zlib.push_back(static_cast<std::uint8_t>((len >> 8) & 0xffu));
        zlib.push_back(static_cast<std::uint8_t>(nlen & 0xffu));
        zlib.push_back(static_cast<std::uint8_t>((nlen >> 8) & 0xffu));
        zlib.insert(zlib.end(), raw.begin() + static_cast<std::ptrdiff_t>(cursor), raw.begin() + static_cast<std::ptrdiff_t>(cursor + block_size));
        cursor += block_size;
    }
    append_uint32_be(zlib, adler32(raw));

    std::vector<std::uint8_t> png = {
        0x89u, 0x50u, 0x4eu, 0x47u, 0x0du, 0x0au, 0x1au, 0x0au,
    };

    std::vector<std::uint8_t> ihdr;
    append_uint32_be(ihdr, static_cast<std::uint32_t>(image.width));
    append_uint32_be(ihdr, static_cast<std::uint32_t>(image.height));
    ihdr.push_back(8u);
    ihdr.push_back(2u);
    ihdr.push_back(0u);
    ihdr.push_back(0u);
    ihdr.push_back(0u);
    append_chunk(png, {'I', 'H', 'D', 'R'}, ihdr);
    append_chunk(png, {'I', 'D', 'A', 'T'}, zlib);
    append_chunk(png, {'I', 'E', 'N', 'D'}, {});
    return png;
}

void write_binary(const std::filesystem::path& path, const std::vector<std::uint8_t>& bytes) {
    std::ofstream output(path, std::ios::binary);
    if (!output) {
        throw std::runtime_error("Failed to open plot output " + path.string());
    }
    output.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
}

struct PlotBounds {
    double x_min = 0.0;
    double x_max = 1.0;
    double y_min = 0.0;
    double y_max = 1.0;
};

PlotBounds compute_bounds(const HistogramSeries& primary, const HistogramSeries* overlay) {
    if (primary.bins.empty()) {
        throw std::runtime_error("Cannot plot an empty histogram");
    }

    PlotBounds bounds;
    bounds.x_min = primary.bins.front().low;
    bounds.x_max = primary.bins.back().high;
    bounds.y_min = 0.0;
    bounds.y_max = 0.0;

    for (const HistogramBin& bin : primary.bins) {
        bounds.y_min = std::min(bounds.y_min, bin.value - bin.error);
        bounds.y_max = std::max(bounds.y_max, bin.value + bin.error);
    }
    if (overlay != nullptr) {
        for (const HistogramBin& bin : overlay->bins) {
            bounds.y_min = std::min(bounds.y_min, bin.value - bin.error);
            bounds.y_max = std::max(bounds.y_max, bin.value + bin.error);
        }
    }
    if (bounds.y_max <= bounds.y_min) {
        bounds.y_max = bounds.y_min + 1.0;
    }
    const double padding = 0.08 * (bounds.y_max - bounds.y_min);
    bounds.y_min -= padding;
    bounds.y_max += padding;
    return bounds;
}

double x_map(double x, const PlotBounds& bounds) {
    const double width = static_cast<double>(kWidth - kLeftMargin - kRightMargin);
    return static_cast<double>(kLeftMargin) + (x - bounds.x_min) / (bounds.x_max - bounds.x_min) * width;
}

double y_map(double y, const PlotBounds& bounds) {
    const double height = static_cast<double>(kHeight - kTopMargin - kBottomMargin);
    return static_cast<double>(kHeight - kBottomMargin) - (y - bounds.y_min) / (bounds.y_max - bounds.y_min) * height;
}

std::string pdf_escape(const std::string& text) {
    std::string escaped;
    escaped.reserve(text.size());
    for (char ch : text) {
        if (ch == '(' || ch == ')' || ch == '\\') {
            escaped.push_back('\\');
        }
        escaped.push_back(ch);
    }
    return escaped;
}

}  // namespace

void write_histogram_plot_png(const HistogramSeries& primary,
                              const std::filesystem::path& path,
                              const std::string& title,
                              const HistogramSeries* overlay) {
    (void)title;
    const PlotBounds bounds = compute_bounds(primary, overlay);
    Image image;
    image.fill_rect(0, 0, kWidth - 1, kHeight - 1, kWhite);

    const int plot_left = kLeftMargin;
    const int plot_right = kWidth - kRightMargin;
    const int plot_top = kTopMargin;
    const int plot_bottom = kHeight - kBottomMargin;
    const int zero_y = static_cast<int>(std::lround(y_map(0.0, bounds)));

    for (int tick = 0; tick <= 5; ++tick) {
        const double fraction = static_cast<double>(tick) / 5.0;
        const int y = static_cast<int>(std::lround(plot_top + fraction * static_cast<double>(plot_bottom - plot_top)));
        image.draw_line(plot_left, y, plot_right, y, kGrid);
    }

    image.draw_line(plot_left, plot_top, plot_left, plot_bottom, kAxis);
    image.draw_line(plot_left, zero_y, plot_right, zero_y, kAxis);

    for (const HistogramBin& bin : primary.bins) {
        const int x0 = static_cast<int>(std::lround(x_map(bin.low, bounds)));
        const int x1 = static_cast<int>(std::lround(x_map(bin.high, bounds)));
        const int y0 = static_cast<int>(std::lround(y_map(0.0, bounds)));
        const int y1 = static_cast<int>(std::lround(y_map(bin.value, bounds)));
        image.fill_rect(x0 + 1, std::min(y0, y1), x1 - 1, std::max(y0, y1), kPrimary);
    }

    if (overlay != nullptr && !overlay->bins.empty()) {
        for (std::size_t index = 1; index < overlay->bins.size(); ++index) {
            const double x_prev = 0.5 * (overlay->bins[index - 1].low + overlay->bins[index - 1].high);
            const double x_curr = 0.5 * (overlay->bins[index].low + overlay->bins[index].high);
            const double y_prev = overlay->bins[index - 1].value;
            const double y_curr = overlay->bins[index].value;
            image.draw_line(static_cast<int>(std::lround(x_map(x_prev, bounds))),
                            static_cast<int>(std::lround(y_map(y_prev, bounds))),
                            static_cast<int>(std::lround(x_map(x_curr, bounds))),
                            static_cast<int>(std::lround(y_map(y_curr, bounds))),
                            kOverlay);
        }
    }

    write_binary(path, build_png_bytes(image));
}

void write_histogram_plot_pdf(const HistogramSeries& primary,
                              const std::filesystem::path& path,
                              const std::string& title,
                              const HistogramSeries* overlay) {
    const PlotBounds bounds = compute_bounds(primary, overlay);
    std::ostringstream content;
    content.setf(std::ios::fixed);
    content.precision(3);

    content << "0.862 0.886 0.910 RG\n";
    for (int tick = 0; tick <= 5; ++tick) {
        const double fraction = static_cast<double>(tick) / 5.0;
        const double y = static_cast<double>(kTopMargin) + fraction * static_cast<double>(kHeight - kTopMargin - kBottomMargin);
        content << kLeftMargin << ' ' << y << " m " << (kWidth - kRightMargin) << ' ' << y << " l S\n";
    }

    const double zero_y = y_map(0.0, bounds);
    content << "0.235 0.282 0.345 RG\n";
    content << kLeftMargin << ' ' << kTopMargin << " m " << kLeftMargin << ' ' << (kHeight - kBottomMargin) << " l S\n";
    content << kLeftMargin << ' ' << zero_y << " m " << (kWidth - kRightMargin) << ' ' << zero_y << " l S\n";

    content << "0.106 0.286 0.396 rg\n";
    for (const HistogramBin& bin : primary.bins) {
        const double x0 = x_map(bin.low, bounds);
        const double x1 = x_map(bin.high, bounds);
        const double bar_y = y_map(bin.value, bounds);
        const double baseline = y_map(0.0, bounds);
        const double left = x0 + 1.0;
        const double width = std::max(0.0, x1 - x0 - 2.0);
        const double low_y = std::min(bar_y, baseline);
        const double height = std::abs(bar_y - baseline);
        content << left << ' ' << low_y << ' ' << width << ' ' << height << " re f\n";
    }

    if (overlay != nullptr && overlay->bins.size() >= 2) {
        content << "0.863 0.369 0.118 RG\n1.5 w\n";
        for (std::size_t index = 1; index < overlay->bins.size(); ++index) {
            const double x_prev = 0.5 * (overlay->bins[index - 1].low + overlay->bins[index - 1].high);
            const double x_curr = 0.5 * (overlay->bins[index].low + overlay->bins[index].high);
            const double y_prev = overlay->bins[index - 1].value;
            const double y_curr = overlay->bins[index].value;
            content << x_map(x_prev, bounds) << ' ' << y_map(y_prev, bounds)
                    << " m " << x_map(x_curr, bounds) << ' ' << y_map(y_curr, bounds) << " l S\n";
        }
    }

    content << "BT /F1 18 Tf 60 770 Td (" << pdf_escape(title) << ") Tj ET\n";
    content << "BT /F1 11 Tf 60 748 Td (" << pdf_escape(primary.y_label + " vs " + primary.x_label) << ") Tj ET\n";

    const std::string content_string = content.str();

    std::ostringstream output;
    output << "%PDF-1.4\n";

    std::vector<long> offsets;
    auto write_object = [&](int id, const std::string& body) {
        offsets.push_back(static_cast<long>(output.tellp()));
        output << id << " 0 obj\n" << body << "\nendobj\n";
    };

    write_object(1, "<< /Type /Catalog /Pages 2 0 R >>");
    write_object(2, "<< /Type /Pages /Kids [3 0 R] /Count 1 >>");
    write_object(3,
                 "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 1200 800] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>");
    write_object(4, "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>");
    write_object(5, "<< /Length " + std::to_string(content_string.size()) + " >>\nstream\n" + content_string + "endstream");

    const long xref_offset = static_cast<long>(output.tellp());
    output << "xref\n0 6\n";
    output << "0000000000 65535 f \n";
    for (long offset : offsets) {
        output << std::setw(10) << std::setfill('0') << offset << " 00000 n \n";
    }
    output << "trailer\n<< /Size 6 /Root 1 0 R >>\nstartxref\n" << xref_offset << "\n%%EOF\n";

    std::ofstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Failed to open plot output " + path.string());
    }
    file << output.str();
}

}  // namespace cs_ppttb