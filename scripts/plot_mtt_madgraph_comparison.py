#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import AutoMinorLocator, MaxNLocator, MultipleLocator


MG_LINE_COLOR = "#1E5A8A"
MG_BAND_COLOR = "#A9C7E8"
MG_STAT_COLOR = "#D8E7F5"
LOCAL_LINE_COLOR = "#B43C4A"
LOCAL_BAND_COLOR = "#F0B7BF"
FRAME_COLOR = "#2D2A26"
GRID_COLOR = "#D9D3C8"
TEXT_MUTED_COLOR = "#655C52"


mpl.rcParams.update(
    {
        "figure.dpi": 160,
        "savefig.dpi": 300,
        "font.family": "DejaVu Sans",
        "axes.labelsize": 18,
        "axes.titlesize": 22,
        "axes.linewidth": 1.0,
        "axes.edgecolor": FRAME_COLOR,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
        "legend.fontsize": 12.5,
        "axes.unicode_minus": False,
    }
)


@dataclass
class BinnedHistogram:
    name: str
    edges: list[float]
    values: list[float]
    errors: list[float]
    band_low: list[float]
    band_high: list[float]
    band_label: str = "Scale band"


@dataclass(frozen=True)
class ObservableSpec:
    key: str
    mg_title: str
    stem: str
    axis_label: str
    x_major_step: float


OBSERVABLES = {
    "mtt": ObservableSpec("mtt", "tt inv m", "mtt", "mtt [GeV]", 100.0),
    "top_pt": ObservableSpec("top_pt", "t pt", "top_pt", "pT(t) [GeV]", 50.0),
    "tb_pt": ObservableSpec("tb_pt", "tb pt", "tb_pt", "pT(tbar) [GeV]", 50.0),
    "ytt": ObservableSpec("ytt", "y_tt", "ytt", "y_tt", 1.0),
}


def observable_spec(key: str) -> ObservableSpec:
    normalized = key.strip().lower()
    aliases = {
        "t_pt": "top_pt",
        "tpt": "top_pt",
        "antitop_pt": "tb_pt",
        "tbar_pt": "tb_pt",
        "tbpt": "tb_pt",
        "y_tt": "ytt",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in OBSERVABLES:
        raise RuntimeError(f"Unsupported observable key: {key}")
    return OBSERVABLES[normalized]


def step_values(values: list[float]) -> list[float]:
    if not values:
        return []
    return values + [values[-1]]


def positive_step_values(values: list[float]) -> list[float]:
    if not values:
        return []
    masked = [value if math.isfinite(value) and value > 0.0 else float("nan") for value in values]
    return masked + [masked[-1]]


def clipped_positive_step_values(values: list[float], floor: float) -> list[float]:
    if not values:
        return []
    clipped = [max(value, floor) if math.isfinite(value) else float("nan") for value in values]
    return clipped + [clipped[-1]]


def format_pdf_label(pdf_set: str) -> str:
    if pdf_set == "NNPDF31_nlo_as_0118":
        return r"NNPDF3.1 NLO, $\alpha_s(m_Z)=0.118$"
    return pdf_set.replace("_", " ")


def format_parameter_value(value: float) -> str:
    rounded = round(value)
    if math.isclose(value, rounded, rel_tol=0.0, abs_tol=1.0e-12):
        return str(int(rounded))
    return f"{value:.1f}"


def normalize_histogram_title(title: str) -> str:
    return " ".join(title.split())


def parse_hwu_histogram(path: Path, observable: str, allow_cuts: bool = False) -> BinnedHistogram:
    expected_title = normalize_histogram_title(observable)
    lines = path.read_text().splitlines()
    for index, line in enumerate(lines):
        if not line.startswith("<histogram>"):
            continue
        raw_title = line.split('"', 2)[1]
        title = normalize_histogram_title(raw_title.split("|", 1)[0])
        if title != expected_title:
            continue
        if ("cuts" in normalize_histogram_title(raw_title).lower()) and not allow_cuts:
            continue

        edges: list[float] = []
        values: list[float] = []
        errors: list[float] = []
        band_low: list[float] = []
        band_high: list[float] = []

        cursor = index + 1
        while cursor < len(lines):
            current = lines[cursor].strip()
            if not current:
                cursor += 1
                continue
            if current.startswith("<"):
                break
            if current.startswith("#"):
                cursor += 1
                continue

            columns = [float(token) for token in current.split()]
            if len(columns) < 7:
                break

            if not edges:
                edges.append(columns[0])
            edges.append(columns[1])
            values.append(columns[2])
            errors.append(columns[3])

            scale_columns = columns[4:]
            band_low.append(min(scale_columns))
            band_high.append(max(scale_columns))
            cursor += 1

        if not values:
            raise RuntimeError(f"Histogram '{observable}' was found in {path} but had no data rows")

        return BinnedHistogram(title, edges, values, errors, band_low, band_high, "MadGraph scale band")

    raise RuntimeError(f"Could not find histogram '{observable}' in {path}")


def parse_local_csv(path: Path, name: str) -> BinnedHistogram:
    edges: list[float] = []
    values: list[float] = []
    errors: list[float] = []
    band_low: list[float] = []
    band_high: list[float] = []
    band_label = "Local statistical band"

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        has_scale_band = reader.fieldnames is not None and "scale_low" in reader.fieldnames and "scale_high" in reader.fieldnames
        for row in reader:
            low = float(row["bin_low_gev"])
            high = float(row["bin_high_gev"])
            width = high - low
            value_fb_per_gev = float(row["value"])
            error_fb_per_gev = float(row["error"])

            if not edges:
                edges.append(low)
            edges.append(high)

            values.append(value_fb_per_gev * width / 1000.0)
            errors.append(error_fb_per_gev * width / 1000.0)
            if has_scale_band and row.get("scale_low") and row.get("scale_high"):
                band_low.append(float(row["scale_low"]) * width / 1000.0)
                band_high.append(float(row["scale_high"]) * width / 1000.0)
                band_label = "Local scale band"
            else:
                band_low.append((value_fb_per_gev - error_fb_per_gev) * width / 1000.0)
                band_high.append((value_fb_per_gev + error_fb_per_gev) * width / 1000.0)

    if not values:
        raise RuntimeError(f"Local histogram CSV is empty: {path}")

    return BinnedHistogram(name, edges, values, errors, band_low, band_high, band_label)


def rebin_integrated(source: BinnedHistogram, target_edges: list[float], name: str) -> BinnedHistogram:
    values: list[float] = []
    errors: list[float] = []
    band_low: list[float] = []
    band_high: list[float] = []

    for target_index in range(len(target_edges) - 1):
        target_low = target_edges[target_index]
        target_high = target_edges[target_index + 1]
        value = 0.0
        error2 = 0.0
        lower = 0.0
        upper = 0.0

        for source_index in range(len(source.values)):
            source_low = source.edges[source_index]
            source_high = source.edges[source_index + 1]
            overlap = max(0.0, min(target_high, source_high) - max(target_low, source_low))
            if overlap <= 0.0:
                continue

            source_width = source_high - source_low
            if source_width <= 0.0:
                continue

            fraction = overlap / source_width
            value += source.values[source_index] * fraction
            error2 += (source.errors[source_index] * fraction) ** 2
            lower += source.band_low[source_index] * fraction
            upper += source.band_high[source_index] * fraction

        values.append(value)
        errors.append(math.sqrt(error2))
        band_low.append(lower)
        band_high.append(upper)

    return BinnedHistogram(name, target_edges, values, errors, band_low, band_high)


def overlapping_range(first: BinnedHistogram, second: BinnedHistogram) -> tuple[float, float]:
    low = max(first.edges[0], second.edges[0])
    high = min(first.edges[-1], second.edges[-1])
    if high <= low:
        raise RuntimeError("The local and MadGraph histograms do not overlap")
    return low, high


def trim_histogram(histogram: BinnedHistogram, low: float, high: float, name: str) -> BinnedHistogram:
    return rebin_integrated(histogram, clipped_edges(histogram.edges, low, high), name)


def clipped_edges(edges: list[float], low: float, high: float) -> list[float]:
    clipped: list[float] = [low]
    for edge in edges:
        if low < edge < high:
            clipped.append(edge)
    clipped.append(high)
    clipped = sorted(set(clipped))
    if len(clipped) < 2:
        raise RuntimeError("No bins remain after clipping the histogram range")
    return clipped


def histogram_has_positive_content(histogram: BinnedHistogram) -> bool:
    return any(value > 0.0 for value in histogram.values)


def positive_floor(*histograms: BinnedHistogram) -> float:
    positives = [value for histogram in histograms for value in histogram.values if value > 0.0]
    if not positives:
        return 1.0e-3
    return min(positives)


def ratio_histogram(numerator: BinnedHistogram, denominator: BinnedHistogram, name: str) -> BinnedHistogram:
    rebinned_denominator = rebin_integrated(denominator, numerator.edges, f"{name}_denominator")
    values: list[float] = []
    errors: list[float] = []
    band_low: list[float] = []
    band_high: list[float] = []

    for num_value, num_error, denom_value in zip(numerator.values, numerator.errors, rebinned_denominator.values):
        if denom_value == 0.0:
            values.append(float("nan"))
            errors.append(float("nan"))
            band_low.append(float("nan"))
            band_high.append(float("nan"))
            continue

        ratio = num_value / denom_value
        ratio_error = num_error / abs(denom_value)
        values.append(ratio)
        errors.append(ratio_error)
        band_low.append(ratio - ratio_error)
        band_high.append(ratio + ratio_error)

    return BinnedHistogram(name, numerator.edges, values, errors, band_low, band_high, numerator.band_label)


def ratio_histogram_with_band(numerator: BinnedHistogram, denominator: BinnedHistogram, name: str) -> BinnedHistogram:
    rebinned_denominator = rebin_integrated(denominator, numerator.edges, f"{name}_denominator")
    values: list[float] = []
    errors: list[float] = []
    band_low: list[float] = []
    band_high: list[float] = []

    for num_value, num_error, num_low, num_high, denom_value in zip(
        numerator.values,
        numerator.errors,
        numerator.band_low,
        numerator.band_high,
        rebinned_denominator.values,
    ):
        if denom_value == 0.0 or not math.isfinite(denom_value):
            values.append(float("nan"))
            errors.append(float("nan"))
            band_low.append(float("nan"))
            band_high.append(float("nan"))
            continue

        ratio = num_value / denom_value
        ratio_error = num_error / abs(denom_value)
        values.append(ratio)
        errors.append(ratio_error)
        band_low.append(num_low / denom_value if math.isfinite(num_low) else float("nan"))
        band_high.append(num_high / denom_value if math.isfinite(num_high) else float("nan"))

    return BinnedHistogram(name, numerator.edges, values, errors, band_low, band_high, numerator.band_label)


def mg_ratio_band(histogram: BinnedHistogram, kind: str) -> BinnedHistogram:
    values = [1.0 for _ in histogram.values]
    errors = [0.0 for _ in histogram.values]
    band_low: list[float] = []
    band_high: list[float] = []

    for value, error, low, high in zip(histogram.values, histogram.errors, histogram.band_low, histogram.band_high):
        if value == 0.0:
            band_low.append(float("nan"))
            band_high.append(float("nan"))
            continue

        if kind == "scale":
            band_low.append(low / value)
            band_high.append(high / value)
        elif kind == "stat":
            band_low.append((value - error) / value)
            band_high.append((value + error) / value)
        else:
            raise RuntimeError(f"Unsupported ratio band type: {kind}")

    return BinnedHistogram(f"{histogram.name}_{kind}", histogram.edges, values, errors, band_low, band_high, histogram.band_label)


def mask_unstable_ratio_bins(reference: BinnedHistogram,
                            histogram: BinnedHistogram,
                            minimum_fraction: float = 0.01) -> BinnedHistogram:
    peak = max((abs(value) for value in reference.values if math.isfinite(value)), default=0.0)
    threshold = peak * minimum_fraction

    masked_values: list[float] = []
    masked_errors: list[float] = []
    masked_low: list[float] = []
    masked_high: list[float] = []

    for reference_value, value, error, low, high in zip(
        reference.values, histogram.values, histogram.errors, histogram.band_low, histogram.band_high
    ):
        stable = math.isfinite(reference_value) and reference_value > 0.0 and abs(reference_value) >= threshold
        if stable:
            masked_values.append(value)
            masked_errors.append(error)
            masked_low.append(low)
            masked_high.append(high)
        else:
            masked_values.append(float("nan"))
            masked_errors.append(float("nan"))
            masked_low.append(float("nan"))
            masked_high.append(float("nan"))

    return BinnedHistogram(histogram.name, histogram.edges, masked_values, masked_errors, masked_low, masked_high, histogram.band_label)


def finite_ratio_bounds(*histograms: BinnedHistogram) -> tuple[float, float]:
    values = []
    for histogram in histograms:
        values.extend(value for value in histogram.band_low if math.isfinite(value))
        values.extend(value for value in histogram.band_high if math.isfinite(value))
        values.extend(value for value in histogram.values if math.isfinite(value))

    if not values:
        return 0.82, 1.18

    minimum = min(values)
    maximum = max(values)
    span = max(abs(maximum - 1.0), abs(1.0 - minimum))
    span = max(span * 1.12, 0.12)
    return max(0.0, 1.0 - span), 1.0 + span


def mostly_positive(*histograms: BinnedHistogram) -> bool:
    finite_values = [value for histogram in histograms for value in histogram.values if math.isfinite(value)]
    if not finite_values:
        return False
    positive_count = sum(1 for value in finite_values if value > 0.0)
    return positive_count / len(finite_values) >= 0.9


def configure_main_axis(axis, mg_histogram: BinnedHistogram, local_histogram: BinnedHistogram) -> tuple[bool, float]:
    floor = positive_floor(mg_histogram, local_histogram)
    positive_maxima = [
        value
        for histogram in (mg_histogram, local_histogram)
        for value in (histogram.values + histogram.band_high)
        if math.isfinite(value) and value > 0.0
    ]
    top = max(positive_maxima, default=1.0) * 1.18
    if histogram_has_positive_content(mg_histogram) and histogram_has_positive_content(local_histogram) and mostly_positive(mg_histogram, local_histogram):
        axis.set_yscale("log")
        axis.set_ylim(bottom=floor * 0.55, top=top)
        return True, floor
    else:
        axis.set_yscale("symlog", linthresh=floor * 0.5)
        axis.set_ylim(top=top)
        return False, floor


def style_axis(axis, *, ratio: bool, x_major_step: float) -> None:
    axis.set_facecolor("white")
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color(FRAME_COLOR)
    axis.spines["bottom"].set_color(FRAME_COLOR)
    axis.tick_params(colors=FRAME_COLOR, which="both")
    axis.grid(True, axis="y", which="major", color=GRID_COLOR, linewidth=0.8, alpha=0.75)
    axis.grid(True, axis="x", which="major", color=GRID_COLOR, linewidth=0.45, alpha=0.35)
    if ratio:
        axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
        axis.xaxis.set_major_locator(MultipleLocator(x_major_step))
        axis.xaxis.set_minor_locator(AutoMinorLocator(2))
    else:
        axis.xaxis.set_major_locator(MultipleLocator(x_major_step))
        axis.xaxis.set_minor_locator(AutoMinorLocator(2))


def plot_histogram_comparison(order: str,
                              observable: ObservableSpec,
                              mg_histogram: BinnedHistogram,
                              local_histogram: BinnedHistogram,
                              output_prefix: Path,
                              sqrt_s_tev: float,
                              pdf_set: str,
                              mt_gev: float,
                              mu_gev: float) -> None:
    overlap_low, overlap_high = overlapping_range(mg_histogram, local_histogram)
    mg_trimmed = trim_histogram(mg_histogram, overlap_low, overlap_high, mg_histogram.name)
    local_trimmed = trim_histogram(local_histogram, overlap_low, overlap_high, local_histogram.name)

    mg_scale_ratio = mask_unstable_ratio_bins(mg_trimmed, mg_ratio_band(mg_trimmed, "scale"))
    mg_stat_ratio = mask_unstable_ratio_bins(mg_trimmed, mg_ratio_band(mg_trimmed, "stat"))
    local_ratio = mask_unstable_ratio_bins(mg_trimmed, ratio_histogram_with_band(local_trimmed, mg_trimmed, f"{local_trimmed.name}_ratio"))

    figure, (axis, ratio_axis) = plt.subplots(
        2,
        1,
        figsize=(9, 7),
        sharex=True,
        facecolor="white",
        gridspec_kw={"height_ratios": (3, 1), "hspace": 0.05},
    )
    style_axis(axis, ratio=False, x_major_step=observable.x_major_step)
    style_axis(ratio_axis, ratio=True, x_major_step=observable.x_major_step)

    mg_edges = mg_trimmed.edges
    local_edges = local_trimmed.edges

    use_log_scale, floor = configure_main_axis(axis, mg_trimmed, local_trimmed)

    if use_log_scale:
        mg_values_step = positive_step_values(mg_trimmed.values)
        mg_low_step = clipped_positive_step_values(mg_trimmed.band_low, floor * 0.55)
        mg_high_step = clipped_positive_step_values(mg_trimmed.band_high, floor * 0.55)
        local_values_step = positive_step_values(local_trimmed.values)
        local_low_step = clipped_positive_step_values(local_trimmed.band_low, floor * 0.55)
        local_high_step = clipped_positive_step_values(local_trimmed.band_high, floor * 0.55)
    else:
        mg_values_step = step_values(mg_trimmed.values)
        mg_low_step = step_values(mg_trimmed.band_low)
        mg_high_step = step_values(mg_trimmed.band_high)
        local_values_step = step_values(local_trimmed.values)
        local_low_step = step_values(local_trimmed.band_low)
        local_high_step = step_values(local_trimmed.band_high)

    axis.fill_between(mg_edges, mg_low_step, mg_high_step, step="post", color=MG_BAND_COLOR, alpha=0.35, linewidth=0.0, zorder=1)
    axis.plot(mg_edges, mg_values_step, color=MG_LINE_COLOR, linewidth=2.2, drawstyle="steps-post", zorder=3)
    axis.fill_between(local_edges, local_low_step, local_high_step, step="post", color=LOCAL_BAND_COLOR, alpha=0.28, linewidth=0.0, zorder=2)
    axis.plot(local_edges, local_values_step, color=LOCAL_LINE_COLOR, linewidth=2.2, drawstyle="steps-post", zorder=4)

    axis.set_ylabel("sigma_bin [pb]")
    axis.set_xlim(overlap_low, overlap_high)
    axis.tick_params(labelbottom=False)

    axis.text(
        0.012,
        0.985,
        "(a)",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
        color=FRAME_COLOR,
    )

    main_handles = [
        Patch(facecolor=MG_BAND_COLOR, edgecolor="none", alpha=0.35, label="MadGraph scale band"),
        Line2D([0], [0], color=MG_LINE_COLOR, linewidth=2.2, label="MadGraph central"),
        Patch(facecolor=LOCAL_BAND_COLOR, edgecolor="none", alpha=0.28, label=local_trimmed.band_label),
        Line2D([0], [0], color=LOCAL_LINE_COLOR, linewidth=2.2, label=f"Local {order}"),
    ]
    axis.legend(
        handles=main_handles,
        loc="upper right",
        frameon=True,
        framealpha=0.96,
        fancybox=False,
        edgecolor=GRID_COLOR,
        borderpad=0.55,
        handlelength=2.1,
        fontsize=9,
    )

    mg_scale_low_step = step_values(mg_scale_ratio.band_low)
    mg_scale_high_step = step_values(mg_scale_ratio.band_high)
    mg_stat_low_step = step_values(mg_stat_ratio.band_low)
    mg_stat_high_step = step_values(mg_stat_ratio.band_high)
    local_ratio_step = step_values(local_ratio.values)
    local_ratio_low_step = step_values(local_ratio.band_low)
    local_ratio_high_step = step_values(local_ratio.band_high)

    ratio_axis.axhspan(0.95, 1.05, color="#F5F0E7", alpha=0.9, zorder=0)
    ratio_axis.fill_between(mg_edges, mg_scale_low_step, mg_scale_high_step, step="post", color=MG_BAND_COLOR, alpha=0.35, linewidth=0.0)
    ratio_axis.fill_between(mg_edges, mg_stat_low_step, mg_stat_high_step, step="post", color=MG_STAT_COLOR, alpha=0.95, linewidth=0.0)
    ratio_axis.fill_between(local_edges, local_ratio_low_step, local_ratio_high_step, step="post", color=LOCAL_BAND_COLOR, alpha=0.25, linewidth=0.0)
    ratio_axis.plot(local_edges, local_ratio_step, color=LOCAL_LINE_COLOR, linewidth=2.0, drawstyle="steps-post")
    ratio_axis.axhline(1.0, color=FRAME_COLOR, linestyle=(0, (4, 2)), linewidth=1.2, alpha=0.9)
    ratio_axis.set_ylabel("Local / MG")
    ratio_axis.set_xlabel(observable.axis_label)
    ratio_axis.set_ylim(*finite_ratio_bounds(mg_scale_ratio, mg_stat_ratio, local_ratio))
    ratio_axis.text(
        0.012,
        0.96,
        "(b)",
        transform=ratio_axis.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
        color=FRAME_COLOR,
    )
    ratio_handles = [
        Patch(facecolor=MG_BAND_COLOR, edgecolor="none", alpha=0.35, label="MG scale"),
        Patch(facecolor=MG_STAT_COLOR, edgecolor="none", alpha=0.95, label="MG stat."),
        Patch(facecolor=LOCAL_BAND_COLOR, edgecolor="none", alpha=0.25, label="Local scale" if local_trimmed.band_label == "Local scale band" else "Local stat."),
        Line2D([0], [0], color=LOCAL_LINE_COLOR, linewidth=2.0, label="Local / MG"),
    ]
    ratio_axis.legend(
        handles=ratio_handles,
        loc="best",
        ncol=2,
        frameon=True,
        framealpha=0.96,
        fancybox=False,
        edgecolor=GRID_COLOR,
        borderpad=0.45,
        handlelength=1.8,
        columnspacing=1.3,
        fontsize=8,
    )

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(output_prefix.with_suffix(".png"), bbox_inches="tight", dpi=300)
    figure.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    plt.close(figure)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot local vs MadGraph observable comparisons from cs_ppttb outputs")
    parser.add_argument("--observable-key", default="mtt", help="Observable key: mtt, top_pt, tb_pt, or ytt")
    parser.add_argument("--observable", default=None, help="Exact observable title inside the .hwu files")
    parser.add_argument("--mg-lo-hwu", default="data/references/mg_lo_ppttx_mtt_8tev_mu173.hwu")
    parser.add_argument("--mg-nlo-hwu", default="data/references/mg_nlo_ppttx_mtt_8tev_mu173.hwu")
    parser.add_argument("--local-lo-csv", default=None)
    parser.add_argument("--local-nlo-csv", default=None)
    parser.add_argument("--local-lo-envelope-csv", default=None)
    parser.add_argument("--local-nlo-envelope-csv", default=None)
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--sqrt-s-tev", type=float, default=8.0)
    parser.add_argument("--pdf-set", default="NNPDF31_nlo_as_0118")
    parser.add_argument("--mt-gev", type=float, default=173.0)
    parser.add_argument("--mu-gev", type=float, default=173.0)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    spec = observable_spec(arguments.observable_key)

    script_dir = Path(__file__).resolve().parent
    package_root = script_dir.parent
    result_root = Path(f"results/{spec.stem}_madgraph_reference_8tev")
    outdir = package_root / (arguments.outdir or result_root / "plots/python")

    mg_lo_hwu = package_root / arguments.mg_lo_hwu
    mg_nlo_hwu = package_root / arguments.mg_nlo_hwu
    local_lo_csv = package_root / (arguments.local_lo_csv or result_root / f"csv/{spec.stem}_lo.csv")
    local_nlo_csv = package_root / (arguments.local_nlo_csv or result_root / f"csv/{spec.stem}_nlo.csv")
    local_lo_envelope_csv = package_root / (arguments.local_lo_envelope_csv or result_root / f"csv/{spec.stem}_lo_scale_envelope.csv")
    local_nlo_envelope_csv = package_root / (arguments.local_nlo_envelope_csv or result_root / f"csv/{spec.stem}_nlo_scale_envelope.csv")
    mg_observable = arguments.observable or spec.mg_title

    if not mg_lo_hwu.exists():
        raise RuntimeError(f"MadGraph LO reference does not exist: {mg_lo_hwu}")
    local_lo_source = local_lo_envelope_csv if local_lo_envelope_csv.exists() else local_lo_csv
    if not local_lo_source.exists():
        raise RuntimeError(f"Local LO histogram does not exist: {local_lo_source}")

    mg_lo = parse_hwu_histogram(mg_lo_hwu, mg_observable)
    local_lo = parse_local_csv(local_lo_source, "local_lo")
    plot_histogram_comparison("LO", spec, mg_lo, local_lo, outdir / f"{spec.stem}_lo_vs_madgraph", arguments.sqrt_s_tev, arguments.pdf_set, arguments.mt_gev, arguments.mu_gev)
    print(f"Wrote LO comparison plots to {outdir}")

    local_nlo_source = local_nlo_envelope_csv if local_nlo_envelope_csv.exists() else local_nlo_csv
    if mg_nlo_hwu.exists() and local_nlo_source.exists():
        mg_nlo = parse_hwu_histogram(mg_nlo_hwu, mg_observable)
        local_nlo = parse_local_csv(local_nlo_source, "local_nlo")
        plot_histogram_comparison("NLO", spec, mg_nlo, local_nlo, outdir / f"{spec.stem}_nlo_vs_madgraph", arguments.sqrt_s_tev, arguments.pdf_set, arguments.mt_gev, arguments.mu_gev)
        print(f"Wrote NLO comparison plots to {outdir}")
    else:
        print(
            f"Skipping NLO comparison because a local {spec.stem}_nlo.csv is not available yet. "
            "The script is ready to plot it as soon as the packaged RA differential input exists.",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())