#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import subprocess
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class VariationSpec:
    label: str
    factor: float


@dataclass(frozen=True)
class ObservableSpec:
    key: str
    config_observable: str
    stem: str


OBSERVABLES = {
    "mtt": ObservableSpec("mtt", "mtt", "mtt"),
    "top_pt": ObservableSpec("top_pt", "top_pt", "top_pt"),
    "tb_pt": ObservableSpec("tb_pt", "tb_pt", "tb_pt"),
    "ytt": ObservableSpec("ytt", "y_tt", "ytt"),
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


VARIATIONS = (
    VariationSpec("0p5mt", 0.5),
    VariationSpec("1p0mt", 1.0),
    VariationSpec("2p0mt", 2.0),
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate local observable scale-variation data for MadGraph comparison plots")
    parser.add_argument("--observable-key", default="mtt")
    parser.add_argument("--base-config", default=None)
    parser.add_argument("--mt-gev", type=float, default=173.0)
    parser.add_argument("--ra-samples", type=int, default=100000)
    parser.add_argument("--lo-samples", type=int, default=0)
    parser.add_argument("--ca-samples", type=int, default=0)
    parser.add_argument("--va-samples", type=int, default=0)
    parser.add_argument("--plot-after", action="store_true", help="Rerun the Python comparison plot after generating the scale-variation data")
    return parser.parse_args()


def load_config_entries(path: Path) -> OrderedDict[str, str]:
    entries: OrderedDict[str, str] = OrderedDict()
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        key, value = line.split("=", 1)
        entries[key.strip()] = value.strip()
    return entries


def write_config(path: Path, entries: OrderedDict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "\n".join(f"{key} = {value}" for key, value in entries.items()) + "\n"
    path.write_text(text)


def configure_runtime_environment() -> dict[str, str]:
    env = os.environ.copy()
    conda_prefix = env.get("CONDA_PREFIX", "/opt/homebrew/Caskroom/miniconda/base/envs/lhapdf")
    lhapdf_data = Path(conda_prefix) / "share" / "LHAPDF"
    lhapdf_lib = Path(conda_prefix) / "lib"
    if "LHAPDF_DATA_PATH" not in env and lhapdf_data.exists():
        env["LHAPDF_DATA_PATH"] = str(lhapdf_data)
    if lhapdf_lib.exists():
        existing = env.get("DYLD_LIBRARY_PATH", "")
        env["DYLD_LIBRARY_PATH"] = f"{lhapdf_lib}:{existing}" if existing else str(lhapdf_lib)
    return env


def run_command(command: list[str], *, cwd: Path, env: dict[str, str]) -> None:
    print("Running:", " ".join(command))
    subprocess.run(command, cwd=cwd, env=env, check=True)


def read_histogram(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "bin_low_gev": float(row["bin_low_gev"]),
                    "bin_high_gev": float(row["bin_high_gev"]),
                    "value": float(row["value"]),
                    "error": float(row["error"]),
                }
            )
    if not rows:
        raise RuntimeError(f"Empty histogram CSV: {path}")
    return rows


def write_scale_envelope(path: Path,
                         central: list[dict[str, float]],
                         variations: list[list[dict[str, float]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["bin_low_gev", "bin_high_gev", "value", "error", "scale_low", "scale_high"],
        )
        writer.writeheader()
        for index, central_row in enumerate(central):
            values = [histogram[index]["value"] for histogram in variations]
            writer.writerow(
                {
                    "bin_low_gev": f"{central_row['bin_low_gev']:.17e}",
                    "bin_high_gev": f"{central_row['bin_high_gev']:.17e}",
                    "value": f"{central_row['value']:.17e}",
                    "error": f"{central_row['error']:.17e}",
                    "scale_low": f"{min(values):.17e}",
                    "scale_high": f"{max(values):.17e}",
                }
            )


def format_scale_label(spec: VariationSpec) -> str:
    return spec.label


def main() -> int:
    arguments = parse_arguments()
    spec = observable_spec(arguments.observable_key)

    script_dir = Path(__file__).resolve().parent
    package_root = script_dir.parent
    repo_root = package_root.parent
    base_config_path = package_root / (arguments.base_config or f"config/{spec.stem}_madgraph_reference_8tev.cfg")
    if not base_config_path.exists():
        raise RuntimeError(f"Base config does not exist: {base_config_path}")

    cs_binary = package_root / "build" / "bin" / "cs_ppttb"
    ra_binary = repo_root / "build" / "run_ra_integral_independent"
    if not cs_binary.exists():
        raise RuntimeError(f"Package binary does not exist: {cs_binary}")
    if not ra_binary.exists():
        raise RuntimeError(f"RA binary does not exist: {ra_binary}")

    env = configure_runtime_environment()
    entries = load_config_entries(base_config_path)
    sqrt_s_gev = float(entries.get("sqrt_s_gev", "8000.0"))
    output_dir_relative = Path(entries["output_dir"])
    central_output_dir = package_root / output_dir_relative
    scale_root = central_output_dir / "scale_variations"
    generated_config_dir = package_root / "config" / "generated" / "scale_variations"
    validation_root = package_root / "validation" / "scale_variations"
    validation_root.mkdir(parents=True, exist_ok=True)

    bin_count = int(entries["mtt_bins"])
    observable_min = float(entries["mtt_min_gev"])
    observable_max = float(entries["mtt_max_gev"])

    lo_histograms: dict[str, list[dict[str, float]]] = {}
    nlo_histograms: dict[str, list[dict[str, float]]] = {}

    for variation in VARIATIONS:
        mu = variation.factor * arguments.mt_gev
        label = format_scale_label(variation)
        if variation.factor == 1.0:
            output_dir = central_output_dir
        else:
            output_dir = scale_root / label

        ra_artifact = validation_root / f"local_ra_integral_{spec.stem}_mt173_mu_{label}_{arguments.ra_samples // 1000}k.txt"
        ra_convergence = validation_root / f"local_ra_integral_{spec.stem}_mt173_mu_{label}_{arguments.ra_samples // 1000}k_convergence.txt"
        ra_histogram = validation_root / f"local_ra_{spec.stem}_8tev_mu_{label}_{arguments.ra_samples // 1000}k.csv"

        run_command(
            [
                str(ra_binary),
                "--samples",
                str(arguments.ra_samples),
                "--output",
                str(ra_artifact),
                "--convergence-output",
                str(ra_convergence),
                "--histogram-output",
                str(ra_histogram),
                "--observable",
                spec.config_observable,
                "--sqrt-s-gev",
                str(sqrt_s_gev),
                "--mtt-bins",
                str(bin_count),
                "--mtt-min",
                str(observable_min),
                "--mtt-max",
                str(observable_max),
                "--mt-gev",
                str(arguments.mt_gev),
                "--muF-gev",
                str(mu),
                "--muR-gev",
                str(mu),
            ],
            cwd=repo_root,
            env=env,
        )

        run_entries = OrderedDict(entries)
        run_entries["observable"] = spec.config_observable
        run_entries["mt_gev"] = f"{arguments.mt_gev:.1f}"
        run_entries["muF_gev"] = f"{mu:.1f}"
        run_entries["muR_gev"] = f"{mu:.1f}"
        run_entries["ra_artifact"] = os.path.relpath(ra_artifact, package_root)
        run_entries["ra_histogram_csv"] = os.path.relpath(ra_histogram, package_root)
        run_entries["output_dir"] = os.path.relpath(output_dir, package_root)
        run_entries["write_plots"] = "false"
        if arguments.lo_samples > 0:
            run_entries["lo_samples"] = str(arguments.lo_samples)
        if arguments.ca_samples > 0:
            run_entries["ca_samples"] = str(arguments.ca_samples)
        if arguments.va_samples > 0:
            run_entries["va_samples"] = str(arguments.va_samples)

        generated_config = generated_config_dir / f"{spec.stem}_madgraph_reference_8tev_mu_{label}.cfg"
        write_config(generated_config, run_entries)
        run_command([str(cs_binary), "--config", str(generated_config)], cwd=package_root, env=env)

        csv_dir = output_dir / "csv"
        lo_histograms[label] = read_histogram(csv_dir / f"{spec.stem}_lo.csv")
        nlo_histograms[label] = read_histogram(csv_dir / f"{spec.stem}_nlo.csv")

    central_lo = lo_histograms["1p0mt"]
    central_nlo = nlo_histograms["1p0mt"]
    lo_variations = [lo_histograms[variation.label] for variation in VARIATIONS]
    nlo_variations = [nlo_histograms[variation.label] for variation in VARIATIONS]

    central_csv_dir = central_output_dir / "csv"
    write_scale_envelope(central_csv_dir / f"{spec.stem}_lo_scale_envelope.csv", central_lo, lo_variations)
    write_scale_envelope(central_csv_dir / f"{spec.stem}_nlo_scale_envelope.csv", central_nlo, nlo_variations)
    print(f"Wrote local scale envelopes to {central_csv_dir}")

    if arguments.plot_after:
        plot_script = script_dir / "plot_mtt_madgraph_comparison.py"
        run_command(
            [
                os.fspath(Path(env.get("CONDA_PREFIX", "/opt/homebrew/Caskroom/miniconda/base/envs/lhapdf")) / "bin" / "python"),
                str(plot_script),
                "--observable-key",
                spec.key,
            ],
            cwd=package_root,
            env=env,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())