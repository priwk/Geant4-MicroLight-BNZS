#!/usr/bin/env python3
"""Generate Stage B effective ZnS(Ag) scintillation-source figures."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

DEPENDENCY_ERROR = """
Failed to import the plotting dependencies.

Use a clean virtual environment from the project root:

    cd ~/g4work/B2
    python3 -m venv .venv
    . .venv/bin/activate
    python3 -m pip install --upgrade pip
    python3 -m pip install -r plots/stageB_scintillation_strength/requirements.txt
    python3 plots/stageB_scintillation_strength/generate_stageB_scintillation_strength.py
""".strip()

try:
    import matplotlib

    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import font_manager
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.lines import Line2D
except ImportError as exc:
    raise SystemExit(f"{DEPENDENCY_ERROR}\n\nOriginal error: {exc}") from None
except AttributeError as exc:
    if "_ARRAY_API" in str(exc):
        raise SystemExit(f"{DEPENDENCY_ERROR}\n\nOriginal error: {exc}") from None
    raise


SCRIPT_DIR = Path(__file__).resolve().parent
FONT_DIR = SCRIPT_DIR / "fonts"

DESIRED_RATIO_ORDER = ["2-1", "1-1", "1-1.5", "1-2", "1-2.5", "1-3"]
DISCRETE_RATIO_COLORS = [
    "#0b3c79",
    "#c0392b",
    "#188977",
    "#8e44ad",
    "#d68910",
    "#2d3436",
]
RAINBOW_CMAP = LinearSegmentedColormap.from_list(
    "red_orange_yellow_green_cyan_blue",
    ["#d73027", "#fc8d59", "#fee08b", "#91cf60", "#41b6c4", "#225ea8"],
)

RATIO_RE = re.compile(
    r"^(?P<bn>[+]?(?:\d+(?:\.\d*)?|\.\d+))-(?P<zns>[+]?(?:\d+(?:\.\d*)?|\.\d+))$"
)
ANCHOR_RE = re.compile(
    r"^(?P<thickness>[+]?(?:\d+(?:\.\d*)?|\.\d+))_capture_anchors\.csv$"
)
STAGEA_CAPTURE_RE = re.compile(
    r"^(?P<thickness>[+]?(?:\d+(?:\.\d*)?|\.\d+))_neutron_capture_positions\.csv$"
)


@dataclass(frozen=True)
class RatioKey:
    tag: str
    bn_wt: float
    zns_wt: float

    @property
    def display_tag(self) -> str:
        return self.tag.replace("-", ":")


def configure_fonts() -> None:
    font_files = [
        FONT_DIR / "Arial.ttf",
        FONT_DIR / "Arial-Bold.ttf",
        FONT_DIR / "Arial-Italic.ttf",
        FONT_DIR / "Arial-BoldItalic.ttf",
    ]
    for path in font_files:
        if path.is_file():
            font_manager.fontManager.addfont(str(path))

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "legend.fontsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.unicode_minus": False,
        }
    )


configure_fonts()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Stage B ZnS(Ag) effective energy-deposition metrics."
    )
    parser.add_argument(
        "--project-root",
        default=None,
        help="Project root. Defaults to the repository containing this script.",
    )
    parser.add_argument(
        "--stageb-root",
        default=None,
        help="Stage B output root. Defaults to <project-root>/Output/stageB.",
    )
    parser.add_argument(
        "--stagea-root",
        default=None,
        help="Stage A input root. Defaults to <project-root>/Input/stageA.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory. Defaults to this plotting folder's output/.",
    )
    parser.add_argument(
        "--birks-kb",
        type=float,
        default=0.0,
        help=(
            "Birks kB in um/keV when dE/dx is computed in keV/um. "
            "Use 0 to disable quenching. Default: 0."
        ),
    )
    parser.add_argument(
        "--effective-threshold-kev",
        type=float,
        default=0.0,
        help="Minimum effective ZnS(Ag) edep to count as a nonzero event.",
    )
    return parser.parse_args()


def parse_ratio_tag(tag: str) -> RatioKey:
    match = RATIO_RE.fullmatch(tag)
    if match is None:
        raise ValueError(f"Unsupported ratio folder name: {tag}")
    return RatioKey(
        tag=tag,
        bn_wt=float(match.group("bn")),
        zns_wt=float(match.group("zns")),
    )


def ratio_display_sort_key(ratio: RatioKey) -> tuple[int, float, float]:
    try:
        return (DESIRED_RATIO_ORDER.index(ratio.tag), ratio.bn_wt, ratio.zns_wt)
    except ValueError:
        return (len(DESIRED_RATIO_ORDER), ratio.bn_wt, ratio.zns_wt)


def thickness_sort_key(value: float) -> tuple[int, float]:
    rounded = round(value)
    is_int_like = math.isclose(value, rounded, rel_tol=0.0, abs_tol=1.0e-9)
    return (0 if is_int_like else 1, value)


def parse_anchor_thickness(path: Path) -> float:
    match = ANCHOR_RE.fullmatch(path.name)
    if match is None:
        raise ValueError(f"Unsupported Stage B anchor CSV name: {path.name}")
    return float(match.group("thickness"))


def parse_stagea_capture_thickness(path: Path) -> float:
    match = STAGEA_CAPTURE_RE.fullmatch(path.name)
    if match is None:
        raise ValueError(f"Unsupported Stage A capture CSV name: {path.name}")
    return float(match.group("thickness"))


def mu_m_text() -> str:
    return "μm"


def style_axes(ax) -> None:
    ax.tick_params(direction="in", which="both", top=True, right=True, width=1.3)
    for spine in ax.spines.values():
        spine.set_linewidth(1.6)


def safe_float(value: str, *, default: float | None = None) -> float:
    if value is None or value == "":
        if default is not None:
            return default
        raise ValueError("empty numeric value")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"non-finite numeric value: {value!r}")
    return number


def fmt(value: object) -> object:
    if isinstance(value, float):
        return f"{value:.12g}"
    return value


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key, "")) for key in fieldnames})


def visible_edep(edep_kev: float, step_len_um: float, birks_kb: float) -> float:
    if edep_kev <= 0.0:
        return 0.0
    if birks_kb == 0.0:
        return edep_kev
    if step_len_um <= 0.0:
        return 0.0
    dedx = edep_kev / step_len_um
    return edep_kev / (1.0 + birks_kb * dedx)


def read_source_edep(
    track_csv: Path,
    birks_kb: float,
) -> tuple[dict[str, tuple[float, float, int]], int]:
    source_raw: defaultdict[str, float] = defaultdict(float)
    source_visible: defaultdict[str, float] = defaultdict(float)
    source_steps: defaultdict[str, int] = defaultdict(int)
    step_count = 0

    if not track_csv.is_file():
        return {}, 0

    with track_csv.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"source_event_uid", "edep_keV", "step_len_um"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"{track_csv} is missing columns: {', '.join(sorted(missing))}")

        for row in reader:
            source_uid = row["source_event_uid"]
            edep_kev = safe_float(row["edep_keV"])
            step_len_um = safe_float(row["step_len_um"], default=0.0)
            source_raw[source_uid] += edep_kev
            source_visible[source_uid] += visible_edep(edep_kev, step_len_um, birks_kb)
            source_steps[source_uid] += 1
            step_count += 1

    return {
        source_uid: (
            source_raw[source_uid],
            source_visible[source_uid],
            source_steps[source_uid],
        )
        for source_uid in source_raw
    }, step_count


def summarize_stageb_pair(
    ratio: RatioKey,
    thickness_um: float,
    anchor_csv: Path,
    track_csv: Path,
    birks_kb: float,
    threshold_kev: float,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    source_edep, zns_step_count = read_source_edep(track_csv, birks_kb)
    physical_visible: defaultdict[str, float] = defaultdict(float)
    physical_raw: defaultdict[str, float] = defaultdict(float)
    physical_weight: defaultdict[str, float] = defaultdict(float)
    physical_source_count: defaultdict[str, int] = defaultdict(int)

    placement_total_weight: defaultdict[str, float] = defaultdict(float)
    placement_effective_weight: defaultdict[str, float] = defaultdict(float)
    placement_raw_sum: defaultdict[str, float] = defaultdict(float)
    placement_visible_sum: defaultdict[str, float] = defaultdict(float)
    placement_anchor_rows: defaultdict[str, int] = defaultdict(int)

    anchor_rows = 0
    total_weight = 0.0
    source_effective_weight = 0.0
    weighted_raw_sum = 0.0
    weighted_visible_sum = 0.0

    with anchor_csv.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {
            "physical_event_uid",
            "source_event_uid",
            "trajectory_weight",
            "placement_file",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"{anchor_csv} is missing columns: {', '.join(sorted(missing))}")

        for row in reader:
            anchor_rows += 1
            physical_uid = row["physical_event_uid"]
            source_uid = row["source_event_uid"]
            placement_file = row["placement_file"]
            weight = safe_float(row["trajectory_weight"], default=1.0)
            raw_kev, visible_kev, _step_count = source_edep.get(source_uid, (0.0, 0.0, 0))
            effective = visible_kev > threshold_kev

            physical_raw[physical_uid] += weight * raw_kev
            physical_visible[physical_uid] += weight * visible_kev
            physical_weight[physical_uid] += weight
            physical_source_count[physical_uid] += 1

            placement_total_weight[placement_file] += weight
            placement_raw_sum[placement_file] += weight * raw_kev
            placement_visible_sum[placement_file] += weight * visible_kev
            placement_anchor_rows[placement_file] += 1
            if effective:
                placement_effective_weight[placement_file] += weight
                source_effective_weight += weight

            total_weight += weight
            weighted_raw_sum += weight * raw_kev
            weighted_visible_sum += weight * visible_kev

    if anchor_rows == 0:
        raise SystemExit(f"No anchor rows found in {anchor_csv}")

    physical_expected_visible = []
    physical_expected_raw = []
    for physical_uid, weight_sum in physical_weight.items():
        if weight_sum <= 0.0:
            continue
        physical_expected_visible.append(physical_visible[physical_uid] / weight_sum)
        physical_expected_raw.append(physical_raw[physical_uid] / weight_sum)

    if not physical_expected_visible:
        raise SystemExit(f"No weighted physical captures found in {anchor_csv}")

    physical_visible_values = np.asarray(physical_expected_visible, dtype=float)
    physical_raw_values = np.asarray(physical_expected_raw, dtype=float)
    n_physical = int(physical_visible_values.size)
    physical_effective_count = int(np.count_nonzero(physical_visible_values > threshold_kev))

    placement_rows: list[dict[str, object]] = []
    placement_mean_values = []
    placement_fraction_values = []
    for placement_file in sorted(placement_total_weight):
        denominator = placement_total_weight[placement_file]
        if denominator <= 0.0:
            continue
        mean_visible = placement_visible_sum[placement_file] / denominator
        effective_fraction = placement_effective_weight[placement_file] / denominator
        placement_mean_values.append(mean_visible)
        placement_fraction_values.append(effective_fraction)
        placement_rows.append(
            {
                "ratio_tag": ratio.tag,
                "bn_wt": ratio.bn_wt,
                "zns_wt": ratio.zns_wt,
                "thickness_um": thickness_um,
                "placement_file": placement_file,
                "n_anchor_rows": placement_anchor_rows[placement_file],
                "capture_weight_sum": denominator,
                "effective_deposition_fraction_weighted": effective_fraction,
                "mean_raw_zns_edep_per_capture_keV": placement_raw_sum[placement_file] / denominator,
                "mean_effective_zns_edep_per_capture_keV": mean_visible,
            }
        )

    placement_mean_arr = np.asarray(placement_mean_values, dtype=float)
    placement_fraction_arr = np.asarray(placement_fraction_values, dtype=float)
    placement_mean_std = (
        float(np.std(placement_mean_arr, ddof=1)) if placement_mean_arr.size > 1 else 0.0
    )
    placement_fraction_std = (
        float(np.std(placement_fraction_arr, ddof=1)) if placement_fraction_arr.size > 1 else 0.0
    )

    record = {
        "ratio_tag": ratio.tag,
        "bn_wt": ratio.bn_wt,
        "zns_wt": ratio.zns_wt,
        "thickness_um": thickness_um,
        "anchor_rows": anchor_rows,
        "zns_step_rows": zns_step_count,
        "n_physical_captures": n_physical,
        "capture_weight_sum": total_weight,
        "n_effective_physical_captures": physical_effective_count,
        "effective_deposition_event_fraction": physical_effective_count / n_physical,
        "effective_deposition_fraction_weighted": source_effective_weight / total_weight,
        "mean_raw_zns_edep_per_capture_keV": float(np.mean(physical_raw_values)),
        "mean_effective_zns_edep_per_capture_keV": float(np.mean(physical_visible_values)),
        "median_effective_zns_edep_per_capture_keV": float(np.median(physical_visible_values)),
        "q25_effective_zns_edep_per_capture_keV": float(np.quantile(physical_visible_values, 0.25)),
        "q75_effective_zns_edep_per_capture_keV": float(np.quantile(physical_visible_values, 0.75)),
        "weighted_raw_zns_edep_sum_keV": weighted_raw_sum,
        "weighted_effective_zns_edep_sum_keV": weighted_visible_sum,
        "n_placements": len(placement_mean_values),
        "placement_mean_effective_edep_std_keV": placement_mean_std,
        "placement_effective_fraction_std": placement_fraction_std,
    }
    return record, placement_rows


def collect_stageb_metrics(
    stageb_root: Path,
    birks_kb: float,
    threshold_kev: float,
) -> tuple[list[RatioKey], list[float], list[dict[str, object]], list[dict[str, object]]]:
    records: list[dict[str, object]] = []
    placement_records: list[dict[str, object]] = []
    ratio_map: dict[str, RatioKey] = {}
    thickness_values: set[float] = set()

    if not stageb_root.is_dir():
        raise SystemExit(f"Stage B root not found: {stageb_root}")

    for ratio_dir in sorted(stageb_root.iterdir()):
        if not ratio_dir.is_dir():
            continue
        try:
            ratio = parse_ratio_tag(ratio_dir.name)
        except ValueError:
            continue

        for anchor_csv in sorted(ratio_dir.glob("*_capture_anchors.csv")):
            thickness_um = parse_anchor_thickness(anchor_csv)
            track_csv = ratio_dir / f"{thickness_um:g}_zns_track_steps.csv"
            if not track_csv.is_file():
                print(f"Skipping {anchor_csv}: missing {track_csv.name}")
                continue

            record, placement_rows = summarize_stageb_pair(
                ratio,
                thickness_um,
                anchor_csv,
                track_csv,
                birks_kb,
                threshold_kev,
            )
            records.append(record)
            placement_records.extend(placement_rows)
            ratio_map[ratio.tag] = ratio
            thickness_values.add(thickness_um)
            print(
                "Processed "
                f"{ratio.tag} t={thickness_um:g} um: "
                f"{record['n_physical_captures']} captures, "
                f"{record['effective_deposition_fraction_weighted']:.4g} weighted hit fraction"
            )

    if not records:
        raise SystemExit(f"No Stage B slim file pairs found in {stageb_root}")

    ratios = sorted(ratio_map.values(), key=ratio_display_sort_key)
    thicknesses = sorted(thickness_values, key=thickness_sort_key)
    records.sort(
        key=lambda row: (
            ratio_display_sort_key(parse_ratio_tag(str(row["ratio_tag"]))),
            thickness_sort_key(float(row["thickness_um"])),
        )
    )
    return ratios, thicknesses, records, placement_records


def read_stagea_absorption(stagea_root: Path) -> dict[tuple[str, float], float]:
    absorption: dict[tuple[str, float], float] = {}
    if not stagea_root.is_dir():
        return absorption

    for ratio_dir in sorted(stagea_root.iterdir()):
        if not ratio_dir.is_dir():
            continue
        ratio_tag = ratio_dir.name

        absorption_csv = ratio_dir / "neutron_capture_absorption" / "neutron_capture_absorption.csv"
        if absorption_csv.is_file():
            with absorption_csv.open(newline="", encoding="utf-8-sig") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    if "thickness_um" not in row:
                        continue
                    thickness_um = safe_float(row["thickness_um"])
                    if "absorb_efficiency" in row and row["absorb_efficiency"] != "":
                        absorption[(ratio_tag, thickness_um)] = safe_float(row["absorb_efficiency"])
                    elif "n_absorb" in row and "n_incident" in row:
                        n_incident = safe_float(row["n_incident"])
                        if n_incident > 0.0:
                            absorption[(ratio_tag, thickness_um)] = safe_float(row["n_absorb"]) / n_incident

        capture_dir = ratio_dir / "neutron_capture_positions"
        if not capture_dir.is_dir():
            continue
        for capture_csv in sorted(capture_dir.glob("*_neutron_capture_positions.csv")):
            thickness_um = parse_stagea_capture_thickness(capture_csv)
            key = (ratio_tag, thickness_um)
            if key in absorption:
                continue

            count = 0
            max_event_id = -1
            with capture_csv.open(newline="", encoding="utf-8-sig") as handle:
                reader = csv.DictReader(handle)
                if "eventID" not in (reader.fieldnames or []):
                    continue
                for row in reader:
                    count += 1
                    max_event_id = max(max_event_id, int(float(row["eventID"])))
            if count > 0 and max_event_id >= 0:
                absorption[key] = count / (max_event_id + 1)

    return absorption


def add_initial_strength(
    records: list[dict[str, object]],
    stagea_absorption: dict[tuple[str, float], float],
) -> None:
    products: list[float] = []
    for row in records:
        ratio_tag = str(row["ratio_tag"])
        thickness_um = float(row["thickness_um"])
        absorption_rate = stagea_absorption.get((ratio_tag, thickness_um), math.nan)
        row["stageA_absorption_rate"] = absorption_rate
        if math.isfinite(absorption_rate):
            product = absorption_rate * float(row["mean_effective_zns_edep_per_capture_keV"])
        else:
            product = math.nan
        row["initial_scintillation_strength_absorption_x_edep"] = product
        if math.isfinite(product):
            products.append(product)

    max_product = max(products) if products else math.nan
    for row in records:
        product = float(row["initial_scintillation_strength_absorption_x_edep"])
        if math.isfinite(product) and math.isfinite(max_product) and max_product > 0.0:
            row["relative_initial_scintillation_strength"] = product / max_product
        else:
            row["relative_initial_scintillation_strength"] = math.nan


def records_to_matrix(
    ratios: list[RatioKey],
    thicknesses: list[float],
    records: list[dict[str, object]],
    field: str,
) -> np.ndarray:
    matrix = np.full((len(ratios), len(thicknesses)), np.nan, dtype=float)
    ratio_index = {ratio.tag: idx for idx, ratio in enumerate(ratios)}
    thickness_index = {value: idx for idx, value in enumerate(thicknesses)}
    for row in records:
        value = row.get(field, math.nan)
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(number):
            continue
        matrix[ratio_index[str(row["ratio_tag"])], thickness_index[float(row["thickness_um"])]] = number
    return matrix


def make_log_axis_edges(values: np.ndarray) -> np.ndarray:
    coords = np.asarray(values, dtype=float)
    if coords.size == 1:
        return np.asarray([coords[0] / 1.5, coords[0] * 1.5], dtype=float)
    logs = np.log10(coords)
    edge_logs = np.empty(coords.size + 1, dtype=float)
    edge_logs[1:-1] = 0.5 * (logs[:-1] + logs[1:])
    edge_logs[0] = logs[0] - 0.5 * (logs[1] - logs[0])
    edge_logs[-1] = logs[-1] + 0.5 * (logs[-1] - logs[-2])
    return np.power(10.0, edge_logs)


def interp_line(x: np.ndarray, y: np.ndarray, n: int = 800) -> tuple[np.ndarray, np.ndarray]:
    if x.size < 2:
        return x, y
    if np.all(x > 0.0):
        x_dense = np.geomspace(float(np.min(x)), float(np.max(x)), n)
        y_dense = np.interp(np.log10(x_dense), np.log10(x), y)
    else:
        x_dense = np.linspace(float(np.min(x)), float(np.max(x)), n)
        y_dense = np.interp(x_dense, x, y)
    return x_dense, y_dense


def plot_metric_lines(
    ratios: list[RatioKey],
    thicknesses: list[float],
    matrix: np.ndarray,
    output_path: Path,
    ylabel: str,
    title: str,
    *,
    y_min: float = 0.0,
    y_max: float | None = None,
    placement_std_matrix: np.ndarray | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 6.0))
    x_all = np.asarray(thicknesses, dtype=float)
    handles: list[Line2D] = []

    for idx, ratio in enumerate(ratios):
        row = matrix[idx, :]
        valid = np.isfinite(row)
        if np.count_nonzero(valid) == 0:
            continue
        x = x_all[valid]
        y = row[valid]
        color = DISCRETE_RATIO_COLORS[idx % len(DISCRETE_RATIO_COLORS)]
        x_dense, y_dense = interp_line(x, y)
        ax.plot(x_dense, y_dense, color=color, linewidth=2.0, alpha=0.86)
        ax.scatter(x, y, s=24, color=color, alpha=0.95, zorder=5)

        if placement_std_matrix is not None:
            std = placement_std_matrix[idx, :][valid]
            if np.isfinite(std).any() and x.size >= 2:
                _, std_dense = interp_line(x, std)
                ax.fill_between(
                    x_dense,
                    y_dense - std_dense,
                    y_dense + std_dense,
                    color=color,
                    alpha=0.13,
                    linewidth=0.0,
                )
        handles.append(Line2D([0], [0], color=color, linewidth=2.0, label=ratio.display_tag))

    if np.nanmin(x_all) > 0.0:
        ax.set_xscale("log")
        common_ticks = [1, 10, 30, 100, 300, 1000]
        ticks = [tick for tick in common_ticks if np.nanmin(x_all) <= tick <= np.nanmax(x_all)]
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(tick) for tick in ticks])
    ax.set_xlim(float(np.nanmin(x_all)), float(np.nanmax(x_all)))

    finite = matrix[np.isfinite(matrix)]
    if finite.size:
        if y_max is None:
            y_max = max(y_min + 1.0e-12, float(np.max(finite)) * 1.08)
        ax.set_ylim(y_min, y_max)
    ax.set_xlabel(f"Thickness t ({mu_m_text()})")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.18, linewidth=0.7)
    ax.legend(handles=handles, frameon=False, ncol=2)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def interpolate_rows_along_log_thickness(
    matrix: np.ndarray,
    thicknesses: list[float],
    sample_count: int = 800,
) -> tuple[np.ndarray, np.ndarray]:
    x_known = np.asarray(thicknesses, dtype=float)
    if x_known.size == 1:
        return x_known, matrix.copy()

    x_dense = np.geomspace(float(np.min(x_known)), float(np.max(x_known)), sample_count)
    dense = np.full((matrix.shape[0], sample_count), np.nan, dtype=float)
    log_x_known = np.log10(x_known)
    log_x_dense = np.log10(x_dense)

    for row_index in range(matrix.shape[0]):
        row = matrix[row_index, :]
        valid = np.isfinite(row)
        if np.count_nonzero(valid) == 1:
            dense[row_index, :] = row[valid][0]
        elif np.count_nonzero(valid) >= 2:
            dense[row_index, :] = np.interp(log_x_dense, log_x_known[valid], row[valid])
    return x_dense, dense


def plot_metric_heatmap(
    ratios: list[RatioKey],
    thicknesses: list[float],
    matrix: np.ndarray,
    output_path: Path,
    colorbar_label: str,
    title: str,
    *,
    vmin: float = 0.0,
    vmax: float | None = None,
) -> None:
    x_dense, matrix_dense = interpolate_rows_along_log_thickness(matrix, thicknesses)
    x_edges = make_log_axis_edges(x_dense)
    gap = 0.18
    bar_height = 0.82

    fig, ax = plt.subplots(figsize=(8.0, 6.4))
    cmap = RAINBOW_CMAP.copy()
    cmap.set_bad(color="#ececec")
    if vmax is None:
        finite = matrix_dense[np.isfinite(matrix_dense)]
        vmax = float(np.max(finite)) if finite.size else 1.0
    image = None
    y_centers: list[float] = []

    for idx, ratio in enumerate(ratios):
        y_bottom = idx * (bar_height + gap)
        y_top = y_bottom + bar_height
        y_centers.append(0.5 * (y_bottom + y_top))
        image = ax.pcolormesh(
            x_edges,
            np.asarray([y_bottom, y_top], dtype=float),
            matrix_dense[idx : idx + 1, :],
            shading="flat",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )

    ax.set_xscale("log")
    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(len(ratios) * (bar_height + gap) - gap, 0.0)
    common_ticks = [1, 10, 30, 100, 300, 1000]
    ticks = [tick for tick in common_ticks if x_edges[0] <= tick <= x_edges[-1]]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(tick) for tick in ticks])
    ax.set_yticks(y_centers)
    ax.set_yticklabels([ratio.display_tag for ratio in ratios])
    ax.set_xlabel(f"Thickness t ({mu_m_text()})")
    ax.set_ylabel("BN/ZnS(Ag) ratio")
    ax.set_title(title)
    style_axes(ax)

    colorbar = fig.colorbar(image, ax=ax, pad=0.06)
    colorbar.set_label(colorbar_label)
    colorbar.ax.tick_params(direction="in", which="both", width=1.3)
    colorbar.outline.set_linewidth(1.6)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def summarize_ratio_bars(
    ratios: list[RatioKey],
    records: list[dict[str, object]],
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    means = []
    stds = []
    counts = []
    for ratio in ratios:
        values = [
            float(row["mean_effective_zns_edep_per_capture_keV"])
            for row in records
            if row["ratio_tag"] == ratio.tag
            and math.isfinite(float(row["mean_effective_zns_edep_per_capture_keV"]))
        ]
        arr = np.asarray(values, dtype=float)
        counts.append(int(arr.size))
        if arr.size == 0:
            means.append(math.nan)
            stds.append(math.nan)
        elif arr.size == 1:
            means.append(float(arr[0]))
            stds.append(0.0)
        else:
            means.append(float(np.mean(arr)))
            stds.append(float(np.std(arr, ddof=1)))
    return np.asarray(means, dtype=float), np.asarray(stds, dtype=float), counts


def plot_ratio_bar_with_dispersion(
    ratios: list[RatioKey],
    records: list[dict[str, object]],
    placement_records: list[dict[str, object]],
    output_path: Path,
) -> None:
    means, stds, counts = summarize_ratio_bars(ratios, records)
    fig, ax = plt.subplots(figsize=(7.4, 5.6))
    x = np.arange(len(ratios), dtype=float)
    colors = [DISCRETE_RATIO_COLORS[idx % len(DISCRETE_RATIO_COLORS)] for idx in range(len(ratios))]
    ax.bar(x, means, yerr=stds, capsize=4, color=colors, alpha=0.82, edgecolor="#222222", linewidth=0.7)

    for idx, ratio in enumerate(ratios):
        values = [
            float(row["mean_effective_zns_edep_per_capture_keV"])
            for row in placement_records
            if row["ratio_tag"] == ratio.tag
        ]
        if not values:
            continue
        arr = np.asarray(values, dtype=float)
        if arr.size > 250:
            rng = np.random.default_rng(20260612 + idx)
            arr = rng.choice(arr, size=250, replace=False)
        jitter = np.linspace(-0.18, 0.18, arr.size) if arr.size > 1 else np.asarray([0.0])
        ax.scatter(
            np.full(arr.size, x[idx]) + jitter,
            arr,
            s=8,
            color="#1f1f1f",
            alpha=0.18,
            linewidth=0.0,
            zorder=5,
        )
        ax.text(
            x[idx],
            max(0.0, means[idx]) * 1.02 if math.isfinite(means[idx]) else 0.0,
            f"n={counts[idx]}",
            ha="center",
            va="bottom",
            fontsize=9,
            color="#333333",
        )

    finite = means[np.isfinite(means)]
    ax.set_ylim(0.0, max(1.0, float(np.max(finite)) * 1.18) if finite.size else 1.0)
    ax.set_xticks(x)
    ax.set_xticklabels([ratio.display_tag for ratio in ratios])
    ax.set_xlabel("BN/ZnS(Ag) ratio")
    ax.set_ylabel("Mean effective ZnS(Ag) edep per capture (keV)")
    ax.set_title("Mean effective ZnS(Ag) edep by ratio")
    ax.grid(True, axis="y", alpha=0.18, linewidth=0.7)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_combined_fig3(
    ratios: list[RatioKey],
    thicknesses: list[float],
    fraction_matrix: np.ndarray,
    edep_matrix: np.ndarray,
    strength_matrix: np.ndarray,
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), constrained_layout=True)
    x_all = np.asarray(thicknesses, dtype=float)
    specs = [
        (
            axes[0],
            fraction_matrix,
            "(a)",
            "Effective deposition fraction",
            "Weighted fraction",
            1.0,
        ),
        (
            axes[1],
            edep_matrix,
            "(b)",
            "Mean effective ZnS(Ag) edep",
            "Mean edep per capture (keV)",
            None,
        ),
        (
            axes[2],
            strength_matrix,
            "(c)",
            "Relative initial scintillation strength",
            "Relative strength",
            1.0,
        ),
    ]

    for ax, matrix, label, title, ylabel, ymax in specs:
        handles = []
        for idx, ratio in enumerate(ratios):
            row = matrix[idx, :]
            valid = np.isfinite(row)
            if np.count_nonzero(valid) == 0:
                continue
            x = x_all[valid]
            y = row[valid]
            color = DISCRETE_RATIO_COLORS[idx % len(DISCRETE_RATIO_COLORS)]
            x_dense, y_dense = interp_line(x, y, n=500)
            ax.plot(x_dense, y_dense, color=color, linewidth=2.0, alpha=0.86)
            ax.scatter(x, y, s=18, color=color, alpha=0.95, zorder=5)
            handles.append(Line2D([0], [0], color=color, linewidth=2.0, label=ratio.display_tag))

        if np.nanmin(x_all) > 0.0:
            ax.set_xscale("log")
            ticks = [tick for tick in [10, 30, 100, 300, 1000] if np.nanmin(x_all) <= tick <= np.nanmax(x_all)]
            ax.set_xticks(ticks)
            ax.set_xticklabels([str(tick) for tick in ticks])
        ax.set_xlim(float(np.nanmin(x_all)), float(np.nanmax(x_all)))
        finite = matrix[np.isfinite(matrix)]
        if finite.size:
            upper = ymax if ymax is not None else max(1.0e-12, float(np.max(finite)) * 1.12)
            ax.set_ylim(0.0, upper)
        ax.set_xlabel(f"Thickness t ({mu_m_text()})")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{label} {title}", fontsize=12)
        ax.grid(True, which="both", alpha=0.18, linewidth=0.7)
        style_axes(ax)
        if handles:
            ax.legend(handles=handles, frameon=False, fontsize=8.8, loc="best")

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    if args.birks_kb < 0.0:
        raise SystemExit("--birks-kb must be >= 0")
    if args.effective_threshold_kev < 0.0:
        raise SystemExit("--effective-threshold-kev must be >= 0")

    project_root = Path(args.project_root).resolve() if args.project_root else SCRIPT_DIR.parent.parent
    stageb_root = Path(args.stageb_root).resolve() if args.stageb_root else project_root / "Output" / "stageB"
    stagea_root = Path(args.stagea_root).resolve() if args.stagea_root else project_root / "Input" / "stageA"
    output_dir = Path(args.output_dir).resolve() if args.output_dir else SCRIPT_DIR / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    ratios, thicknesses, records, placement_records = collect_stageb_metrics(
        stageb_root,
        args.birks_kb,
        args.effective_threshold_kev,
    )
    stagea_absorption = read_stagea_absorption(stagea_root)
    add_initial_strength(records, stagea_absorption)

    metric_fields = [
        "ratio_tag",
        "bn_wt",
        "zns_wt",
        "thickness_um",
        "anchor_rows",
        "zns_step_rows",
        "n_physical_captures",
        "capture_weight_sum",
        "n_effective_physical_captures",
        "effective_deposition_event_fraction",
        "effective_deposition_fraction_weighted",
        "mean_raw_zns_edep_per_capture_keV",
        "mean_effective_zns_edep_per_capture_keV",
        "median_effective_zns_edep_per_capture_keV",
        "q25_effective_zns_edep_per_capture_keV",
        "q75_effective_zns_edep_per_capture_keV",
        "weighted_raw_zns_edep_sum_keV",
        "weighted_effective_zns_edep_sum_keV",
        "n_placements",
        "placement_mean_effective_edep_std_keV",
        "placement_effective_fraction_std",
        "stageA_absorption_rate",
        "initial_scintillation_strength_absorption_x_edep",
        "relative_initial_scintillation_strength",
    ]
    placement_fields = [
        "ratio_tag",
        "bn_wt",
        "zns_wt",
        "thickness_um",
        "placement_file",
        "n_anchor_rows",
        "capture_weight_sum",
        "effective_deposition_fraction_weighted",
        "mean_raw_zns_edep_per_capture_keV",
        "mean_effective_zns_edep_per_capture_keV",
    ]

    metrics_csv = output_dir / "stageB_effective_deposition_metrics_by_ratio_and_thickness.csv"
    placements_csv = output_dir / "stageB_effective_deposition_metrics_by_placement.csv"
    write_csv(metrics_csv, records, metric_fields)
    write_csv(placements_csv, placement_records, placement_fields)

    fraction_matrix = records_to_matrix(ratios, thicknesses, records, "effective_deposition_fraction_weighted")
    event_fraction_matrix = records_to_matrix(ratios, thicknesses, records, "effective_deposition_event_fraction")
    edep_matrix = records_to_matrix(ratios, thicknesses, records, "mean_effective_zns_edep_per_capture_keV")
    strength_matrix = records_to_matrix(ratios, thicknesses, records, "relative_initial_scintillation_strength")
    placement_edep_std_matrix = records_to_matrix(
        ratios,
        thicknesses,
        records,
        "placement_mean_effective_edep_std_keV",
    )

    plot_metric_lines(
        ratios,
        thicknesses,
        fraction_matrix,
        output_dir / "fig3a_effective_deposition_fraction_vs_thickness.png",
        "Weighted fraction with ZnS(Ag) edep",
        "Fig. 3(a) Effective ZnS(Ag) deposition fraction",
        y_min=0.0,
        y_max=1.0,
    )
    plot_metric_heatmap(
        ratios,
        thicknesses,
        fraction_matrix,
        output_dir / "fig3a_effective_deposition_fraction_heatmap.png",
        "Weighted fraction",
        "Fig. 3(a) Effective ZnS(Ag) deposition fraction",
        vmin=0.0,
        vmax=min(1.0, max(1.0e-12, float(np.nanmax(fraction_matrix)))),
    )
    plot_metric_lines(
        ratios,
        thicknesses,
        edep_matrix,
        output_dir / "fig3b_mean_effective_edep_vs_thickness.png",
        "Mean effective ZnS(Ag) edep per capture (keV)",
        "Fig. 3(b) Mean effective ZnS(Ag) edep per capture",
        y_min=0.0,
        placement_std_matrix=placement_edep_std_matrix,
    )
    plot_ratio_bar_with_dispersion(
        ratios,
        records,
        placement_records,
        output_dir / "fig3b_mean_effective_edep_by_ratio_with_placement_dispersion.png",
    )
    plot_metric_lines(
        ratios,
        thicknesses,
        strength_matrix,
        output_dir / "fig3c_relative_initial_scintillation_strength_vs_thickness.png",
        "Relative initial scintillation strength",
        "Fig. 3(c) Relative initial scintillation-source strength",
        y_min=0.0,
        y_max=1.05,
    )
    plot_metric_heatmap(
        ratios,
        thicknesses,
        strength_matrix,
        output_dir / "fig3c_relative_initial_scintillation_strength_heatmap.png",
        "Relative strength",
        "Fig. 3(c) Relative initial scintillation-source strength",
        vmin=0.0,
        vmax=1.0,
    )
    plot_combined_fig3(
        ratios,
        thicknesses,
        fraction_matrix,
        edep_matrix,
        strength_matrix,
        output_dir / "fig3_stageB_scintillation_metrics.png",
    )

    # Keep the physical-event fraction in the CSV and write a small diagnostic plot
    # so the replay-weighted and grouped-event definitions can be compared directly.
    plot_metric_lines(
        ratios,
        thicknesses,
        event_fraction_matrix,
        output_dir / "diagnostic_physical_event_effective_fraction_vs_thickness.png",
        "Physical-event fraction with ZnS(Ag) edep",
        "Diagnostic physical-event effective fraction",
        y_min=0.0,
        y_max=min(1.0, max(1.0e-12, float(np.nanmax(event_fraction_matrix)) * 1.12)),
    )

    print(f"Wrote: {metrics_csv}")
    print(f"Wrote: {placements_csv}")
    print(f"Wrote: {output_dir / 'fig3_stageB_scintillation_metrics.png'}")
    print(f"Wrote: {output_dir / 'fig3a_effective_deposition_fraction_vs_thickness.png'}")
    print(f"Wrote: {output_dir / 'fig3b_mean_effective_edep_vs_thickness.png'}")
    print(f"Wrote: {output_dir / 'fig3c_relative_initial_scintillation_strength_vs_thickness.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
