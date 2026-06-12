#!/usr/bin/env python3
"""Generate Stage A absorption figures."""

from __future__ import annotations

import csv
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

DEPENDENCY_ERROR = """
Failed to import the plotting dependencies.

This usually means NumPy and Matplotlib were installed from incompatible
environments. For example, an old system Matplotlib compiled against NumPy 1.x
cannot run with NumPy 2.x.

Use a clean virtual environment from the project root:

    python3 -m venv .venv
    . .venv/bin/activate
    python3 -m pip install -r plots/stageA_capture_heatmaps/requirements.txt
    python3 plots/stageA_capture_heatmaps/generate_stageA_capture_heatmaps.py
""".strip()

try:
    import matplotlib

    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import font_manager
    from matplotlib.colors import LinearSegmentedColormap, PowerNorm
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
except ImportError as exc:
    raise SystemExit(f"{DEPENDENCY_ERROR}\n\nOriginal error: {exc}") from None
except AttributeError as exc:
    if "_ARRAY_API" in str(exc):
        raise SystemExit(f"{DEPENDENCY_ERROR}\n\nOriginal error: {exc}") from None
    raise

SCRIPT_DIR = Path(__file__).resolve().parent
FONT_DIR = SCRIPT_DIR / "fonts"


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
            "legend.fontsize": 11,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "axes.unicode_minus": False,
        }
    )


configure_fonts()


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
CAPTURE_CSV_RE = re.compile(
    r"^(?P<thickness>[+]?(?:\d+(?:\.\d*)?|\.\d+))_neutron_capture_positions\.csv$"
)


@dataclass(frozen=True)
class RatioKey:
    tag: str
    bn_wt: float
    zns_wt: float

    @property
    def bn_over_zns(self) -> float:
        return self.bn_wt / self.zns_wt

    @property
    def display_tag(self) -> str:
        return self.tag.replace("-", ":")


def parse_ratio_tag(tag: str) -> RatioKey:
    match = RATIO_RE.fullmatch(tag)
    if match is None:
        raise ValueError(f"Unsupported ratio folder name: {tag}")
    return RatioKey(
        tag=tag,
        bn_wt=float(match.group("bn")),
        zns_wt=float(match.group("zns")),
    )


def parse_thickness_from_name(path: Path) -> float:
    match = CAPTURE_CSV_RE.fullmatch(path.name)
    if match is None:
        raise ValueError(f"Unsupported capture CSV name: {path.name}")
    return float(match.group("thickness"))


def thickness_sort_key(value: float) -> tuple[int, float]:
    rounded = round(value)
    is_int_like = math.isclose(value, rounded, rel_tol=0.0, abs_tol=1.0e-9)
    return (0 if is_int_like else 1, value)


def ratio_display_sort_key(ratio: RatioKey) -> tuple[int, float, float]:
    try:
        return (DESIRED_RATIO_ORDER.index(ratio.tag), ratio.bn_wt, ratio.zns_wt)
    except ValueError:
        return (len(DESIRED_RATIO_ORDER), ratio.bn_wt, ratio.zns_wt)


def style_axes(ax) -> None:
    ax.tick_params(direction="in", which="both", top=True, right=True, width=1.4)
    ax.set_box_aspect(1)
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)


def mu_m_text() -> str:
    return "μm"


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_overlay_points_from_xlsx(path: Path) -> dict[str, list[tuple[float, float]]]:
    if not path.is_file():
        return {}

    ns = {
        "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
        "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
        "pr": "http://schemas.openxmlformats.org/package/2006/relationships",
    }

    def col_to_idx(ref: str) -> int:
        letters = ""
        for ch in ref:
            if ch.isalpha():
                letters += ch
            else:
                break
        idx = 0
        for ch in letters:
            idx = idx * 26 + (ord(ch.upper()) - 64)
        return idx - 1

    overlay: dict[str, list[tuple[float, float]]] = {}
    with ZipFile(path) as zf:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall("a:si", ns):
                texts = [t.text or "" for t in si.iterfind(".//a:t", ns)]
                shared.append("".join(texts))

        workbook = ET.fromstring(zf.read("xl/workbook.xml"))
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rel_map = {
            rel.attrib["Id"]: rel.attrib["Target"]
            for rel in rels.findall("pr:Relationship", ns)
        }
        first_sheet = workbook.find("a:sheets", ns)[0]
        rel_id = first_sheet.attrib[
            "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
        ]
        sheet_path = "xl/" + rel_map[rel_id]
        sheet = ET.fromstring(zf.read(sheet_path))
        data = sheet.find("a:sheetData", ns)

        rows: list[list[str]] = []
        for row in data.findall("a:row", ns):
            cells = []
            max_idx = -1
            for cell in row.findall("a:c", ns):
                ref = cell.attrib.get("r", "A1")
                idx = col_to_idx(ref)
                max_idx = max(max_idx, idx)
                cell_type = cell.attrib.get("t")
                value = ""
                value_node = cell.find("a:v", ns)
                if value_node is not None and value_node.text is not None:
                    if cell_type == "s":
                        value = shared[int(value_node.text)]
                    else:
                        value = value_node.text
                cells.append((idx, value))
            if cells:
                arr = [""] * (max_idx + 1)
                for idx, value in cells:
                    arr[idx] = value
                rows.append(arr)

    if not rows:
        return overlay

    header = rows[0]
    for col in range(0, len(header), 2):
        ratio_tag = header[col].strip() if col < len(header) else ""
        if not ratio_tag:
            continue
        points: list[tuple[float, float]] = []
        for row in rows[1:]:
            if col + 1 >= len(row):
                continue
            thickness_raw = row[col].strip()
            rate_raw = row[col + 1].strip()
            if not thickness_raw or not rate_raw:
                continue
            points.append((float(thickness_raw), float(rate_raw)))
        if points:
            overlay[ratio_tag] = sorted(points, key=lambda item: item[0])

    return overlay


def make_axis_edges(values: list[float]) -> np.ndarray:
    coords = np.asarray(values, dtype=float)
    if coords.size == 1:
        half_width = max(abs(coords[0]) * 0.5, 0.5)
        return np.asarray([coords[0] - half_width, coords[0] + half_width], dtype=float)
    deltas = np.diff(coords)
    edges = np.empty(coords.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (coords[:-1] + coords[1:])
    edges[0] = coords[0] - 0.5 * deltas[0]
    edges[-1] = coords[-1] + 0.5 * deltas[-1]
    return edges


def make_log_axis_edges(values: list[float]) -> np.ndarray:
    coords = np.asarray(values, dtype=float)
    if np.any(coords <= 0.0):
        raise ValueError("Log-axis coordinates must be positive.")
    if coords.size == 1:
        return np.asarray([coords[0] / 1.5, coords[0] * 1.5], dtype=float)
    logs = np.log10(coords)
    edge_logs = np.empty(coords.size + 1, dtype=float)
    edge_logs[1:-1] = 0.5 * (logs[:-1] + logs[1:])
    edge_logs[0] = logs[0] - 0.5 * (logs[1] - logs[0])
    edge_logs[-1] = logs[-1] + 0.5 * (logs[-1] - logs[-2])
    return np.power(10.0, edge_logs)


def monotone_cubic_interpolate(
    x_values: list[float] | np.ndarray,
    y_values: list[float] | np.ndarray,
    x_dense: np.ndarray,
) -> np.ndarray:
    x = np.asarray(x_values, dtype=float)
    y = np.asarray(y_values, dtype=float)
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    if x.size < 2:
        raise ValueError("Need at least two points.")

    h = np.diff(x)
    delta = np.diff(y) / h
    m = np.zeros_like(y)

    if x.size == 2:
        m[:] = delta[0]
    else:
        for idx in range(1, x.size - 1):
            if delta[idx - 1] == 0.0 or delta[idx] == 0.0:
                m[idx] = 0.0
            elif np.sign(delta[idx - 1]) != np.sign(delta[idx]):
                m[idx] = 0.0
            else:
                w1 = 2.0 * h[idx] + h[idx - 1]
                w2 = h[idx] + 2.0 * h[idx - 1]
                m[idx] = (w1 + w2) / ((w1 / delta[idx - 1]) + (w2 / delta[idx]))

        m[0] = delta[0]
        m[-1] = delta[-1]

    indices = np.searchsorted(x, x_dense, side="right") - 1
    indices = np.clip(indices, 0, x.size - 2)

    x0 = x[indices]
    x1 = x[indices + 1]
    y0 = y[indices]
    y1 = y[indices + 1]
    m0 = m[indices]
    m1 = m[indices + 1]
    h_seg = x1 - x0
    t = (x_dense - x0) / h_seg

    h00 = 2.0 * t**3 - 3.0 * t**2 + 1.0
    h10 = t**3 - 2.0 * t**2 + t
    h01 = -2.0 * t**3 + 3.0 * t**2
    h11 = t**3 - t**2
    return h00 * y0 + h10 * h_seg * m0 + h01 * y1 + h11 * h_seg * m1


def collect_absorption_and_depth(
    stagea_root: Path,
) -> tuple[list[RatioKey], list[float], np.ndarray, np.ndarray, list[dict[str, object]], list[dict[str, object]]]:
    ratio_dirs = []
    for path in stagea_root.iterdir():
        if path.is_dir() and (path / "neutron_capture_positions").is_dir():
            ratio_dirs.append(parse_ratio_tag(path.name))
    ratios = sorted(ratio_dirs, key=ratio_display_sort_key)

    thickness_values: set[float] = set()
    absorption_records: list[dict[str, object]] = []
    mean_depth_records: list[dict[str, object]] = []
    absorption_maps: dict[str, dict[float, float]] = {}
    depth_maps: dict[str, dict[float, float]] = {}

    for ratio in ratios:
        capture_dir = stagea_root / ratio.tag / "neutron_capture_positions"
        absorption_maps[ratio.tag] = {}
        depth_maps[ratio.tag] = {}

        for capture_csv in sorted(capture_dir.glob("*_neutron_capture_positions.csv")):
            thickness_um = parse_thickness_from_name(capture_csv)
            thickness_values.add(thickness_um)

            count = 0
            max_event_id = -1
            depth_sum = 0.0

            with capture_csv.open(newline="", encoding="utf-8-sig") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    count += 1
                    event_id = int(float(row["eventID"]))
                    depth_um = float(row["depth_um"])
                    max_event_id = max(max_event_id, event_id)
                    depth_sum += depth_um

            if count == 0 or max_event_id < 0:
                continue

            incident_count = max_event_id + 1
            absorption_rate = count / incident_count
            mean_depth_um = depth_sum / count

            absorption_maps[ratio.tag][thickness_um] = absorption_rate
            depth_maps[ratio.tag][thickness_um] = mean_depth_um

            absorption_records.append(
                {
                    "ratio_tag": ratio.tag,
                    "bn_wt": ratio.bn_wt,
                    "zns_wt": ratio.zns_wt,
                    "thickness_um": thickness_um,
                    "incident_count_inferred": incident_count,
                    "n_absorbed": count,
                    "absorption_rate": absorption_rate,
                }
            )
            mean_depth_records.append(
                {
                    "ratio_tag": ratio.tag,
                    "bn_wt": ratio.bn_wt,
                    "zns_wt": ratio.zns_wt,
                    "thickness_um": thickness_um,
                    "mean_absorption_depth_um": mean_depth_um,
                    "n_absorbed": count,
                }
            )

    thicknesses = sorted(thickness_values, key=thickness_sort_key)
    absorption_matrix = np.full((len(ratios), len(thicknesses)), np.nan, dtype=float)
    mean_depth_matrix = np.full((len(ratios), len(thicknesses)), np.nan, dtype=float)

    ratio_index = {ratio.tag: idx for idx, ratio in enumerate(ratios)}
    thickness_index = {value: idx for idx, value in enumerate(thicknesses)}
    for ratio_tag, values in absorption_maps.items():
        for thickness_um, absorption_rate in values.items():
            absorption_matrix[ratio_index[ratio_tag], thickness_index[thickness_um]] = absorption_rate
    for ratio_tag, values in depth_maps.items():
        for thickness_um, mean_depth_um in values.items():
            mean_depth_matrix[ratio_index[ratio_tag], thickness_index[thickness_um]] = mean_depth_um

    return ratios, thicknesses, absorption_matrix, mean_depth_matrix, absorption_records, mean_depth_records


def collect_normalized_depth_panel_data(
    stagea_root: Path,
    ratios: list[RatioKey],
    thicknesses: list[float],
    n_bins: int = 80,
) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, np.ndarray], list[dict[str, object]]]:
    z_edges = np.linspace(0.0, 1.0, n_bins + 1)
    density_by_ratio: dict[str, np.ndarray] = {}
    mean_by_ratio: dict[str, np.ndarray] = {}
    records: list[dict[str, object]] = []

    for ratio in ratios:
        capture_dir = stagea_root / ratio.tag / "neutron_capture_positions"
        density_matrix = np.full((len(thicknesses), n_bins), np.nan, dtype=float)
        mean_values = np.full(len(thicknesses), np.nan, dtype=float)

        for thickness_idx, thickness_um in enumerate(thicknesses):
            capture_csv = capture_dir / f"{thickness_um:g}_neutron_capture_positions.csv"
            if not capture_csv.is_file():
                continue

            normalized_depths: list[float] = []
            with capture_csv.open(newline="", encoding="utf-8-sig") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    depth_um = float(row["depth_um"])
                    normalized_depths.append(min(1.0, max(0.0, depth_um / thickness_um)))

            if not normalized_depths:
                continue

            mean_values[thickness_idx] = float(np.mean(normalized_depths))
            density, _ = np.histogram(normalized_depths, bins=z_edges, density=True)
            density_matrix[thickness_idx, :] = density

            total_count = len(normalized_depths)
            for bin_idx, density_value in enumerate(density):
                records.append(
                    {
                        "ratio_tag": ratio.tag,
                        "thickness_um": thickness_um,
                        "z_over_t_bin_left": z_edges[bin_idx],
                        "z_over_t_bin_right": z_edges[bin_idx + 1],
                        "z_over_t_bin_center": 0.5 * (z_edges[bin_idx] + z_edges[bin_idx + 1]),
                        "density": float(density_value),
                        "n_absorbed": total_count,
                    }
                )

        density_by_ratio[ratio.tag] = density_matrix
        mean_by_ratio[ratio.tag] = mean_values

    return z_edges, density_by_ratio, mean_by_ratio, records


def collect_absolute_depth_panel_data(
    stagea_root: Path,
    ratios: list[RatioKey],
    thicknesses: list[float],
    n_bins: int = 80,
) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, np.ndarray], list[dict[str, object]]]:
    max_thickness = float(max(thicknesses))
    z_edges = np.linspace(0.0, max_thickness, n_bins + 1)
    density_by_ratio: dict[str, np.ndarray] = {}
    mean_by_ratio: dict[str, np.ndarray] = {}
    records: list[dict[str, object]] = []

    for ratio in ratios:
        capture_dir = stagea_root / ratio.tag / "neutron_capture_positions"
        density_matrix = np.full((len(thicknesses), n_bins), np.nan, dtype=float)
        mean_values = np.full(len(thicknesses), np.nan, dtype=float)

        for thickness_idx, thickness_um in enumerate(thicknesses):
            capture_csv = capture_dir / f"{thickness_um:g}_neutron_capture_positions.csv"
            if not capture_csv.is_file():
                continue

            depths_um: list[float] = []
            with capture_csv.open(newline="", encoding="utf-8-sig") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    depths_um.append(float(row["depth_um"]))

            if not depths_um:
                continue

            mean_values[thickness_idx] = float(np.mean(depths_um))
            density, _ = np.histogram(depths_um, bins=z_edges, density=True)
            density_matrix[thickness_idx, :] = density

            total_count = len(depths_um)
            for bin_idx, density_value in enumerate(density):
                records.append(
                    {
                        "ratio_tag": ratio.tag,
                        "thickness_um": thickness_um,
                        "z_um_bin_left": z_edges[bin_idx],
                        "z_um_bin_right": z_edges[bin_idx + 1],
                        "z_um_bin_center": 0.5 * (z_edges[bin_idx] + z_edges[bin_idx + 1]),
                        "density_per_um": float(density_value),
                        "n_absorbed": total_count,
                    }
                )

        density_by_ratio[ratio.tag] = density_matrix
        mean_by_ratio[ratio.tag] = mean_values

    return z_edges, density_by_ratio, mean_by_ratio, records


def collect_absolute_depth_summary_panel_data(
    stagea_root: Path,
    ratios: list[RatioKey],
    thicknesses: list[float],
) -> tuple[dict[str, dict[str, np.ndarray]], list[dict[str, object]]]:
    stats_by_ratio: dict[str, dict[str, np.ndarray]] = {}
    records: list[dict[str, object]] = []

    for ratio in ratios:
        capture_dir = stagea_root / ratio.tag / "neutron_capture_positions"
        stats = {
            "mean": np.full(len(thicknesses), np.nan, dtype=float),
            "q10": np.full(len(thicknesses), np.nan, dtype=float),
            "q25": np.full(len(thicknesses), np.nan, dtype=float),
            "q50": np.full(len(thicknesses), np.nan, dtype=float),
            "q75": np.full(len(thicknesses), np.nan, dtype=float),
            "q90": np.full(len(thicknesses), np.nan, dtype=float),
        }

        for thickness_idx, thickness_um in enumerate(thicknesses):
            capture_csv = capture_dir / f"{thickness_um:g}_neutron_capture_positions.csv"
            if not capture_csv.is_file():
                continue

            depths_um: list[float] = []
            with capture_csv.open(newline="", encoding="utf-8-sig") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    depths_um.append(float(row["depth_um"]))

            if not depths_um:
                continue

            values = np.asarray(depths_um, dtype=float)
            stats["mean"][thickness_idx] = float(np.mean(values))
            stats["q10"][thickness_idx] = float(np.quantile(values, 0.10))
            stats["q25"][thickness_idx] = float(np.quantile(values, 0.25))
            stats["q50"][thickness_idx] = float(np.quantile(values, 0.50))
            stats["q75"][thickness_idx] = float(np.quantile(values, 0.75))
            stats["q90"][thickness_idx] = float(np.quantile(values, 0.90))

            records.append(
                {
                    "ratio_tag": ratio.tag,
                    "thickness_um": thickness_um,
                    "mean_absorption_depth_um": stats["mean"][thickness_idx],
                    "q10_absorption_depth_um": stats["q10"][thickness_idx],
                    "q25_absorption_depth_um": stats["q25"][thickness_idx],
                    "q50_absorption_depth_um": stats["q50"][thickness_idx],
                    "q75_absorption_depth_um": stats["q75"][thickness_idx],
                    "q90_absorption_depth_um": stats["q90"][thickness_idx],
                    "n_absorbed": int(values.size),
                }
            )

        stats_by_ratio[ratio.tag] = stats

    return stats_by_ratio, records


def interpolate_rows_along_log_thickness(
    matrix: np.ndarray,
    thicknesses: list[float],
    sample_count: int = 1000,
) -> tuple[np.ndarray, np.ndarray]:
    x_known = np.asarray(thicknesses, dtype=float)
    x_dense = np.geomspace(float(np.min(x_known)), float(np.max(x_known)), sample_count)
    dense = np.full((matrix.shape[0], sample_count), np.nan, dtype=float)
    log_x_known = np.log10(x_known)
    log_x_dense = np.log10(x_dense)

    for row_index in range(matrix.shape[0]):
        row = matrix[row_index, :]
        valid = np.isfinite(row)
        if np.count_nonzero(valid) < 2:
            continue
        dense[row_index, :] = np.interp(log_x_dense, log_x_known[valid], row[valid])
    return x_dense, dense


def interpolate_panel_along_thickness(
    matrix: np.ndarray,
    thicknesses: list[float],
    sample_count: int = 1200,
) -> tuple[np.ndarray, np.ndarray]:
    x_known = np.asarray(thicknesses, dtype=float)
    x_dense = np.linspace(float(np.min(x_known)), float(np.max(x_known)), sample_count)
    dense = np.full((sample_count, matrix.shape[1]), np.nan, dtype=float)

    for col_index in range(matrix.shape[1]):
        column = matrix[:, col_index]
        valid = np.isfinite(column)
        if np.count_nonzero(valid) < 2:
            continue
        dense[:, col_index] = np.interp(x_dense, x_known[valid], column[valid])

    return x_dense, dense


def interpolate_panel_field(
    matrix: np.ndarray,
    x_values: list[float],
    y_centers: np.ndarray,
    x_samples: int = 1200,
    y_samples: int = 500,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_dense = np.linspace(float(min(x_values)), float(max(x_values)), x_samples)
    x_interp, dense_x = interpolate_panel_along_thickness(
        matrix,
        x_values,
        sample_count=x_samples,
    )
    y_dense = np.linspace(float(np.min(y_centers)), float(np.max(y_centers)), y_samples)
    field = np.full((y_samples, x_samples), np.nan, dtype=float)

    for x_idx in range(x_samples):
        row = dense_x[x_idx, :]
        valid = np.isfinite(row)
        valid_count = int(np.count_nonzero(valid))
        if valid_count >= 2:
            field[:, x_idx] = np.interp(y_dense, y_centers[valid], row[valid])
        elif valid_count == 1:
            field[:, x_idx] = row[valid][0]

    return x_interp, y_dense, field


def plot_absorption_rate_lines(
    ratios: list[RatioKey],
    thicknesses: list[float],
    matrix: np.ndarray,
    overlay_points: dict[str, list[tuple[float, float]]],
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    ratio_handles: list[Line2D] = []
    overlay_handle = None
    overlay_markers = ["o", "s", "^", "D", "P", "v"]

    for idx, ratio in enumerate(ratios):
        row = matrix[idx, :]
        valid = np.isfinite(row)
        if np.count_nonzero(valid) < 2:
            continue
        x = np.asarray(thicknesses, dtype=float)[valid]
        y = row[valid]
        x_dense = np.linspace(x.min(), x.max(), 3000)
        y_dense = monotone_cubic_interpolate(x, y, x_dense)
        color = DISCRETE_RATIO_COLORS[idx % len(DISCRETE_RATIO_COLORS)]
        ax.plot(
            x_dense,
            y_dense,
            color=color,
            linewidth=2.0,
            alpha=0.72,
            solid_capstyle="round",
            solid_joinstyle="round",
            antialiased=True,
        )
        ratio_handles.append(
            Line2D([0], [0], color=color, linewidth=2.0, alpha=0.72, label=ratio.display_tag)
        )
        if ratio.tag in overlay_points:
            points = np.asarray(overlay_points[ratio.tag], dtype=float)
            ax.scatter(
                points[:, 0],
                points[:, 1],
                s=78,
                marker=overlay_markers[idx % len(overlay_markers)],
                facecolor="white",
                edgecolor=color,
                linewidth=2.0,
                zorder=6,
            )
            if overlay_handle is None:
                overlay_handle = Line2D(
                    [0],
                    [0],
                    linestyle="none",
                    marker="o",
                    markersize=8.5,
                    markerfacecolor="white",
                    markeredgecolor="#1f1f1f",
                    markeredgewidth=1.8,
                    label="Experimental data",
                )

    ax.set_xlim(0.0, max(thicknesses))
    ax.set_ylim(0.0, min(1.0, float(np.nanmax(matrix)) * 1.08))
    ax.set_xticks(np.arange(0, int(max(thicknesses)) + 1, 100))
    ax.set_xlabel(f"Thickness t ({mu_m_text()})")
    ax.set_ylabel("Neutron absorption rate")
    ax.set_title("Neutron absorption rate vs thickness")
    ax.grid(True, alpha=0.18, linewidth=0.7)
    legend_handles = ratio_handles + ([overlay_handle] if overlay_handle is not None else [])
    ax.legend(handles=legend_handles, frameon=False, fontsize=10, ncol=2)
    style_axes(ax)
    ax.set_box_aspect(2 / 3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_absorption_rate_heatmap(
    ratios: list[RatioKey],
    thicknesses: list[float],
    matrix: np.ndarray,
    output_path: Path,
) -> None:
    x_dense, matrix_dense = interpolate_rows_along_log_thickness(matrix, thicknesses)
    x_edges = make_log_axis_edges(list(x_dense))
    gap = 0.18
    bar_height = 0.82

    fig, ax = plt.subplots(figsize=(8.0, 7.0))
    cmap = RAINBOW_CMAP.copy()
    cmap.set_bad(color="#ececec")
    vmin = 0.0
    vmax = float(np.nanmax(matrix_dense))
    image = None
    y_centers: list[float] = []

    for idx, ratio in enumerate(ratios):
        y_bottom = idx * (bar_height + gap)
        y_top = y_bottom + bar_height
        y_centers.append(0.5 * (y_bottom + y_top))
        band = matrix_dense[idx : idx + 1, :]
        image = ax.pcolormesh(
            x_edges,
            np.asarray([y_bottom, y_top], dtype=float),
            band,
            shading="flat",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )

    ax.set_xscale("log")
    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(len(ratios) * (bar_height + gap) - gap, 0.0)
    ax.set_xticks([1, 10, 100, 1000])
    ax.set_xticklabels(["1", "10", "100", "1000"])
    ax.set_yticks(y_centers)
    ax.set_yticklabels([ratio.display_tag for ratio in ratios])
    ax.set_xlabel(f"Thickness t ({mu_m_text()})")
    ax.set_ylabel("BN/ZnS(Ag) ratio")
    ax.set_title("Neutron absorption rate heatmap")
    style_axes(ax)

    colorbar = fig.colorbar(image, ax=ax, pad=0.06)
    colorbar.set_label("Absorption rate")
    colorbar.ax.tick_params(direction="in", which="both", width=1.4)
    colorbar.outline.set_linewidth(1.8)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_mean_absorption_depth(
    ratios: list[RatioKey],
    thicknesses: list[float],
    matrix: np.ndarray,
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 7.0))
    valid_values = matrix[np.isfinite(matrix)]
    y_top = max(1.0, float(np.max(valid_values)) * 1.08)

    for idx, ratio in enumerate(ratios):
        row = matrix[idx, :]
        valid = np.isfinite(row)
        if np.count_nonzero(valid) < 2:
            continue
        x = np.asarray(thicknesses, dtype=float)[valid]
        y = row[valid]
        x_dense = np.geomspace(x.min(), x.max(), 1000)
        y_dense = monotone_cubic_interpolate(np.log10(x), y, np.log10(x_dense))
        color = DISCRETE_RATIO_COLORS[idx % len(DISCRETE_RATIO_COLORS)]
        ax.scatter(x, y, s=18, color=color, alpha=0.92)
        ax.plot(x_dense, y_dense, color=color, linewidth=2.0, label=ratio.display_tag)

    ax.set_xscale("log")
    ax.set_xlim(min(thicknesses), max(thicknesses))
    ax.set_ylim(0.0, y_top)
    ax.set_xlabel(f"Thickness t ({mu_m_text()})")
    ax.set_ylabel(f"Mean absorption depth ({mu_m_text()})")
    ax.set_title("Mean absorption depth vs thickness")
    ax.grid(True, which="both", alpha=0.18, linewidth=0.7)
    ax.legend(frameon=False, fontsize=10, ncol=2)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_normalized_depth_distribution_panels(
    ratios: list[RatioKey],
    thicknesses: list[float],
    z_edges: np.ndarray,
    density_by_ratio: dict[str, np.ndarray],
    mean_by_ratio: dict[str, np.ndarray],
    output_path: Path,
) -> None:
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    vmax = max(
        float(np.nanmax(matrix))
        for matrix in density_by_ratio.values()
        if np.isfinite(matrix).any()
    )

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(12.4, 7.8),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    panel_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
    mappable = None

    for idx, (ax, ratio) in enumerate(zip(axes.flat, ratios)):
        x_dense, z_dense, field = interpolate_panel_field(
            density_by_ratio[ratio.tag],
            thicknesses,
            z_centers,
        )
        mappable = ax.imshow(
            field,
            origin="lower",
            extent=[x_dense[0], x_dense[-1], z_dense[0], z_dense[-1]],
            aspect="auto",
            cmap=RAINBOW_CMAP,
            vmin=0.0,
            vmax=vmax,
            interpolation="bilinear",
        )

        mean_values = mean_by_ratio[ratio.tag]
        valid_mean = np.isfinite(mean_values)
        if np.count_nonzero(valid_mean) >= 2:
            x_line = np.asarray(thicknesses, dtype=float)[valid_mean]
            y_line = mean_values[valid_mean]
            x_line_dense = np.linspace(x_line.min(), x_line.max(), 900)
            y_line_dense = monotone_cubic_interpolate(x_line, y_line, x_line_dense)
            ax.plot(x_line_dense, y_line_dense, color="white", linewidth=2.2, zorder=5)
            ax.plot(x_line_dense, y_line_dense, color="black", linewidth=0.8, zorder=6)

        ax.set_xlim(min(thicknesses), max(thicknesses))
        ax.set_ylim(0.0, 1.0)
        ax.set_xticks([0, 100, 200, 400, 700, 1000])
        ax.set_xticklabels(["0", "100", "200", "400", "700", "1000"])
        ax.set_yticks(np.linspace(0.0, 1.0, 6))
        ax.set_title(f"{panel_labels[idx]} BN:ZnS = {ratio.display_tag}", fontsize=11, pad=6)
        ax.tick_params(direction="in", which="both", top=True, right=True, width=1.2)
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)

    for ax in axes[1, :]:
        ax.set_xlabel(f"Thickness t ({mu_m_text()})")
    for ax in axes[:, 0]:
        ax.set_ylabel("Normalized capture depth z/t")

    cbar = fig.colorbar(
        mappable,
        ax=axes,
        pad=0.02,
        shrink=0.92,
        location="right",
    )
    cbar.set_label("Conditional density p(z/t | absorbed)")
    cbar.ax.tick_params(direction="in", which="both", width=1.2)
    cbar.outline.set_linewidth(1.6)

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_absolute_depth_distribution_panels(
    ratios: list[RatioKey],
    thicknesses: list[float],
    stats_by_ratio: dict[str, dict[str, np.ndarray]],
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(12.4, 7.8),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    panel_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
    legend_handles = [
        Patch(facecolor="#9ecae1", edgecolor="none", alpha=0.28, label="10-90%"),
        Patch(facecolor="#3182bd", edgecolor="none", alpha=0.26, label="25-75%"),
        Line2D([0], [0], color="#2d2d2d", linewidth=1.4, label="Median"),
        Line2D([0], [0], color="#d7301f", linewidth=2.0, label="Mean"),
    ]

    for idx, (ax, ratio) in enumerate(zip(axes.flat, ratios)):
        stats = stats_by_ratio[ratio.tag]
        mean_values = stats["mean"]
        q10 = stats["q10"]
        q25 = stats["q25"]
        q50 = stats["q50"]
        q75 = stats["q75"]
        q90 = stats["q90"]
        valid_mean = np.isfinite(mean_values)
        if np.count_nonzero(valid_mean) >= 2:
            x_line = np.asarray(thicknesses, dtype=float)[valid_mean]
            mean_line = mean_values[valid_mean]
            q10_line = q10[valid_mean]
            q25_line = q25[valid_mean]
            q50_line = q50[valid_mean]
            q75_line = q75[valid_mean]
            q90_line = q90[valid_mean]
            x_line_dense = np.linspace(x_line.min(), x_line.max(), 900)
            mean_dense = monotone_cubic_interpolate(x_line, mean_line, x_line_dense)
            q10_dense = monotone_cubic_interpolate(x_line, q10_line, x_line_dense)
            q25_dense = monotone_cubic_interpolate(x_line, q25_line, x_line_dense)
            q50_dense = monotone_cubic_interpolate(x_line, q50_line, x_line_dense)
            q75_dense = monotone_cubic_interpolate(x_line, q75_line, x_line_dense)
            q90_dense = monotone_cubic_interpolate(x_line, q90_line, x_line_dense)

            ax.fill_between(
                x_line_dense,
                q10_dense,
                q90_dense,
                color="#9ecae1",
                alpha=0.28,
                linewidth=0.0,
                zorder=1,
            )
            ax.fill_between(
                x_line_dense,
                q25_dense,
                q75_dense,
                color="#3182bd",
                alpha=0.26,
                linewidth=0.0,
                zorder=2,
            )
            ax.plot(x_line_dense, q50_dense, color="#2d2d2d", linewidth=1.4, zorder=4)
            ax.plot(x_line_dense, mean_dense, color="#d7301f", linewidth=2.0, zorder=5)
            ax.scatter(x_line, mean_line, color="#d7301f", s=10, zorder=6)

        ax.set_xlim(min(thicknesses), max(thicknesses))
        ax.set_ylim(0.0, max(thicknesses))
        ax.set_xticks([0, 100, 200, 400, 700, 1000])
        ax.set_xticklabels(["0", "100", "200", "400", "700", "1000"])
        ax.set_yticks([0, 200, 400, 600, 800, 1000])
        ax.set_title(f"{panel_labels[idx]} BN:ZnS = {ratio.display_tag}", fontsize=11, pad=6)
        ax.tick_params(direction="in", which="both", top=True, right=True, width=1.2)
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
        ax.legend(
            handles=legend_handles,
            loc="upper left",
            frameon=False,
            fontsize=9.5,
            handlelength=2.2,
            borderaxespad=0.35,
            labelspacing=0.35,
        )

    for ax in axes[1, :]:
        ax.set_xlabel(f"Thickness t ({mu_m_text()})")
    for ax in axes[:, 0]:
        ax.set_ylabel(f"Capture depth z ({mu_m_text()})")

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    project_root = Path(__file__).resolve().parent.parent.parent
    stagea_root = project_root / "Input" / "stageA"
    output_dir = Path(__file__).resolve().parent / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    overlay_xlsx = Path(__file__).resolve().parent / "input" / "neutron_absorption.xlsx"

    ratios, thicknesses, absorption_matrix, mean_depth_matrix, absorption_records, mean_depth_records = (
        collect_absorption_and_depth(stagea_root)
    )
    overlay_points = read_overlay_points_from_xlsx(overlay_xlsx)

    absorption_csv = output_dir / "neutron_absorption_rate_by_ratio_and_thickness.csv"
    write_csv(
        absorption_csv,
        absorption_records,
        fieldnames=[
            "ratio_tag",
            "bn_wt",
            "zns_wt",
            "thickness_um",
            "incident_count_inferred",
            "n_absorbed",
            "absorption_rate",
        ],
    )

    mean_depth_csv = output_dir / "mean_absorption_depth_by_ratio_and_thickness.csv"
    write_csv(
        mean_depth_csv,
        mean_depth_records,
        fieldnames=[
            "ratio_tag",
            "bn_wt",
            "zns_wt",
            "thickness_um",
            "mean_absorption_depth_um",
            "n_absorbed",
        ],
    )

    plot_absorption_rate_lines(
        ratios,
        thicknesses,
        absorption_matrix,
        overlay_points,
        output_dir / "neutron_absorption_rate_vs_thickness.png",
    )
    plot_absorption_rate_heatmap(
        ratios,
        thicknesses,
        absorption_matrix,
        output_dir / "neutron_absorption_rate_heatmap.png",
    )
    plot_mean_absorption_depth(
        ratios,
        thicknesses,
        mean_depth_matrix,
        output_dir / "mean_absorption_depth_vs_thickness.png",
    )
    z_edges, density_by_ratio, mean_by_ratio, normalized_depth_records = collect_normalized_depth_panel_data(
        stagea_root,
        ratios,
        thicknesses,
    )
    normalized_depth_csv = output_dir / "normalized_capture_depth_distribution_by_ratio_and_thickness.csv"
    write_csv(
        normalized_depth_csv,
        normalized_depth_records,
        fieldnames=[
            "ratio_tag",
            "thickness_um",
            "z_over_t_bin_left",
            "z_over_t_bin_right",
            "z_over_t_bin_center",
            "density",
            "n_absorbed",
        ],
    )
    plot_normalized_depth_distribution_panels(
        ratios,
        thicknesses,
        z_edges,
        density_by_ratio,
        mean_by_ratio,
        output_dir / "normalized_capture_depth_distribution_panels.png",
    )
    stats_abs_by_ratio, absolute_depth_records = collect_absolute_depth_summary_panel_data(
        stagea_root,
        ratios,
        thicknesses,
    )
    absolute_depth_csv = output_dir / "absolute_capture_depth_distribution_by_ratio_and_thickness.csv"
    write_csv(
        absolute_depth_csv,
        absolute_depth_records,
        fieldnames=[
            "ratio_tag",
            "thickness_um",
            "mean_absorption_depth_um",
            "q10_absorption_depth_um",
            "q25_absorption_depth_um",
            "q50_absorption_depth_um",
            "q75_absorption_depth_um",
            "q90_absorption_depth_um",
            "n_absorbed",
        ],
    )
    plot_absolute_depth_distribution_panels(
        ratios,
        thicknesses,
        stats_abs_by_ratio,
        output_dir / "absolute_capture_depth_distribution_panels.png",
    )

    print(f"Wrote: {absorption_csv}")
    print(f"Wrote: {mean_depth_csv}")
    print(f"Wrote: {normalized_depth_csv}")
    print(f"Wrote: {absolute_depth_csv}")
    print(f"Wrote: {output_dir / 'neutron_absorption_rate_vs_thickness.png'}")
    print(f"Wrote: {output_dir / 'neutron_absorption_rate_heatmap.png'}")
    print(f"Wrote: {output_dir / 'mean_absorption_depth_vs_thickness.png'}")
    print(f"Wrote: {output_dir / 'normalized_capture_depth_distribution_panels.png'}")
    print(f"Wrote: {output_dir / 'absolute_capture_depth_distribution_panels.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
