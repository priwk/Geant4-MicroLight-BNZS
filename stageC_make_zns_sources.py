#!/usr/bin/env python3
import argparse
import csv
import math
import re
from collections import OrderedDict
from pathlib import Path


ALPHA_LI_STEPS_RE = re.compile(
    r"^(?P<thickness>[+]?(?:\d+(?:\.\d*)?|\.\d+))_alpha_li_steps\.csv$"
)
CAPTURE_ANCHORS_RE = re.compile(
    r"^(?P<thickness>[+]?(?:\d+(?:\.\d*)?|\.\d+))_capture_anchors\.csv$"
)
ZNS_TRACK_STEPS_RE = re.compile(
    r"^(?P<thickness>[+]?(?:\d+(?:\.\d*)?|\.\d+))_zns_track_steps\.csv$"
)


STEP_OUTPUT_FIELDS = [
    "physical_event_uid",
    "source_event_uid",
    "source_step_uid",
    "eventID",
    "record_index",
    "trackID",
    "stepID",
    "particle",
    "thickness_um",
    "bn_wt",
    "zns_wt",
    "ratio_label",
    "depth_um",
    "capture_depth_um",
    "placement_file",
    "placement_hash",
    "surface_mode",
    "local_capture_x_um",
    "local_capture_y_um",
    "local_capture_z_um",
    "alphali_replay_index",
    "alphali_replay_count",
    "trajectory_weight",
    "source_model",
    "sampling_model",
    "phase_pre",
    "phase_post",
    "x_pre_um",
    "y_pre_um",
    "z_pre_um",
    "x_post_um",
    "y_post_um",
    "z_post_um",
    "x_mid_um",
    "y_mid_um",
    "z_mid_um",
    "step_len_um",
    "edep_keV",
    "visible_edep_keV",
    "quenching_factor_step",
    "n_photon_step",
]


EVENT_OUTPUT_FIELDS = [
    "physical_event_uid",
    "source_event_uid",
    "eventID",
    "record_index",
    "thickness_um",
    "bn_wt",
    "zns_wt",
    "ratio_label",
    "depth_um",
    "capture_depth_um",
    "capture_x_um",
    "capture_y_um",
    "local_capture_x_um",
    "local_capture_y_um",
    "local_capture_z_um",
    "surface_mode",
    "placement_file",
    "placement_hash",
    "alphali_replay_index",
    "alphali_replay_count",
    "trajectory_weight",
    "n_total_steps",
    "edep_ZnS_keV",
    "edep_total_keV",
    "visible_edep_ZnS_keV",
    "quenching_factor_event",
    "n_photon0",
    "n_zns_steps",
    "has_zns_edep",
]


REQUIRED_INPUT_FIELDS = [
    "eventID",
    "trackID",
    "stepID",
    "particle",
    "thickness_um",
    "bn_wt",
    "zns_wt",
    "depth_um",
    "placement_file",
    "phase_pre",
    "phase_post",
    "x_pre_um",
    "y_pre_um",
    "z_pre_um",
    "x_post_um",
    "y_post_um",
    "z_post_um",
    "step_len_um",
    "edep_keV",
]


REQUIRED_ANCHOR_FIELDS = [
    "physical_event_uid",
    "source_event_uid",
    "eventID",
    "record_index",
    "thickness_um",
    "bn_wt",
    "zns_wt",
    "depth_um",
    "capture_x_um",
    "capture_y_um",
    "local_capture_x_um",
    "local_capture_y_um",
    "local_capture_z_um",
    "surface_mode",
    "placement_file",
    "alphali_replay_index",
    "alphali_replay_count",
    "trajectory_weight",
]


REQUIRED_TRACK_FIELDS = [
    "physical_event_uid",
    "source_event_uid",
    "eventID",
    "record_index",
    "trackID",
    "stepID",
    "particle",
    "phase_post",
    "x_pre_um",
    "y_pre_um",
    "z_pre_um",
    "x_post_um",
    "y_post_um",
    "z_post_um",
    "step_len_um",
    "edep_keV",
    "ekin_pre_keV",
    "ekin_post_keV",
    "alphali_replay_index",
    "alphali_replay_count",
    "trajectory_weight",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Create Stage C ZnS(Ag) scintillation source CSVs from an existing "
            "Stage B alpha_li_steps.csv without modifying the input file."
        )
    )
    parser.add_argument(
        "input_csv",
        nargs="?",
        help="Existing Stage B alpha_li_steps.csv file. This file is opened read-only.",
    )
    parser.add_argument(
        "--project-root",
        default=None,
        help="Project root. Defaults to the directory containing this script.",
    )
    parser.add_argument(
        "--ratio-tag",
        default=None,
        help="Batch mode ratio folder, such as 1-1, 1-1.5, or 1-2.",
    )
    parser.add_argument(
        "--thickness-range",
        default=None,
        help="Batch mode thickness range in um, such as 30-200.",
    )
    parser.add_argument(
        "--batch-output-root",
        default=None,
        help=(
            "Batch mode output root. Defaults to "
            "<project-root>/Input/alpha_li_steps."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for derived CSVs. Defaults to the input CSV directory.",
    )
    parser.add_argument(
        "--step-output",
        default="stageB_zns_step_sources.csv",
        help="Filename or path for ZnS step-level source output.",
    )
    parser.add_argument(
        "--event-output",
        default="stageB_event_light_sources.csv",
        help="Filename or path for event-level light-source summary output.",
    )
    parser.add_argument(
        "--yield-zns",
        type=float,
        default=55000.0,
        help="ZnS(Ag) scintillation yield in photons/MeV. Default: 55000.",
    )
    parser.add_argument(
        "--birks-kb",
        type=float,
        default=0.0,
        help=(
            "Birks kB in um/keV when dE/dx is computed in keV/um. "
            "Use 0 to disable quenching."
        ),
    )
    parser.add_argument(
        "--split-by-placement",
        action="store_true",
        help=(
            "Write one zns_step_sources/event_light_sources pair per placement_file. "
            "This is the recommended Stage C optical input because each Geant4 run "
            "uses one RVE placement. Batch mode writes them under "
            "<output>/<ratio>/by_placement/<thickness>/."
        ),
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help=(
            "Prefix for split output filenames in single-file mode. Defaults to the "
            "input filename with _alpha_li_steps removed when present."
        ),
    )
    return parser.parse_args()


def parse_thickness_range(value):
    match = re.fullmatch(
        r"\s*(?P<start>[+]?(?:\d+(?:\.\d*)?|\.\d+))\s*-\s*"
        r"(?P<end>[+]?(?:\d+(?:\.\d*)?|\.\d+))\s*",
        value,
    )
    if match is None:
        raise SystemExit(
            "Invalid thickness range. Use start-end, for example: 30-200"
        )

    start = float(match.group("start"))
    end = float(match.group("end"))
    if start > end:
        raise SystemExit("Invalid thickness range: start must be <= end.")
    return start, end


def prompt_required(prompt):
    while True:
        value = input(prompt).strip()
        if value:
            return value
        print("Value cannot be empty.")


def resolve_output_path(value, output_dir):
    path = Path(value)
    if path.is_absolute():
        return path
    return output_dir / path


def require_fields(fieldnames, required):
    if fieldnames is None:
        raise SystemExit("Input CSV is empty or missing a header row.")

    missing = [name for name in required if name not in fieldnames]
    if missing:
        raise SystemExit(
            "Input CSV is missing required column(s): " + ", ".join(missing)
        )


def to_float(row, name, line_number):
    value = row.get(name, "")
    try:
        number = float(value)
    except ValueError as exc:
        raise SystemExit(
            f"Invalid numeric value at CSV line {line_number}, column {name}: {value!r}"
        ) from exc
    if not math.isfinite(number):
        raise SystemExit(
            f"Non-finite numeric value at CSV line {line_number}, column {name}: {value!r}"
        )
    return number


def fmt(number):
    return f"{number:.12g}"


def safe_name(value):
    value = Path(value.strip()).stem or "unknown_placement"
    value = re.sub(r"[^A-Za-z0-9._+-]+", "_", value)
    return value.strip("_") or "unknown_placement"


def ratio_label_from_row(row):
    return f"{row.get('bn_wt', '').strip()}-{row.get('zns_wt', '').strip()}"


def placement_hash(placement_file):
    return safe_name(placement_file)


def event_uid(row):
    source = str(row.get("source_event_uid", "")).strip()
    if source:
        return source
    return "|".join(
        [
            ratio_label_from_row(row),
            str(row.get("thickness_um", "")).strip(),
            placement_hash(row.get("placement_file", "")),
            str(row.get("eventID", "")).strip(),
            str(row.get("record_index", "")).strip(),
            str(row.get("alphali_replay_index", "")).strip(),
        ]
    )


def step_uid(row):
    return "|".join(
        [
            event_uid(row),
            str(row.get("trackID", "")).strip(),
            str(row.get("stepID", "")).strip(),
            str(row.get("particle", "")).strip(),
        ]
    )


def input_prefix(input_csv):
    stem = Path(input_csv).stem
    suffix = "_alpha_li_steps"
    if stem.endswith(suffix):
        return stem[: -len(suffix)]
    return stem


def placement_subpath_from_input_layout(path_value):
    parts = Path(str(path_value).strip()).parts
    for index in range(len(parts) - 1):
        if parts[index] == "Input" and parts[index + 1] == "placements":
            tail = parts[index + 2 :]
            if tail:
                return Path(*tail)
            break
    return None


def path_for_record(path_value, project_root):
    raw = str(path_value).strip()
    if not raw:
        return raw

    project_root = Path(project_root).resolve()
    path = Path(raw)
    candidate = path if path.is_absolute() else (project_root / path)
    candidate = candidate.resolve(strict=False)

    try:
        return candidate.relative_to(project_root).as_posix()
    except ValueError:
        return candidate.as_posix()


def normalize_placement_file(path_value, project_root, ratio_tag=None):
    raw = str(path_value).strip()
    if not raw:
        return raw

    project_root = Path(project_root).resolve()
    placements_root = project_root / "Input" / "placements"
    raw_path = Path(raw)
    candidates = []

    direct_path = raw_path if raw_path.is_absolute() else (project_root / raw_path)
    if direct_path.exists():
        return path_for_record(direct_path, project_root)

    layout_tail = placement_subpath_from_input_layout(raw_path)
    if layout_tail is not None:
        candidates.append(placements_root / layout_tail)

    if raw_path.parent.name:
        candidates.append(placements_root / raw_path.parent.name / raw_path.name)

    if ratio_tag:
        candidates.append(placements_root / ratio_tag / raw_path.name)

    seen = set()
    for candidate in candidates:
        candidate = candidate.resolve(strict=False)
        key = candidate.as_posix()
        if key in seen:
            continue
        seen.add(key)
        if candidate.exists():
            return path_for_record(candidate, project_root)

    if candidates:
        return path_for_record(candidates[-1], project_root)

    return path_for_record(raw, project_root)


def visible_edep(edep_kev, step_len_um, birks_kb, line_number):
    if edep_kev < 0.0:
        raise SystemExit(f"Negative edep_keV at CSV line {line_number}: {edep_kev}")

    if birks_kb == 0.0 or edep_kev == 0.0:
        return edep_kev, 1.0

    if step_len_um <= 0.0:
        return 0.0, 0.0

    dedx = edep_kev / step_len_um
    visible = edep_kev / (1.0 + birks_kb * dedx)
    return visible, visible / edep_kev


def make_step_row(row, birks_kb, yield_zns, line_number):
    x_pre = to_float(row, "x_pre_um", line_number)
    y_pre = to_float(row, "y_pre_um", line_number)
    z_pre = to_float(row, "z_pre_um", line_number)
    x_post = to_float(row, "x_post_um", line_number)
    y_post = to_float(row, "y_post_um", line_number)
    z_post = to_float(row, "z_post_um", line_number)
    step_len_um = to_float(row, "step_len_um", line_number)
    edep_kev = to_float(row, "edep_keV", line_number)

    visible_kev, quenching_factor = visible_edep(
        edep_kev, step_len_um, birks_kb, line_number
    )
    n_photon_step = visible_kev / 1000.0 * yield_zns
    ratio_label = ratio_label_from_row(row)
    placement_id = placement_hash(row["placement_file"])

    return {
        "physical_event_uid": row.get("physical_event_uid", ""),
        "source_event_uid": event_uid(row),
        "source_step_uid": step_uid(row),
        "eventID": row["eventID"],
        "record_index": row.get("record_index", ""),
        "trackID": row["trackID"],
        "stepID": row["stepID"],
        "particle": row["particle"],
        "thickness_um": row["thickness_um"],
        "bn_wt": row["bn_wt"],
        "zns_wt": row["zns_wt"],
        "ratio_label": ratio_label,
        "depth_um": row["depth_um"],
        "capture_depth_um": row["depth_um"],
        "placement_file": row["placement_file"],
        "placement_hash": placement_id,
        "surface_mode": row.get("surface_mode", ""),
        "local_capture_x_um": row.get("local_capture_x_um", ""),
        "local_capture_y_um": row.get("local_capture_y_um", ""),
        "local_capture_z_um": row.get("local_capture_z_um", ""),
        "alphali_replay_index": row.get("alphali_replay_index", "0"),
        "alphali_replay_count": row.get("alphali_replay_count", "1"),
        "trajectory_weight": row.get("trajectory_weight", ""),
        "source_model": "trajectory_conditioned_alpha_li",
        "sampling_model": "uniformAlongStep",
        "phase_pre": row["phase_pre"],
        "phase_post": row["phase_post"],
        "x_pre_um": row["x_pre_um"],
        "y_pre_um": row["y_pre_um"],
        "z_pre_um": row["z_pre_um"],
        "x_post_um": row["x_post_um"],
        "y_post_um": row["y_post_um"],
        "z_post_um": row["z_post_um"],
        "x_mid_um": fmt(0.5 * (x_pre + x_post)),
        "y_mid_um": fmt(0.5 * (y_pre + y_post)),
        "z_mid_um": fmt(0.5 * (z_pre + z_post)),
        "step_len_um": row["step_len_um"],
        "edep_keV": row["edep_keV"],
        "visible_edep_keV": fmt(visible_kev),
        "quenching_factor_step": fmt(quenching_factor),
        "n_photon_step": fmt(n_photon_step),
    }, edep_kev, visible_kev


def ensure_event_summary(events, row):
    uid = event_uid(row)
    event = events.get(uid)
    if event is None:
        event = {
            "physical_event_uid": row.get("physical_event_uid", ""),
            "source_event_uid": uid,
            "eventID": row["eventID"],
            "record_index": row.get("record_index", ""),
            "thickness_um": row["thickness_um"],
            "bn_wt": row["bn_wt"],
            "zns_wt": row["zns_wt"],
            "ratio_label": ratio_label_from_row(row),
            "depth_um": row["depth_um"],
            "capture_depth_um": row["depth_um"],
            "capture_x_um": row.get("capture_x_um", ""),
            "capture_y_um": row.get("capture_y_um", ""),
            "local_capture_x_um": row.get("local_capture_x_um", ""),
            "local_capture_y_um": row.get("local_capture_y_um", ""),
            "local_capture_z_um": row.get("local_capture_z_um", ""),
            "surface_mode": row.get("surface_mode", ""),
            "placement_file": row["placement_file"],
            "placement_hash": placement_hash(row["placement_file"]),
            "alphali_replay_index": row.get("alphali_replay_index", "0"),
            "alphali_replay_count": row.get("alphali_replay_count", "1"),
            "trajectory_weight": row.get("trajectory_weight", ""),
            "n_total_steps": 0,
            "edep_total_keV": 0.0,
            "edep_ZnS_keV": 0.0,
            "visible_edep_ZnS_keV": 0.0,
            "n_zns_steps": 0,
        }
        events[uid] = event
    return event


def update_event_all_steps(events, row, edep_kev):
    event = ensure_event_summary(events, row)
    event["n_total_steps"] += 1
    event["edep_total_keV"] += edep_kev


def update_event_zns_summary(events, step_row, edep_kev, visible_kev):
    event = ensure_event_summary(events, step_row)
    event["edep_ZnS_keV"] += edep_kev
    event["visible_edep_ZnS_keV"] += visible_kev
    event["n_zns_steps"] += 1


def event_output_row(event, yield_zns):
    edep = event["edep_ZnS_keV"]
    visible = event["visible_edep_ZnS_keV"]
    quenching_factor = visible / edep if edep > 0.0 else 1.0
    n_photon0 = visible / 1000.0 * yield_zns

    return {
        "physical_event_uid": event["physical_event_uid"],
        "source_event_uid": event["source_event_uid"],
        "eventID": event["eventID"],
        "record_index": event["record_index"],
        "thickness_um": event["thickness_um"],
        "bn_wt": event["bn_wt"],
        "zns_wt": event["zns_wt"],
        "ratio_label": event["ratio_label"],
        "depth_um": event["depth_um"],
        "capture_depth_um": event["capture_depth_um"],
        "capture_x_um": event["capture_x_um"],
        "capture_y_um": event["capture_y_um"],
        "local_capture_x_um": event["local_capture_x_um"],
        "local_capture_y_um": event["local_capture_y_um"],
        "local_capture_z_um": event["local_capture_z_um"],
        "surface_mode": event["surface_mode"],
        "placement_file": event["placement_file"],
        "placement_hash": event["placement_hash"],
        "alphali_replay_index": event["alphali_replay_index"],
        "alphali_replay_count": event["alphali_replay_count"],
        "trajectory_weight": event["trajectory_weight"],
        "n_total_steps": str(event["n_total_steps"]),
        "edep_ZnS_keV": fmt(edep),
        "edep_total_keV": fmt(event["edep_total_keV"]),
        "visible_edep_ZnS_keV": fmt(visible),
        "quenching_factor_event": fmt(quenching_factor),
        "n_photon0": fmt(n_photon0),
        "n_zns_steps": str(event["n_zns_steps"]),
        "has_zns_edep": "1" if visible > 0.0 else "0",
    }


def process_csv(
    input_csv,
    step_output,
    event_output,
    yield_zns,
    birks_kb,
    project_root,
    ratio_tag=None,
):
    input_csv = Path(input_csv).resolve()
    step_output = Path(step_output).resolve()
    event_output = Path(event_output).resolve()

    if not input_csv.is_file():
        raise SystemExit(f"Input CSV not found: {input_csv}")

    if step_output == input_csv or event_output == input_csv:
        raise SystemExit("Refusing to overwrite the input alpha_li_steps.csv.")
    if step_output == event_output:
        raise SystemExit("Step output and event output must be different paths.")

    step_output.parent.mkdir(parents=True, exist_ok=True)
    event_output.parent.mkdir(parents=True, exist_ok=True)

    events = OrderedDict()
    input_rows = 0
    zns_rows = 0

    with open(input_csv, newline="", encoding="utf-8-sig") as src, open(
        step_output, "w", newline="", encoding="utf-8"
    ) as step_dst:
        reader = csv.DictReader(src)
        require_fields(reader.fieldnames, REQUIRED_INPUT_FIELDS)

        step_writer = csv.DictWriter(step_dst, fieldnames=STEP_OUTPUT_FIELDS)
        step_writer.writeheader()

        for input_rows, row in enumerate(reader, start=1):
            line_number = input_rows + 1
            row = dict(row)
            row["placement_file"] = normalize_placement_file(
                row.get("placement_file", ""),
                project_root=project_root,
                ratio_tag=ratio_tag,
            )
            edep_kev_all = to_float(row, "edep_keV", line_number)
            update_event_all_steps(events, row, edep_kev_all)

            # Only the pre-step phase assigns this step's energy deposition.
            if row["phase_pre"].strip() != "ZnS":
                continue

            step_row, edep_kev, visible_kev = make_step_row(
                row, birks_kb, yield_zns, line_number
            )
            step_writer.writerow(step_row)
            update_event_zns_summary(events, step_row, edep_kev, visible_kev)
            zns_rows += 1

    with open(event_output, "w", newline="", encoding="utf-8") as event_dst:
        event_writer = csv.DictWriter(event_dst, fieldnames=EVENT_OUTPUT_FIELDS)
        event_writer.writeheader()
        for event in events.values():
            event_writer.writerow(event_output_row(event, yield_zns))

    return {
        "input_csv": input_csv,
        "step_output": step_output,
        "event_output": event_output,
        "input_rows": input_rows,
        "zns_rows": zns_rows,
        "events": len(events),
        "placements": 1 if zns_rows > 0 else 0,
    }


def load_anchor_events(anchor_csv, project_root, ratio_tag=None):
    events = OrderedDict()
    with open(anchor_csv, newline="", encoding="utf-8-sig") as src:
        reader = csv.DictReader(src)
        require_fields(reader.fieldnames, REQUIRED_ANCHOR_FIELDS)
        for line_idx, row in enumerate(reader, start=2):
            row = dict(row)
            row["placement_file"] = normalize_placement_file(
                row.get("placement_file", ""),
                project_root=project_root,
                ratio_tag=ratio_tag,
            )
            uid = event_uid(row)
            row["source_event_uid"] = uid
            if uid in events:
                raise SystemExit(
                    f"Duplicate source_event_uid in capture anchors at line {line_idx}: {uid}"
                )
            events[uid] = {
                "physical_event_uid": row.get("physical_event_uid", ""),
                "source_event_uid": uid,
                "eventID": row["eventID"],
                "record_index": row.get("record_index", ""),
                "thickness_um": row["thickness_um"],
                "bn_wt": row["bn_wt"],
                "zns_wt": row["zns_wt"],
                "ratio_label": ratio_label_from_row(row),
                "depth_um": row["depth_um"],
                "capture_depth_um": row["depth_um"],
                "capture_x_um": row.get("capture_x_um", ""),
                "capture_y_um": row.get("capture_y_um", ""),
                "local_capture_x_um": row.get("local_capture_x_um", ""),
                "local_capture_y_um": row.get("local_capture_y_um", ""),
                "local_capture_z_um": row.get("local_capture_z_um", ""),
                "surface_mode": row.get("surface_mode", ""),
                "placement_file": row["placement_file"],
                "placement_hash": placement_hash(row["placement_file"]),
                "alphali_replay_index": row.get("alphali_replay_index", "0"),
                "alphali_replay_count": row.get("alphali_replay_count", "1"),
                "trajectory_weight": row.get("trajectory_weight", ""),
                "n_total_steps": 0,
                "edep_total_keV": 0.0,
                "edep_ZnS_keV": 0.0,
                "visible_edep_ZnS_keV": 0.0,
                "n_zns_steps": 0,
            }
    return events


def process_slim_pair(
    anchor_csv,
    track_csv,
    step_output,
    event_output,
    yield_zns,
    birks_kb,
    project_root,
    ratio_tag=None,
):
    anchor_csv = Path(anchor_csv).resolve()
    track_csv = Path(track_csv).resolve()
    step_output = Path(step_output).resolve()
    event_output = Path(event_output).resolve()

    if not anchor_csv.is_file():
        raise SystemExit(f"Capture anchor CSV not found: {anchor_csv}")
    if not track_csv.is_file():
        raise SystemExit(f"ZnS track CSV not found: {track_csv}")

    events = load_anchor_events(anchor_csv, project_root, ratio_tag=ratio_tag)
    input_rows = 0
    zns_rows = 0

    step_output.parent.mkdir(parents=True, exist_ok=True)
    event_output.parent.mkdir(parents=True, exist_ok=True)

    with open(track_csv, newline="", encoding="utf-8-sig") as src, open(
        step_output, "w", newline="", encoding="utf-8"
    ) as step_dst:
        reader = csv.DictReader(src)
        require_fields(reader.fieldnames, REQUIRED_TRACK_FIELDS)

        step_writer = csv.DictWriter(step_dst, fieldnames=STEP_OUTPUT_FIELDS)
        step_writer.writeheader()

        for input_rows, row in enumerate(reader, start=1):
            line_number = input_rows + 1
            row = dict(row)
            uid = event_uid(row)
            anchor = events.get(uid)
            if anchor is None:
                raise SystemExit(
                    f"Track row at line {line_number} references unknown source_event_uid: {uid}"
                )

            merged_row = {
                **anchor,
                **row,
                "source_event_uid": uid,
                "phase_pre": "ZnS",
                "placement_file": anchor["placement_file"],
                "bn_wt": anchor["bn_wt"],
                "zns_wt": anchor["zns_wt"],
                "thickness_um": anchor["thickness_um"],
                "depth_um": anchor["depth_um"],
                "local_capture_x_um": anchor["local_capture_x_um"],
                "local_capture_y_um": anchor["local_capture_y_um"],
                "local_capture_z_um": anchor["local_capture_z_um"],
                "surface_mode": anchor["surface_mode"],
                "capture_x_um": anchor["capture_x_um"],
                "capture_y_um": anchor["capture_y_um"],
                "alphali_replay_index": anchor["alphali_replay_index"],
                "alphali_replay_count": anchor["alphali_replay_count"],
            }

            edep_kev_all = to_float(merged_row, "edep_keV", line_number)
            update_event_all_steps(events, merged_row, edep_kev_all)

            step_row, edep_kev, visible_kev = make_step_row(
                merged_row, birks_kb, yield_zns, line_number
            )
            step_writer.writerow(step_row)
            update_event_zns_summary(events, step_row, edep_kev, visible_kev)
            zns_rows += 1

    with open(event_output, "w", newline="", encoding="utf-8") as event_dst:
        event_writer = csv.DictWriter(event_dst, fieldnames=EVENT_OUTPUT_FIELDS)
        event_writer.writeheader()
        for event in events.values():
            event_writer.writerow(event_output_row(event, yield_zns))

    placements = {
        safe_name(event["placement_file"])
        for event in events.values()
        if event.get("placement_file", "").strip()
    }
    return {
        "anchor_csv": anchor_csv,
        "track_csv": track_csv,
        "step_output": step_output,
        "event_output": event_output,
        "input_rows": input_rows,
        "zns_rows": zns_rows,
        "events": len(events),
        "placements": max(1, len(placements)) if events else 0,
    }


def open_split_sink(sinks, output_dir, prefix, placement_file):
    placement_key = placement_file.strip()
    sink = sinks.get(placement_key)
    if sink is not None:
        return sink

    placement_name = safe_name(placement_file)
    step_output = output_dir / f"{prefix}_{placement_name}_zns_step_sources.csv"
    event_output = output_dir / f"{prefix}_{placement_name}_event_light_sources.csv"

    if step_output == event_output:
        raise SystemExit("Internal error: split step/event output paths are identical.")

    step_file = open(step_output, "w", newline="", encoding="utf-8")
    writer = csv.DictWriter(step_file, fieldnames=STEP_OUTPUT_FIELDS)
    writer.writeheader()

    sink = {
        "placement_file": placement_key,
        "step_output": step_output,
        "event_output": event_output,
        "step_file": step_file,
        "writer": writer,
        "events": OrderedDict(),
        "zns_rows": 0,
    }
    sinks[placement_key] = sink
    return sink


def close_split_sinks(sinks, yield_zns):
    for sink in sinks.values():
        sink["step_file"].flush()
        sink["step_file"].close()

        with open(sink["event_output"], "w", newline="", encoding="utf-8") as event_dst:
            event_writer = csv.DictWriter(event_dst, fieldnames=EVENT_OUTPUT_FIELDS)
            event_writer.writeheader()
            for event in sink["events"].values():
                event_writer.writerow(event_output_row(event, yield_zns))


def process_csv_split_by_placement(
    input_csv,
    output_dir,
    prefix,
    yield_zns,
    birks_kb,
    project_root,
    ratio_tag=None,
):
    input_csv = Path(input_csv).resolve()
    output_dir = Path(output_dir).resolve()

    if not input_csv.is_file():
        raise SystemExit(f"Input CSV not found: {input_csv}")

    split_output_dir = output_dir / "by_placement" / prefix
    split_output_dir.mkdir(parents=True, exist_ok=True)

    sinks = OrderedDict()
    input_rows = 0
    zns_rows = 0

    try:
        with open(input_csv, newline="", encoding="utf-8-sig") as src:
            reader = csv.DictReader(src)
            require_fields(reader.fieldnames, REQUIRED_INPUT_FIELDS)

            for input_rows, row in enumerate(reader, start=1):
                line_number = input_rows + 1
                row = dict(row)
                row["placement_file"] = normalize_placement_file(
                    row.get("placement_file", ""),
                    project_root=project_root,
                    ratio_tag=ratio_tag,
                )
                edep_kev_all = to_float(row, "edep_keV", line_number)
                sink = open_split_sink(
                    sinks,
                    split_output_dir,
                    prefix,
                    row["placement_file"],
                )
                update_event_all_steps(sink["events"], row, edep_kev_all)

                # Only the pre-step phase assigns this step's energy deposition.
                if row["phase_pre"].strip() != "ZnS":
                    continue

                step_row, edep_kev, visible_kev = make_step_row(
                    row, birks_kb, yield_zns, line_number
                )

                sink["writer"].writerow(step_row)
                update_event_zns_summary(sink["events"], step_row, edep_kev, visible_kev)
                sink["zns_rows"] += 1
                zns_rows += 1
    finally:
        close_split_sinks(sinks, yield_zns)

    return {
        "input_csv": input_csv,
        "output_dir": split_output_dir,
        "input_rows": input_rows,
        "zns_rows": zns_rows,
        "events": sum(len(sink["events"]) for sink in sinks.values()),
        "placements": len(sinks),
        "outputs": [
            (sink["placement_file"], sink["step_output"], sink["event_output"], sink["zns_rows"], len(sink["events"]))
            for sink in sinks.values()
        ],
    }


def print_one_summary(stats):
    print(f"Input rows:       {stats['input_rows']}")
    print(f"ZnS source steps: {stats['zns_rows']}")
    print(f"Capture events:   {stats['events']}")
    print(f"Step output:      {stats['step_output']}")
    print(f"Event output:     {stats['event_output']}")
    print("ZnS rule:         phase_pre == 'ZnS' only; phase_post is metadata only.")
    print("Source model:     default Stage C optical photons sample uniformly on pre->post segments.")


def print_split_summary(stats):
    print(f"Input rows:       {stats['input_rows']}")
    print(f"ZnS source steps: {stats['zns_rows']}")
    print(f"Capture events:   {stats['events']}")
    print(f"Placements:       {stats['placements']}")
    print(f"Output dir:       {stats['output_dir']}")
    for placement_file, step_output, event_output, zns_rows, events in stats["outputs"]:
        print(
            f">>> {Path(placement_file).name}: ZnS_steps={zns_rows} "
            f"events={events} -> {step_output.name}, {event_output.name}"
        )
    print("ZnS rule:         phase_pre == 'ZnS' only; phase_post is metadata only.")
    print("Source model:     default Stage C optical photons sample uniformly on pre->post segments.")


def single_file_mode(args):
    project_root = (
        Path(args.project_root).resolve()
        if args.project_root
        else Path(__file__).resolve().parent
    )
    input_csv = Path(args.input_csv).resolve()
    output_dir = Path(args.output_dir).resolve() if args.output_dir else input_csv.parent
    if args.split_by_placement:
        prefix = args.output_prefix or input_prefix(input_csv)
        stats = process_csv_split_by_placement(
            input_csv,
            output_dir,
            prefix,
            yield_zns=args.yield_zns,
            birks_kb=args.birks_kb,
            project_root=project_root,
        )
        print_split_summary(stats)
        return

    step_output = resolve_output_path(args.step_output, output_dir).resolve()
    event_output = resolve_output_path(args.event_output, output_dir).resolve()
    stats = process_csv(
        input_csv,
        step_output,
        event_output,
        yield_zns=args.yield_zns,
        birks_kb=args.birks_kb,
        project_root=project_root,
    )
    print_one_summary(stats)


def stageb_batch_files(stageb_ratio_dir, start_um, end_um):
    matches = []
    slim_anchors = {}
    slim_tracks = {}

    for path in stageb_ratio_dir.glob("*_capture_anchors.csv"):
        match = CAPTURE_ANCHORS_RE.fullmatch(path.name)
        if match is None:
            continue
        thickness_label = match.group("thickness")
        thickness_um = float(thickness_label)
        if start_um <= thickness_um <= end_um:
            slim_anchors[thickness_label] = (thickness_um, path)

    for path in stageb_ratio_dir.glob("*_zns_track_steps.csv"):
        match = ZNS_TRACK_STEPS_RE.fullmatch(path.name)
        if match is None:
            continue
        thickness_label = match.group("thickness")
        thickness_um = float(thickness_label)
        if start_um <= thickness_um <= end_um:
            slim_tracks[thickness_label] = (thickness_um, path)

    for thickness_label, (thickness_um, anchor_path) in slim_anchors.items():
        track_entry = slim_tracks.get(thickness_label)
        if track_entry is None:
            continue
        matches.append((thickness_um, thickness_label, ("slim", anchor_path, track_entry[1])))

    for path in stageb_ratio_dir.glob("*_alpha_li_steps.csv"):
        match = ALPHA_LI_STEPS_RE.fullmatch(path.name)
        if match is None:
            continue
        thickness_label = match.group("thickness")
        thickness_um = float(thickness_label)
        if start_um <= thickness_um <= end_um:
            if thickness_label in slim_anchors and thickness_label in slim_tracks:
                continue
            matches.append((thickness_um, thickness_label, ("full", path)))
    return sorted(matches, key=lambda item: (item[0], item[1]))


def batch_mode(args):
    project_root = Path(args.project_root).resolve() if args.project_root else Path(__file__).resolve().parent
    ratio_tag = args.ratio_tag or prompt_required(
        "请输入配比文件夹名，例如 1-1, 1-1.5, 1-2: "
    )
    thickness_range = args.thickness_range or prompt_required(
        "请输入厚度范围，例如 30-200: "
    )
    start_um, end_um = parse_thickness_range(thickness_range)

    stageb_ratio_dir = project_root / "Output" / "stageB" / ratio_tag
    if not stageb_ratio_dir.is_dir():
        raise SystemExit(f"Input ratio directory not found: {stageb_ratio_dir}")

    output_root = (
        Path(args.batch_output_root).resolve()
        if args.batch_output_root
        else project_root / "Input" / "alpha_li_steps"
    )
    output_ratio_dir = output_root / ratio_tag
    files = stageb_batch_files(stageb_ratio_dir, start_um, end_um)
    if not files:
        raise SystemExit(
            f"No *_alpha_li_steps.csv files found in {stageb_ratio_dir} "
            f"for thickness range {start_um:g}-{end_um:g} um."
        )

    output_ratio_dir.mkdir(parents=True, exist_ok=True)

    print(f"Ratio:           {ratio_tag}")
    print(f"Thickness range: {start_um:g}-{end_um:g} um")
    print(f"Input dir:       {stageb_ratio_dir}")
    print(f"Output dir:      {output_ratio_dir}")
    print(f"Files to process:{len(files)}")
    if args.split_by_placement:
        print("Split mode:      one derived source pair per placement_file")

    totals = {"input_rows": 0, "zns_rows": 0, "events": 0}
    total_placements = 0
    for thickness_um, thickness_label, stageb_input in files:
        if args.split_by_placement:
            raise SystemExit(
                "--split-by-placement is currently supported only for full "
                "Stage B alpha_li_steps.csv inputs."
            )
        else:
            step_output = output_ratio_dir / f"{thickness_label}_zns_step_sources.csv"
            event_output = output_ratio_dir / f"{thickness_label}_event_light_sources.csv"

            if stageb_input[0] == "slim":
                _, anchor_csv, track_csv = stageb_input
                stats = process_slim_pair(
                    anchor_csv,
                    track_csv,
                    step_output,
                    event_output,
                    yield_zns=args.yield_zns,
                    birks_kb=args.birks_kb,
                    project_root=project_root,
                    ratio_tag=ratio_tag,
                )
                input_desc = f"{anchor_csv.name} + {track_csv.name}"
            else:
                _, input_csv = stageb_input
                stats = process_csv(
                    input_csv,
                    step_output,
                    event_output,
                    yield_zns=args.yield_zns,
                    birks_kb=args.birks_kb,
                    project_root=project_root,
                    ratio_tag=ratio_tag,
                )
                input_desc = input_csv.name
        totals["input_rows"] += stats["input_rows"]
        totals["zns_rows"] += stats["zns_rows"]
        totals["events"] += stats["events"]
        total_placements += stats["placements"]
        if args.split_by_placement:
            print(
                f">>> {thickness_um:g} um: rows={stats['input_rows']} "
                f"ZnS_steps={stats['zns_rows']} events={stats['events']} "
                f"placements={stats['placements']}"
            )
            for _, step_output, event_output, zns_rows, events in stats["outputs"]:
                print(
                    f"    {step_output.name}: ZnS_steps={zns_rows} "
                    f"events={events}; event file={event_output.name}"
                )
        else:
            print(
                f">>> {thickness_um:g} um: rows={stats['input_rows']} "
                f"ZnS_steps={stats['zns_rows']} events={stats['events']} "
                f"input={input_desc} -> {step_output.name}, {event_output.name}"
            )

    print()
    print(f"Processed files:  {len(files)}")
    print(f"Total input rows: {totals['input_rows']}")
    print(f"Total ZnS steps:  {totals['zns_rows']}")
    print(f"Total events:     {totals['events']}")
    if args.split_by_placement:
        print(f"Total placements: {total_placements}")
    print("ZnS rule:         phase_pre == 'ZnS' only; phase_post is metadata only.")
    print("Source model:     default Stage C optical photons sample uniformly on pre->post segments.")


def main():
    args = parse_args()

    if args.yield_zns < 0.0:
        raise SystemExit("--yield-zns must be >= 0")
    if args.birks_kb < 0.0:
        raise SystemExit("--birks-kb must be >= 0")

    if args.input_csv:
        single_file_mode(args)
    else:
        batch_mode(args)


if __name__ == "__main__":
    main()
