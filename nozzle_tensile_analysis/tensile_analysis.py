#!/usr/bin/env python3
"""
Tensile test analysis: stress vs. strain across 3D-printed nozzle sizes and
print modes, combining multiple experiment campaigns.

Two folder-naming conventions are supported, both encoding the same physics.

NEW convention (e.g. "Raw Data/20AUG2026"):

    T105N06-1GFR085/
        Test Run 74 8-21-2026 4 21 37 PM/
            DAQ- Crosshead, ... - (Timed).csv

    T105  -> print temperature, 105 C
    N06   -> 0.6 mm nozzle (two digits, tenths implied)
    -1G   -> inverted printer while printing (+1G = upright)
    FR085 -> flow ratio 0.85

OLD convention (e.g. "08_Tensile_Tests"):

    T105N0_6_1G/    -> 0.6 mm nozzle, +1G upright
    T105N0_6(1G)/   -> 0.6 mm nozzle, -1G inverted
    T105N0_6CLS/    -> 0.6 mm nozzle, printed while rotating (CLS)

    N0_6  -> 0.6 mm nozzle (underscore is the decimal point)
    _1G   -> +1G upright
    (1G)  -> -1G inverted
    CLS   -> rotating print mode (RPM not encoded in the folder name)

The old campaign does not encode a flow ratio; those runs get flow_ratio = NaN.

Each DAQ csv has 4 metadata lines, a blank line, a header row, and a units
row before the numeric data starts:

    "Crosshead ","Load ","Time ","Extensometer "
    "mm","N","sec","mm"

Run this directly:

    python3 tensile_analysis.py
    python3 tensile_analysis.py --data-root "/path/to/campaign_a" \
                               --data-root "/path/to/campaign_b"

All the "important variables" -- paths, specimen geometry, cleaning
thresholds, fit windows, QC limits -- live in one place: the CONFIG object
below. Nothing else in the script hardcodes a magic number; edit CONFIG and
rerun.
"""

from __future__ import annotations

import argparse
import dataclasses
import glob
import io
import re
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _running_inside_ipython_kernel() -> bool:
    """True in Spyder/Jupyter (they pre-configure their own matplotlib
    backend via %matplotlib magic before this script's code runs), False for
    a plain `python tensile_analysis.py` invocation."""
    try:
        return get_ipython() is not None  # noqa: F821 -- injected by IPython
    except NameError:
        return False


if not _running_inside_ipython_kernel():
    # Plain script execution: don't let matplotlib auto-select its native
    # backend. On very new Python releases (e.g. 3.14), matplotlib's native
    # macOS backend (a compiled C extension) can crash outright --
    # `SystemError: NULL object passed to Py_BuildValue` then a segfault --
    # before that C code has caught up with the new CPython C-API. TkAgg
    # doesn't touch Apple's Cocoa APIs and has been stable for decades.
    try:
        matplotlib.use("TkAgg")
    except Exception:
        matplotlib.use("Agg")


# Print-mode labels. Used as dict keys, plot legend text and dataframe values,
# so they are defined once here rather than repeated as string literals.
MODE_UPRIGHT = "+1G (upright)"
MODE_INVERTED = "-1G (inverted)"
MODE_ROTATING = "CLS (rotating)"


# ============================================================================
# CONFIG -- the "local environment" for this analysis.
# Every knob that affects parsing, cleaning, or curve-fitting lives here.
# ============================================================================
@dataclasses.dataclass
class Config:
    # --- where the raw data lives / where results get written -------------
    # Every root is scanned recursively; each contributes a "dataset" label
    # (the root folder's own name) so campaigns stay distinguishable.
    data_roots: tuple = (
        Path(
            "/Users/siddharthchandel/Claude/Projects/Research/07_Experiments/"
            "Raw Data/20AUG2026"
        ),
        Path("/Users/siddharthchandel/Claude/Projects/Research/08_Tensile_Tests"),
    )

    # Anchored to this script's own location (not the current working
    # directory) so it always lands in the same place regardless of how or
    # from where you run the script.
    output_dir: Path = Path(__file__).resolve().parent / "analysis_output"

    # --- specimen geometry ---------------------------------------------------
    # Custom dogbone, assumed constant across every nozzle size, print mode
    # AND campaign. If an older campaign used different dimensions, split this
    # into a per-dataset lookup -- a wrong area silently rescales that whole
    # campaign's stress axis.
    specimen_width_mm: float = 20.0
    specimen_thickness_mm: float = 4.0

    @property
    def area_mm2(self) -> float:
        """Gauge-section area (mm^2). Load[N]/Area = Stress[MPa]."""
        return self.specimen_width_mm * self.specimen_thickness_mm

    # --- CSV parsing -----------------------------------------------------
    header_marker: str = '"Crosshead'
    expected_columns: tuple = ("Crosshead", "Load", "Time", "Extensometer")
    csv_glob: str = "**/*.csv"
    # macOS zip archives carry a parallel __MACOSX/ tree of resource forks
    # that look like csvs but are not.
    skip_path_substrings: tuple = ("__MACOSX", "/._", "\\._")

    # --- data cleaning -----------------------------------------------------
    # Pre-pull load-cell noise: many runs record small (and occasionally very
    # large) offsets before the crosshead starts moving. Zero it out using the
    # mean Load measured while Crosshead is still below this threshold (mm).
    crosshead_zero_threshold_mm: float = 0.01
    min_baseline_samples: int = 5

    # Extensometer occasionally faults, logging one physically-impossible
    # instantaneous jump (tens of mm in a single 0.01 s sample -- normal
    # motion is ~0.01-0.05 mm/sample). Trim that sample and everything after.
    max_extensometer_step_mm: float = 1.0

    # --- strain reference ---------------------------------------------------
    # Gauge length L0 is not a fixed nominal value in this rig -- Extensometer
    # reports an absolute position starting a little differently per specimen
    # (~23-25.5 mm across both campaigns). L0 = each run's first valid reading.
    strain_from: str = "Extensometer"

    # --- elastic-modulus fit window (engineering strain) --------------------
    modulus_strain_lo: float = 0.0005
    modulus_strain_hi: float = 0.01
    modulus_min_points: int = 20

    # --- quality control / flagging -----------------------------------------
    # Runs are FLAGGED, never silently dropped. Review flagged runs before
    # treating them as equivalent to their replicates.
    #
    # A baseline offset this large means the load cell was badly zeroed, so
    # every load reading in that run is suspect (not just the first samples).
    max_abs_load_offset_n: float = 50.0
    # Fraction of rows lost to sensor dropout before the run is untrustworthy.
    max_nan_fraction: float = 0.05
    # A run whose recorded strain range exceeds this multiple of its *other*
    # replicates (leave-one-out, so one outlier can't inflate its own
    # reference) was probably stopped much later than they were.
    outlier_strain_ratio: float = 1.75
    # Run labels to exclude entirely, e.g. ("Test Run 52",). Empty by default
    # so nothing disappears without an explicit decision recorded here.
    exclude_runs: tuple = ()

    # --- plotting -------------------------------------------------------------
    dpi: int = 150
    max_curve_points: int = 4000  # decimate long curves for lighter plot files
    show_plots: bool = True


CFG = Config()


# ============================================================================
# Condition parsing -- both naming conventions
# ============================================================================
# T105N06-1GFR085 : nozzle as two digits (tenths implied), explicit +/- and
# flow ratio.
NEW_CONDITION_RE = re.compile(
    r"^T(?P<temp_c>\d+)"
    r"N(?P<nozzle_x10>\d{2})"
    r"(?P<gravity>[+-])1G"
    r"FR(?P<flow_ratio_x100>\d+)$"
)

# T105N0_6_1G / T105N0_6(1G) / T105N0_6CLS : nozzle as digit_digit (the
# underscore is the decimal point), print mode carried by a trailing suffix.
OLD_CONDITION_RE = re.compile(
    r"^T(?P<temp_c>\d+)"
    r"N(?P<nozzle_int>\d)_(?P<nozzle_dec>\d)"
    r"(?P<suffix>.*)$"
)


def _mode_from_old_suffix(suffix: str) -> str:
    """Map an old-convention trailing suffix to a print mode.

    Order matters: "(1G)" must be tested before the bare "1G", since the
    inverted marker contains the upright one as a substring.
    """
    s = suffix.upper()
    if "CLS" in s:
        return MODE_ROTATING
    if "(1G)" in s:
        return MODE_INVERTED
    if "1G" in s:
        return MODE_UPRIGHT
    raise ValueError(f"Unrecognized print-mode suffix: {suffix!r}")


def parse_condition(path: Path) -> dict:
    """Pull temperature / nozzle / print mode / flow ratio out of whichever
    path component matches a known naming convention."""
    for part in path.parts:
        m = NEW_CONDITION_RE.match(part)
        if m:
            g = m.groupdict()
            return {
                "condition": part,
                "temp_c": int(g["temp_c"]),
                "nozzle_mm": int(g["nozzle_x10"]) / 10.0,
                "print_mode": MODE_UPRIGHT if g["gravity"] == "+" else MODE_INVERTED,
                "flow_ratio": int(g["flow_ratio_x100"]) / 100.0,
            }

        m = OLD_CONDITION_RE.match(part)
        if m:
            g = m.groupdict()
            return {
                "condition": part,
                "temp_c": int(g["temp_c"]),
                "nozzle_mm": float(f"{g['nozzle_int']}.{g['nozzle_dec']}"),
                "print_mode": _mode_from_old_suffix(g["suffix"]),
                # Old campaign does not encode a flow ratio in the folder name.
                "flow_ratio": np.nan,
            }

    raise ValueError(f"Could not parse condition nomenclature from path: {path}")


def parse_run_label(path: Path) -> str:
    for part in path.parts:
        if part.lower().startswith("test run"):
            m = re.match(r"(Test Run \d+)", part)
            return m.group(1) if m else part
    return path.parent.name


# ============================================================================
# Loading
# ============================================================================
def load_daq_csv(path: Path) -> pd.DataFrame:
    """Tolerant parser for the MTS DAQ export: skips the metadata preamble,
    strips whitespace/quotes from headers, coerces everything to numeric."""
    with open(path, encoding="utf-8-sig", errors="replace") as fh:
        lines = fh.readlines()

    header_idx = next(
        (i for i, l in enumerate(lines) if l.strip().startswith(CFG.header_marker)),
        None,
    )
    if header_idx is None:
        raise ValueError(f"Header row not found in {path}")

    header = [c.strip().strip('"').strip() for c in lines[header_idx].strip().split(",")]
    data_str = "".join(lines[header_idx + 2 :])  # skip header row + units row
    df = pd.read_csv(io.StringIO(data_str), header=None, names=header)

    missing = [c for c in CFG.expected_columns if c not in df.columns]
    if missing:
        raise ValueError(f"{path} missing expected columns {missing}, got {list(df.columns)}")

    for c in CFG.expected_columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    return df[list(CFG.expected_columns)]


# ============================================================================
# Cleaning / polishing
# ============================================================================
@dataclasses.dataclass
class CleanResult:
    df: pd.DataFrame
    n_raw: int
    n_nan_dropped: int
    n_glitch_trimmed: int
    load_offset_n: float
    n_baseline_samples: int
    gauge_length_mm: float


def clean_run(df_raw: pd.DataFrame) -> CleanResult:
    n_raw = len(df_raw)

    # 1) Baseline-correct Load using samples taken before the crosshead
    #    actually starts moving.
    pre_motion = df_raw["Crosshead"] < CFG.crosshead_zero_threshold_mm
    n_baseline = int(pre_motion.sum())
    load_offset = (
        float(df_raw.loc[pre_motion, "Load"].mean())
        if n_baseline >= CFG.min_baseline_samples
        else 0.0
    )
    df = df_raw.copy()
    df["Load"] = df["Load"] - load_offset

    # 2) Drop rows with any NaN (sensor dropout).
    n_before = len(df)
    df = df.dropna().reset_index(drop=True)
    n_nan_dropped = n_before - len(df)

    # 3) Trim trailing sensor-glitch garbage.
    n_before_glitch = len(df)
    ext_step = df["Extensometer"].diff().abs()
    glitch_idx = ext_step[ext_step > CFG.max_extensometer_step_mm].index
    if len(glitch_idx):
        df = df.iloc[: glitch_idx[0]].reset_index(drop=True)
    n_glitch_trimmed = n_before_glitch - len(df)

    if df.empty:
        raise ValueError("No usable rows left after cleaning")

    gauge_length = float(df[CFG.strain_from].iloc[0])

    return CleanResult(
        df=df,
        n_raw=n_raw,
        n_nan_dropped=n_nan_dropped,
        n_glitch_trimmed=n_glitch_trimmed,
        load_offset_n=load_offset,
        n_baseline_samples=n_baseline,
        gauge_length_mm=gauge_length,
    )


def compute_stress_strain(clean: CleanResult) -> pd.DataFrame:
    df = clean.df.copy()
    df["Stress_MPa"] = df["Load"] / CFG.area_mm2
    ref = df[CFG.strain_from]
    df["Strain"] = (ref - clean.gauge_length_mm) / clean.gauge_length_mm
    return df


# ============================================================================
# Per-run metrics
# ============================================================================
def modulus_fit(df: pd.DataFrame) -> float | None:
    seg = df[
        (df["Strain"] >= CFG.modulus_strain_lo) & (df["Strain"] <= CFG.modulus_strain_hi)
    ]
    if len(seg) < CFG.modulus_min_points:
        return None
    slope, _intercept = np.polyfit(seg["Strain"], seg["Stress_MPa"], 1)
    return float(slope)


def summarize_run(df: pd.DataFrame, clean: CleanResult) -> dict:
    peak_i = int(df["Stress_MPa"].values.argmax())
    peak_stress = float(df["Stress_MPa"].iloc[peak_i])
    final_stress = float(df["Stress_MPa"].iloc[-1])
    return {
        "n_points": len(df),
        "n_raw_points": clean.n_raw,
        "n_nan_dropped": clean.n_nan_dropped,
        "n_glitch_trimmed": clean.n_glitch_trimmed,
        "nan_fraction": clean.n_nan_dropped / clean.n_raw if clean.n_raw else np.nan,
        "load_offset_n": clean.load_offset_n,
        "gauge_length_mm": clean.gauge_length_mm,
        "uts_mpa": peak_stress,
        "strain_at_uts": float(df["Strain"].iloc[peak_i]),
        "modulus_mpa": modulus_fit(df),
        "final_recorded_strain": float(df["Strain"].iloc[-1]),
        "final_recorded_stress_mpa": final_stress,
        "final_stress_frac_of_peak": final_stress / peak_stress if peak_stress else np.nan,
    }


# ============================================================================
# Driver: discover, load, clean, summarize every run across every root
# ============================================================================
def discover_runs(data_root: Path) -> list[Path]:
    files = sorted(Path(p) for p in glob.glob(str(data_root / CFG.csv_glob), recursive=True))
    files = [f for f in files if not any(s in str(f) for s in CFG.skip_path_substrings)]
    if not files:
        raise FileNotFoundError(
            f"No CSV files found under {data_root}. Check --data-root / CONFIG.data_roots."
        )
    return files


def build_dataset(data_roots):
    curves, summary_rows, quality_notes = [], [], []

    for root in data_roots:
        root = Path(root)
        dataset = root.name
        for f in discover_runs(root):
            cond = parse_condition(f.relative_to(root))
            run_label = parse_run_label(f)

            if run_label in CFG.exclude_runs:
                quality_notes.append(
                    f"[{dataset} / {cond['condition']} / {run_label}] EXCLUDED via "
                    f"CONFIG.exclude_runs."
                )
                continue

            clean = clean_run(load_daq_csv(f))
            df = compute_stress_strain(clean)
            summ = summarize_run(df, clean)

            summary_rows.append(
                {"dataset": dataset, **cond, "run": run_label, "file": str(f), **summ}
            )

            if (
                clean.n_baseline_samples >= CFG.min_baseline_samples
                and abs(clean.load_offset_n) > 0.5
            ):
                quality_notes.append(
                    f"[{dataset} / {cond['condition']} / {run_label}] Load zero-offset "
                    f"corrected by {clean.load_offset_n:+.2f} N "
                    f"(from {clean.n_baseline_samples} pre-motion samples)."
                )
            if clean.n_nan_dropped:
                quality_notes.append(
                    f"[{dataset} / {cond['condition']} / {run_label}] Dropped "
                    f"{clean.n_nan_dropped} row(s) with NaN values (sensor dropout)."
                )
            if clean.n_glitch_trimmed:
                quality_notes.append(
                    f"[{dataset} / {cond['condition']} / {run_label}] Trimmed "
                    f"{clean.n_glitch_trimmed} trailing row(s) after a "
                    f">{CFG.max_extensometer_step_mm}mm single-sample Extensometer "
                    f"jump (sensor glitch)."
                )

            curve = df[["Strain", "Stress_MPa"]].copy()
            if len(curve) > CFG.max_curve_points:
                curve = curve.iloc[:: max(1, len(curve) // CFG.max_curve_points)]
            curve["dataset"] = dataset
            curve["condition"] = cond["condition"]
            curve["nozzle_mm"] = cond["nozzle_mm"]
            curve["print_mode"] = cond["print_mode"]
            curve["run"] = run_label
            curves.append(curve)

    summary = pd.DataFrame(summary_rows)
    if summary.empty:
        raise RuntimeError("No runs survived loading -- nothing to analyse.")

    # --- QC flags (flag, never silently drop) ---
    summary["flag_load_offset"] = summary["load_offset_n"].abs() > CFG.max_abs_load_offset_n
    summary["flag_sensor_dropout"] = summary["nan_fraction"] > CFG.max_nan_fraction

    # Over-extension, leave-one-out within each dataset+condition group.
    summary["flag_over_extended"] = False
    summary["extension_reference"] = np.nan
    group_keys = ["dataset", "condition"]
    for _, idx in summary.groupby(group_keys).groups.items():
        idx = list(idx)
        if len(idx) < 2:
            continue
        for i in idx:
            others = [j for j in idx if j != i]
            ref = summary.loc[others, "final_recorded_strain"].max()
            summary.loc[i, "extension_reference"] = ref
            if ref > 0 and summary.loc[i, "final_recorded_strain"] > CFG.outlier_strain_ratio * ref:
                summary.loc[i, "flag_over_extended"] = True

    summary["flagged"] = (
        summary["flag_load_offset"]
        | summary["flag_sensor_dropout"]
        | summary["flag_over_extended"]
    )

    for _, r in summary[summary["flagged"]].iterrows():
        reasons = []
        if r["flag_load_offset"]:
            reasons.append(
                f"load-cell baseline offset {r['load_offset_n']:+.1f} N exceeds "
                f"+/-{CFG.max_abs_load_offset_n} N (every load reading in this run "
                f"is suspect, not just the first samples)"
            )
        if r["flag_sensor_dropout"]:
            reasons.append(
                f"{r['nan_fraction']*100:.1f}% of rows lost to sensor dropout "
                f"(limit {CFG.max_nan_fraction*100:.0f}%)"
            )
        if r["flag_over_extended"]:
            reasons.append(
                f"recorded strain {r['final_recorded_strain']:.2f} is "
                f">{CFG.outlier_strain_ratio}x its other replicate(s) "
                f"({r['extension_reference']:.2f})"
            )
        quality_notes.append(
            f"[{r['dataset']} / {r['condition']} / {r['run']}] FLAGGED: "
            + "; ".join(reasons)
            + ". Included in plots/summary -- review manually."
        )

    return pd.concat(curves, ignore_index=True), summary, quality_notes


# ============================================================================
# Plotting
# ============================================================================
MODE_COLOR = {
    MODE_UPRIGHT: "#2b6cb0",
    MODE_INVERTED: "#c0392b",
    MODE_ROTATING: "#2e8b57",
}
# Datasets are distinguished by line style so print mode keeps its colour.
DATASET_LINESTYLES = ["-", "--", ":", "-."]


def dataset_linestyles(all_curves: pd.DataFrame) -> dict:
    names = sorted(all_curves["dataset"].unique())
    return {n: DATASET_LINESTYLES[i % len(DATASET_LINESTYLES)] for i, n in enumerate(names)}


def _legend_handles(all_curves: pd.DataFrame, ls_map: dict):
    handles = [
        plt.Line2D([], [], color=c, ls="-", label=m)
        for m, c in MODE_COLOR.items()
        if m in set(all_curves["print_mode"])
    ]
    handles += [
        plt.Line2D([], [], color="0.35", ls=ls, label=f"dataset: {name}")
        for name, ls in ls_map.items()
    ]
    return handles


def plot_per_nozzle(all_curves: pd.DataFrame, summary: pd.DataFrame, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    ls_map = dataset_linestyles(all_curves)

    for nz in sorted(all_curves["nozzle_mm"].unique()):
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        sub = all_curves[all_curves["nozzle_mm"] == nz]
        for (dataset, mode, run), g in sub.groupby(["dataset", "print_mode", "run"]):
            flagged = summary.loc[
                (summary["run"] == run)
                & (summary["dataset"] == dataset)
                & (summary["nozzle_mm"] == nz),
                "flagged",
            ]
            is_flagged = bool(flagged.iloc[0]) if len(flagged) else False
            ax.plot(
                g["Strain"],
                g["Stress_MPa"],
                color=MODE_COLOR.get(mode, "0.4"),
                ls=ls_map[dataset],
                alpha=0.45 if is_flagged else 0.9,
                lw=1.0 if is_flagged else 1.6,
            )

        ax.set_xlabel("Engineering strain [mm/mm]")
        ax.set_ylabel("Engineering stress [MPa]")
        ax.set_title(f"Stress vs. Strain -- Nozzle {nz:.1f} mm (T105)")
        ax.grid(alpha=0.3)
        ax.legend(handles=_legend_handles(sub, ls_map), fontsize=8, loc="best")
        fig.tight_layout()
        fig.savefig(outdir / f"stress_strain_nozzle_{nz:.1f}mm.png", dpi=CFG.dpi)
        if CFG.show_plots:
            plt.show()
        plt.close(fig)


def plot_grid(all_curves: pd.DataFrame, summary: pd.DataFrame, outdir: Path):
    ls_map = dataset_linestyles(all_curves)
    nozzles = sorted(all_curves["nozzle_mm"].unique())
    n = len(nozzles)
    ncols = 2
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.5 * ncols, 4.5 * nrows), squeeze=False)

    for ax, nz in zip(axes.flat, nozzles):
        sub = all_curves[all_curves["nozzle_mm"] == nz]
        for (dataset, mode, _run), g in sub.groupby(["dataset", "print_mode", "run"]):
            ax.plot(
                g["Strain"],
                g["Stress_MPa"],
                color=MODE_COLOR.get(mode, "0.4"),
                ls=ls_map[dataset],
                alpha=0.85,
                lw=1.3,
            )
        ax.set_title(f"Nozzle {nz:.1f} mm")
        ax.set_xlabel("Strain [mm/mm]")
        ax.set_ylabel("Stress [MPa]")
        ax.grid(alpha=0.3)

    for ax in axes.flat[n:]:
        ax.axis("off")

    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.suptitle("Stress-Strain by Nozzle Size, Print Mode and Campaign (T105)", y=0.985)
    fig.legend(
        handles=_legend_handles(all_curves, ls_map),
        loc="upper center",
        ncol=3,
        fontsize=9,
        bbox_to_anchor=(0.5, 0.945),
    )
    fig.savefig(outdir / "stress_strain_grid_all.png", dpi=CFG.dpi, bbox_inches="tight")
    if CFG.show_plots:
        plt.show()
    plt.close(fig)


def plot_summary_bars(summary: pd.DataFrame, outdir: Path):
    """UTS and modulus by nozzle size, grouped by print mode + campaign."""
    summary = summary.copy()
    summary["group"] = summary["print_mode"] + "  [" + summary["dataset"] + "]"

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    nozzles = sorted(summary["nozzle_mm"].unique())
    groups = sorted(summary["group"].unique())
    x = np.arange(len(nozzles))
    width = 0.8 / max(len(groups), 1)

    for ax, metric, ylabel in [
        (axes[0], "uts_mpa", "UTS [MPa]"),
        (axes[1], "modulus_mpa", "Modulus (fit) [MPa]"),
    ]:
        agg = (
            summary.groupby(["nozzle_mm", "group"])[metric]
            .agg(["mean", "std"])
            .reset_index()
        )
        for i, grp in enumerate(groups):
            g = agg[agg["group"] == grp].set_index("nozzle_mm").reindex(nozzles)
            mode = grp.split("  [")[0]
            ax.bar(
                x + (i - (len(groups) - 1) / 2) * width,
                g["mean"],
                width,
                yerr=g["std"].fillna(0),
                label=grp,
                color=MODE_COLOR.get(mode, "0.4"),
                # Hatch separates campaigns sharing a print-mode colour.
                hatch="" if grp.endswith(f"[{sorted(summary['dataset'].unique())[0]}]") else "//",
                edgecolor="white",
                alpha=0.9,
                capsize=3,
            )
        ax.set_xticks(x)
        ax.set_xticklabels([f"{v:.1f}mm" for v in nozzles])
        ax.set_xlabel("Nozzle size")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3, axis="y")
        ax.legend(fontsize=7)

    fig.suptitle("Summary: UTS and Modulus by Nozzle Size, Print Mode and Campaign")
    fig.tight_layout()
    fig.savefig(outdir / "summary_uts_modulus.png", dpi=CFG.dpi)
    if CFG.show_plots:
        plt.show()
    plt.close(fig)


def plot_mode_comparison(summary: pd.DataFrame, outdir: Path):
    """Print-mode effect on UTS, per nozzle size, one panel per campaign."""
    datasets = sorted(summary["dataset"].unique())
    fig, axes = plt.subplots(1, len(datasets), figsize=(6.5 * len(datasets), 5), squeeze=False)

    for ax, dataset in zip(axes.flat, datasets):
        sub = summary[summary["dataset"] == dataset]
        for mode in sorted(sub["print_mode"].unique()):
            g = (
                sub[sub["print_mode"] == mode]
                .groupby("nozzle_mm")["uts_mpa"]
                .agg(["mean", "std"])
                .reset_index()
            )
            ax.errorbar(
                g["nozzle_mm"],
                g["mean"],
                yerr=g["std"].fillna(0),
                marker="o",
                capsize=4,
                lw=1.8,
                color=MODE_COLOR.get(mode, "0.4"),
                label=mode,
            )
        ax.set_title(f"campaign: {dataset}")
        ax.set_xlabel("Nozzle size [mm]")
        ax.set_ylabel("UTS [MPa]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)

    fig.suptitle("Print-mode effect on UTS across nozzle sizes")
    fig.tight_layout()
    fig.savefig(outdir / "print_mode_comparison.png", dpi=CFG.dpi)
    if CFG.show_plots:
        plt.show()
    plt.close(fig)


# ============================================================================
# Main
# ============================================================================
def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--data-root",
        action="append",
        default=None,
        help="Campaign folder to scan. Repeat for multiple; overrides CONFIG.data_roots.",
    )
    parser.add_argument("--output-dir", default=None, help="Override CONFIG.output_dir")
    args = parser.parse_args()

    if args.data_root:
        CFG.data_roots = tuple(Path(p) for p in args.data_root)
    if args.output_dir:
        CFG.output_dir = Path(args.output_dir)

    outdir = CFG.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    print("Reading campaigns:")
    for r in CFG.data_roots:
        print(f"  - {r}")
    print(
        f"Specimen area: {CFG.area_mm2:.1f} mm^2 "
        f"({CFG.specimen_width_mm} x {CFG.specimen_thickness_mm} mm)\n"
    )

    all_curves, summary, quality_notes = build_dataset(CFG.data_roots)

    summary.drop(columns=["file"]).to_csv(outdir / "run_summary.csv", index=False)
    all_curves.to_csv(outdir / "all_curves_tidy.csv", index=False)

    report_path = outdir / "data_quality_report.txt"
    with open(report_path, "w") as fh:
        fh.write(f"Tensile data quality report -- {len(summary)} runs processed\n")
        fh.write("=" * 78 + "\n\n")
        for note in quality_notes or ["No cleaning actions were needed.\n"]:
            fh.write(f"- {note}\n")

    flagged = summary[summary["flagged"]]
    print(f"Data quality report -> {report_path}")
    print(f"  {len(quality_notes)} note(s); {len(flagged)} run(s) flagged for review.")
    for _, r in flagged.iterrows():
        print(f"  ! FLAGGED: {r['dataset']} / {r['condition']} / {r['run']}")

    plot_per_nozzle(all_curves, summary, outdir)
    plot_grid(all_curves, summary, outdir)
    plot_summary_bars(summary, outdir)
    plot_mode_comparison(summary, outdir)

    print("\nPer-run summary:")
    cols = [
        "dataset", "condition", "run", "nozzle_mm", "print_mode",
        "uts_mpa", "strain_at_uts", "modulus_mpa", "final_recorded_strain", "flagged",
    ]
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(summary[cols].round(3).to_string(index=False))

    print("\nGroup means (UTS / modulus):")
    grp = (
        summary.groupby(["dataset", "nozzle_mm", "print_mode"])[["uts_mpa", "modulus_mpa"]]
        .agg(["mean", "std", "count"])
        .round(2)
    )
    with pd.option_context("display.width", 200):
        print(grp.to_string())

    print(f"\nAll outputs written to: {outdir.resolve()}")


if __name__ == "__main__":
    main()
