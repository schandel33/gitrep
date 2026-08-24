#!/usr/bin/env python3
"""
Tensile test analysis: stress vs. strain curves across 3D-printed nozzle
sizes and print-gravity orientations.

Folder / filename nomenclature produced by the MTS DAQ export, e.g.:

    T105N06-1GFR085/
        Test Run 74 8-21-2026 4 21 37 PM/
            DAQ- Crosshead, ... - (Timed).csv

    T105  -> print temperature, 105 C
    N06   -> 0.6 mm nozzle
    -1G   -> inverted printer orientation while printing (+1G = upright)
    FR085 -> flow ratio 0.85

Each DAQ csv has 4 metadata lines, a blank line, a header row, and a units
row before the numeric data starts:

    "Crosshead ","Load ","Time ","Extensometer "
    "mm","N","sec","mm"

Run this directly:

    python3 tensile_analysis.py
    python3 tensile_analysis.py --data-root "/path/to/Raw Data/20AUG2026"

All the "important variables" -- paths, specimen geometry, cleaning
thresholds, fit windows -- live in one place: the CONFIG object below.
Nothing else in the script hardcodes a magic number; edit CONFIG and rerun.
"""

from __future__ import annotations

import argparse
import dataclasses
import glob
import io
import re
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def _running_inside_ipython_kernel() -> bool:
    """True in Spyder/Jupyter (they pre-configure their own matplotlib
    backend via %matplotlib magic before this script's code runs), False for
    a plain `python tensile_analysis.py` invocation (a terminal, VS Code's
    integrated terminal, a cloud/headless session)."""
    try:
        return get_ipython() is not None  # noqa: F821 -- injected by IPython at runtime
    except NameError:
        return False


if not _running_inside_ipython_kernel():
    # Plain script execution: don't let matplotlib auto-select its native
    # backend. On very new Python releases (e.g. 3.14), matplotlib's native
    # macOS backend (a compiled C extension) can crash outright --
    # `SystemError: NULL object passed to Py_BuildValue` followed by a
    # segfault -- before the backend's C code has caught up with the new
    # CPython C-API. TkAgg doesn't touch Apple's Cocoa APIs at all and has
    # been stable across Python versions for decades, so try it first; fall
    # back to the non-interactive, save-only Agg backend if Tk isn't present
    # at all (headless CI/cloud sessions, minimal Python builds).
    try:
        matplotlib.use("TkAgg")
    except Exception:
        matplotlib.use("Agg")
# else: running inside Spyder's/Jupyter's own kernel -- leave the backend it
# already configured alone (its Graphics preference / %matplotlib magic).


# ============================================================================
# CONFIG -- the "local environment" for this analysis.
# Every knob that affects parsing, cleaning, or curve-fitting lives here.
# ============================================================================
@dataclasses.dataclass
class Config:
    # --- where the raw data lives / where results get written -------------
    data_root: Path = Path(
        "/Users/siddharthchandel/Claude/Projects/Research/07_Experiments/"
        "Raw Data/20AUG2026"
    )
    # Anchored to this script's own location (not the current working
    # directory) so it always lands in the same place regardless of how or
    # from where you run the script -- VS Code's "Run Python File" button,
    # its integrated terminal, and a Spyder console can all start with a
    # different cwd, which otherwise silently creates a different
    # analysis_output/ folder each time.
    output_dir: Path = Path(__file__).resolve().parent / "analysis_output"

    # --- specimen geometry ---------------------------------------------------
    # Custom dogbone, constant across every nozzle size / gravity condition.
    specimen_width_mm: float = 20.0
    specimen_thickness_mm: float = 4.0

    @property
    def area_mm2(self) -> float:
        """Cross-sectional area of the gauge section (mm^2). Load[N]/Area = Stress[MPa]."""
        return self.specimen_width_mm * self.specimen_thickness_mm

    # --- CSV parsing -----------------------------------------------------
    # The numeric header row is the first line that starts with this marker
    # (tolerant to the stray BOM / unterminated-quote line MTS puts at the top).
    header_marker: str = '"Crosshead'
    expected_columns: tuple = ("Crosshead", "Load", "Time", "Extensometer")

    # --- data cleaning -----------------------------------------------------
    # Pre-pull load-cell noise: several runs record small negative Load
    # values for a handful of samples before the crosshead actually starts
    # moving. We zero this out by subtracting the mean Load measured while
    # Crosshead is still below this threshold (mm).
    crosshead_zero_threshold_mm: float = 0.01
    min_baseline_samples: int = 5  # need at least this many pre-motion samples to baseline-correct

    # --- strain reference ---------------------------------------------------
    # Gauge length L0 is *not* a fixed nominal value in this rig -- the
    # Extensometer column reports an absolute position that starts a little
    # differently for every specimen (~24.6-25.3 mm here). L0 is taken as
    # each run's first valid Extensometer reading.
    strain_from: str = "Extensometer"  # "Extensometer" (recommended) or "Crosshead"

    # --- elastic-modulus fit window (engineering strain) --------------------
    modulus_strain_lo: float = 0.0005
    modulus_strain_hi: float = 0.01
    modulus_min_points: int = 20

    # --- sensor-glitch trimming ------------------------------------------
    # Extensometer occasionally faults right before it drops out to NaN,
    # logging one physically-impossible instantaneous jump (tens of mm in a
    # single 0.01s sample -- normal motion is ~0.01-0.05 mm/sample). Any run
    # with a step larger than this (mm) has itself and everything after it
    # trimmed off as garbage.
    max_extensometer_step_mm: float = 1.0

    # --- outlier flagging ----------------------------------------------------
    # A run is flagged (not dropped) if its total recorded strain range is
    # this many times its condition group's median -- catches runs like
    # "Test Run 73" that were left running far longer than their replicates.
    outlier_strain_ratio: float = 1.75

    # --- filename / folder parsing -------------------------------------------
    condition_pattern: str = (
        r"T(?P<temp_c>\d+)N(?P<nozzle_mm_x10>\d+)(?P<gravity>[+-])1GFR(?P<flow_ratio_x100>\d+)"
    )
    csv_glob: str = "**/*.csv"

    # --- plotting -------------------------------------------------------------
    dpi: int = 150
    max_curve_points: int = 4000  # decimate long curves for lighter plot files
    show_plots: bool = True  # call plt.show() (pops up / renders inline in Spyder, VS Code, Jupyter, etc.)


CFG = Config()


# ============================================================================
# Parsing
# ============================================================================
def parse_condition(path: Path) -> dict:
    """Pull temp/nozzle/gravity/flow-ratio out of whichever path component matches."""
    for part in path.parts:
        m = re.search(CFG.condition_pattern, part)
        if m:
            g = m.groupdict()
            return {
                "condition": m.group(0),
                "temp_c": int(g["temp_c"]),
                "nozzle_mm": int(g["nozzle_mm_x10"]) / 10.0,
                "gravity": "+1G (upright)" if g["gravity"] == "+" else "-1G (inverted)",
                "flow_ratio": int(g["flow_ratio_x100"]) / 100.0,
            }
    raise ValueError(f"Could not parse condition nomenclature from path: {path}")


def parse_run_label(path: Path) -> str:
    for part in path.parts:
        if part.lower().startswith("test run"):
            # "Test Run 65 8-21-2026 4 26 45 PM" -> "Test Run 65"
            m = re.match(r"(Test Run \d+)", part)
            return m.group(1) if m else part
    return path.parent.name


def load_daq_csv(path: Path) -> pd.DataFrame:
    """Tolerant parser for the MTS DAQ export: skips the metadata preamble,
    strips whitespace from headers, and coerces everything to numeric."""
    with open(path, encoding="utf-8-sig", errors="replace") as fh:
        lines = fh.readlines()

    header_idx = next(
        (i for i, l in enumerate(lines) if l.strip().startswith(CFG.header_marker)), None
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

    # 1) Baseline-correct Load using samples before the crosshead actually
    #    starts moving (removes the small negative pre-pull noise).
    pre_motion = df_raw["Crosshead"] < CFG.crosshead_zero_threshold_mm
    n_baseline = int(pre_motion.sum())
    load_offset = float(df_raw.loc[pre_motion, "Load"].mean()) if n_baseline >= CFG.min_baseline_samples else 0.0

    df = df_raw.copy()
    df["Load"] = df["Load"] - load_offset

    # 2) Drop rows with any NaN (edge artifacts: extensometer dropout at the
    #    very start/end of a handful of runs).
    n_before = len(df)
    df = df.dropna().reset_index(drop=True)
    n_nan_dropped = n_before - len(df)

    # 3) Trim trailing sensor-glitch garbage: a physically-impossible
    #    instantaneous jump in Extensometer (see CONFIG.max_extensometer_step_mm)
    #    means that sample and everything after it is not trustworthy.
    n_before_glitch = len(df)
    ext_step = df["Extensometer"].diff().abs()
    glitch_idx = ext_step[ext_step > CFG.max_extensometer_step_mm].index
    if len(glitch_idx):
        df = df.iloc[: glitch_idx[0]].reset_index(drop=True)
    n_glitch_trimmed = n_before_glitch - len(df)

    # 4) Gauge length reference for strain.
    ref_col = CFG.strain_from
    gauge_length = float(df[ref_col].iloc[0])

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
# Per-run summary metrics
# ============================================================================
def modulus_fit(df: pd.DataFrame) -> float | None:
    lo, hi = CFG.modulus_strain_lo, CFG.modulus_strain_hi
    seg = df[(df["Strain"] >= lo) & (df["Strain"] <= hi)]
    if len(seg) < CFG.modulus_min_points:
        return None
    slope, _intercept = np.polyfit(seg["Strain"], seg["Stress_MPa"], 1)
    return float(slope)


def summarize_run(df: pd.DataFrame, clean: CleanResult) -> dict:
    peak_i = int(df["Stress_MPa"].values.argmax())
    peak_stress = float(df["Stress_MPa"].iloc[peak_i])
    strain_at_peak = float(df["Strain"].iloc[peak_i])
    final_strain = float(df["Strain"].iloc[-1])
    final_stress = float(df["Stress_MPa"].iloc[-1])
    return {
        "n_points": len(df),
        "n_raw_points": clean.n_raw,
        "n_nan_dropped": clean.n_nan_dropped,
        "n_glitch_trimmed": clean.n_glitch_trimmed,
        "load_offset_n": clean.load_offset_n,
        "gauge_length_mm": clean.gauge_length_mm,
        "uts_mpa": peak_stress,
        "strain_at_uts": strain_at_peak,
        "modulus_mpa": modulus_fit(df),
        "final_recorded_strain": final_strain,
        "final_recorded_stress_mpa": final_stress,
        "final_stress_frac_of_peak": final_stress / peak_stress if peak_stress else np.nan,
    }


# ============================================================================
# Driver: discover, load, clean, summarize every run
# ============================================================================
def discover_runs(data_root: Path) -> list[Path]:
    files = sorted(Path(p) for p in glob.glob(str(data_root / CFG.csv_glob), recursive=True))
    if not files:
        raise FileNotFoundError(
            f"No CSV files found under {data_root}. Check --data-root / CONFIG.data_root."
        )
    return files


def build_dataset(data_root: Path):
    curves = []       # list of tidy per-point DataFrames (decimated for plotting)
    summary_rows = []
    quality_notes = []

    for f in discover_runs(data_root):
        cond = parse_condition(f)
        run_label = parse_run_label(f)

        raw = load_daq_csv(f)
        clean = clean_run(raw)
        df = compute_stress_strain(clean)
        summ = summarize_run(df, clean)

        row = {**cond, "run": run_label, "file": str(f), **summ}
        summary_rows.append(row)

        if clean.n_baseline_samples >= CFG.min_baseline_samples and abs(clean.load_offset_n) > 0.5:
            quality_notes.append(
                f"[{cond['condition']} / {run_label}] Load zero-offset corrected by "
                f"{clean.load_offset_n:+.2f} N (from {clean.n_baseline_samples} pre-motion samples)."
            )
        if clean.n_nan_dropped:
            quality_notes.append(
                f"[{cond['condition']} / {run_label}] Dropped {clean.n_nan_dropped} row(s) with "
                f"NaN values (sensor dropout)."
            )
        if clean.n_glitch_trimmed:
            quality_notes.append(
                f"[{cond['condition']} / {run_label}] Trimmed {clean.n_glitch_trimmed} trailing row(s) "
                f"after a >{CFG.max_extensometer_step_mm}mm single-sample Extensometer jump "
                f"(sensor glitch, physically impossible motion)."
            )

        curve = df[["Strain", "Stress_MPa"]].copy()
        if len(curve) > CFG.max_curve_points:
            curve = curve.iloc[:: max(1, len(curve) // CFG.max_curve_points)]
        curve["condition"] = cond["condition"]
        curve["nozzle_mm"] = cond["nozzle_mm"]
        curve["gravity"] = cond["gravity"]
        curve["run"] = run_label
        curves.append(curve)

    summary = pd.DataFrame(summary_rows)

    # Flag runs whose recorded strain range is way beyond its *other*
    # replicates in the same condition (leave-one-out, so a single outlier
    # replicate can't drag its own reference up -- group medians/means with
    # only 2-3 replicates per condition are too easily pulled by the outlier
    # itself). Groups with only one run can't be checked this way.
    summary["strain_range"] = summary["final_recorded_strain"]
    summary["outlier_extension"] = False
    summary["outlier_reference"] = np.nan
    for cond, idx in summary.groupby("condition").groups.items():
        idx = list(idx)
        if len(idx) < 2:
            continue
        for i in idx:
            others = [j for j in idx if j != i]
            ref = summary.loc[others, "strain_range"].max()
            summary.loc[i, "outlier_reference"] = ref
            if ref > 0 and summary.loc[i, "strain_range"] > CFG.outlier_strain_ratio * ref:
                summary.loc[i, "outlier_extension"] = True

    for _, r in summary[summary["outlier_extension"]].iterrows():
        quality_notes.append(
            f"[{r['condition']} / {r['run']}] FLAGGED: recorded strain range "
            f"({r['strain_range']:.2f}) is >{CFG.outlier_strain_ratio}x its other replicate(s) "
            f"({r['outlier_reference']:.2f}) -- likely stopped much later than its replicates. "
            f"Included in plots/summary but review manually."
        )

    all_curves = pd.concat(curves, ignore_index=True)
    return all_curves, summary, quality_notes


# ============================================================================
# Plotting
# ============================================================================
GRAVITY_STYLE = {
    "+1G (upright)": dict(color="#2b6cb0", ls="-"),
    "-1G (inverted)": dict(color="#c0392b", ls="--"),
}


def plot_per_nozzle(all_curves: pd.DataFrame, summary: pd.DataFrame, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    nozzles = sorted(all_curves["nozzle_mm"].unique())

    for nz in nozzles:
        fig, ax = plt.subplots(figsize=(7, 5.5))
        sub = all_curves[all_curves["nozzle_mm"] == nz]
        for (gravity, run), g in sub.groupby(["gravity", "run"]):
            style = GRAVITY_STYLE.get(gravity, {})
            is_outlier = summary.loc[
                (summary["run"] == run) & (summary["nozzle_mm"] == nz), "outlier_extension"
            ]
            alpha = 0.5 if (len(is_outlier) and bool(is_outlier.iloc[0])) else 0.9
            lw = 1.0 if (len(is_outlier) and bool(is_outlier.iloc[0])) else 1.6
            label = f"{gravity} - {run}" + (" (flagged: long run)" if alpha < 0.9 else "")
            ax.plot(g["Strain"], g["Stress_MPa"], label=label, alpha=alpha, lw=lw, **style)

        ax.set_xlabel("Engineering strain [mm/mm]")
        ax.set_ylabel("Engineering stress [MPa]")
        ax.set_title(f"Stress vs. Strain -- Nozzle {nz:.1f} mm (T105, FR0.85)")
        ax.legend(fontsize=8, loc="best")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fname = outdir / f"stress_strain_nozzle_{nz:.1f}mm.png"
        fig.savefig(fname, dpi=CFG.dpi)
        if CFG.show_plots:
            plt.show()
        plt.close(fig)


def plot_grid(all_curves: pd.DataFrame, summary: pd.DataFrame, outdir: Path):
    nozzles = sorted(all_curves["nozzle_mm"].unique())
    n = len(nozzles)
    ncols = 2
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.5 * ncols, 4.5 * nrows), squeeze=False)

    for ax, nz in zip(axes.flat, nozzles):
        sub = all_curves[all_curves["nozzle_mm"] == nz]
        for (gravity, run), g in sub.groupby(["gravity", "run"]):
            style = GRAVITY_STYLE.get(gravity, {})
            ax.plot(g["Strain"], g["Stress_MPa"], alpha=0.8, lw=1.3, **style)
        ax.set_title(f"Nozzle {nz:.1f} mm")
        ax.set_xlabel("Strain [mm/mm]")
        ax.set_ylabel("Stress [MPa]")
        ax.grid(alpha=0.3)

    for ax in axes.flat[n:]:
        ax.axis("off")

    handles = [plt.Line2D([], [], label=k, **v) for k, v in GRAVITY_STYLE.items()]
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.suptitle("Stress-Strain by Nozzle Size and Print-Gravity Orientation (T105, FR0.85)", y=0.99)
    fig.legend(handles=handles, loc="upper center", ncol=2, fontsize=10, bbox_to_anchor=(0.5, 0.955))
    fig.savefig(outdir / "stress_strain_grid_all.png", dpi=CFG.dpi, bbox_inches="tight")
    if CFG.show_plots:
        plt.show()
    plt.close(fig)


def plot_summary_bars(summary: pd.DataFrame, outdir: Path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, metric, ylabel in [
        (axes[0], "uts_mpa", "UTS [MPa]"),
        (axes[1], "modulus_mpa", "Modulus (fit) [MPa]"),
    ]:
        agg = summary.groupby(["nozzle_mm", "gravity"])[metric].agg(["mean", "std", "count"]).reset_index()
        nozzles = sorted(agg["nozzle_mm"].unique())
        width = 0.35
        x = np.arange(len(nozzles))
        for i, gravity in enumerate(sorted(agg["gravity"].unique())):
            g = agg[agg["gravity"] == gravity].set_index("nozzle_mm").reindex(nozzles)
            style = GRAVITY_STYLE.get(gravity, {})
            ax.bar(
                x + (i - 0.5) * width,
                g["mean"],
                width,
                yerr=g["std"].fillna(0),
                label=gravity,
                color=style.get("color"),
                alpha=0.85,
                capsize=4,
            )
        ax.set_xticks(x)
        ax.set_xticklabels([f"{n:.1f}mm" for n in nozzles])
        ax.set_xlabel("Nozzle size")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3, axis="y")
        ax.legend(fontsize=8)

    fig.suptitle("Summary: UTS and Modulus by Nozzle Size / Gravity Orientation")
    fig.tight_layout()
    fig.savefig(outdir / "summary_uts_modulus.png", dpi=CFG.dpi)
    if CFG.show_plots:
        plt.show()
    plt.close(fig)


# ============================================================================
# Main
# ============================================================================
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-root", type=str, default=None, help="Override CONFIG.data_root")
    parser.add_argument("--output-dir", type=str, default=None, help="Override CONFIG.output_dir")
    args = parser.parse_args()

    if args.data_root:
        CFG.data_root = Path(args.data_root)
    if args.output_dir:
        CFG.output_dir = Path(args.output_dir)

    data_root = CFG.data_root
    outdir = CFG.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading data from: {data_root}")
    print(f"Specimen area: {CFG.area_mm2:.1f} mm^2 ({CFG.specimen_width_mm} x {CFG.specimen_thickness_mm} mm)")

    all_curves, summary, quality_notes = build_dataset(data_root)

    # --- save tidy outputs ---
    summary_path = outdir / "run_summary.csv"
    curves_path = outdir / "all_curves_tidy.csv"
    summary.drop(columns=["file"]).to_csv(summary_path, index=False)
    all_curves.to_csv(curves_path, index=False)

    # --- data quality report ---
    report_path = outdir / "data_quality_report.txt"
    with open(report_path, "w") as fh:
        fh.write(f"Tensile data quality report -- {len(summary)} runs processed\n")
        fh.write("=" * 70 + "\n\n")
        if quality_notes:
            for note in quality_notes:
                fh.write(f"- {note}\n")
        else:
            fh.write("No cleaning actions were needed.\n")
    print(f"\nData quality report -> {report_path}")
    for note in quality_notes:
        print(f"  - {note}")

    # --- plots ---
    plot_per_nozzle(all_curves, summary, outdir)
    plot_grid(all_curves, summary, outdir)
    plot_summary_bars(summary, outdir)

    # --- console summary table ---
    print("\nPer-run summary:")
    cols = ["condition", "run", "uts_mpa", "strain_at_uts", "modulus_mpa",
            "final_recorded_strain", "outlier_extension"]
    with pd.option_context("display.width", 140, "display.max_columns", 20):
        print(summary[cols].round(3).to_string(index=False))

    print(f"\nAll outputs written to: {outdir.resolve()}")


if __name__ == "__main__":
    main()
