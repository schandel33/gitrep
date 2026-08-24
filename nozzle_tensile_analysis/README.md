# Nozzle Tensile Test Analysis

Stress-vs-strain analysis of tensile pulls across 4 nozzle sizes (0.6, 1.4,
2.2, 2.6 mm) x 2 print-gravity orientations (+1G upright / -1G inverted),
printed at 105 C with a 0.85 flow ratio.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run

```bash
python3 tensile_analysis.py \
  --data-root "/Users/siddharthchandel/Claude/Projects/Research/07_Experiments/Raw Data/20AUG2026"
```

`--data-root` and `--output-dir` are optional; both default to the values in
the `Config` dataclass at the top of `tensile_analysis.py`. `output_dir`
defaults to an `analysis_output/` folder next to the script itself
(anchored to the script's location, not whatever directory you happened to
launch it from) — so it always lands in the same place no matter how you
run it.

### Viewing plots: separate popup windows vs. an inline Plots panel

Plots always save as PNGs to the output directory regardless of how you run
the script. Whether they *also* pop up as windows or render inline in one
panel (like Spyder's Plots pane) depends on how you launch it:

- **Plain script execution** (VS Code's "Run Python File" button, or typing
  `python3 tensile_analysis.py` in a terminal) → each of the 6 figures pops
  up as its own separate OS window (via the Tk backend). This is the
  crash-safe path (see the code comment above the backend-selection block
  for why), but it's noisy if you're iterating a lot.
- **VS Code's Interactive Window** → get a Spyder-like experience instead:
  right-click anywhere in `tensile_analysis.py` and choose **"Run Current
  File in Interactive Window"** (or use the dropdown arrow next to the
  ▷ Run button at the top-right of the editor and pick it from there). This
  runs the whole script through a Jupyter kernel inside VS Code, and every
  figure the script produces collects in VS Code's **Plots** panel (the
  small icon in the Interactive Window's toolbar, or open via the Command
  Palette: "Jupyter: Focus on Plots View") — one place to scroll through
  every generated figure, exactly like Spyder's Plots pane. No code changes
  needed; the script already detects this (`get_ipython()`) and defers to
  VS Code's own backend instead of forcing separate windows.
  - First use may prompt you to install `ipykernel` in the active
    environment — accept that, it's what powers the Interactive Window.
  - Note: `--data-root`/`--output-dir` CLI flags aren't easily passed this
    way, so it runs with `Config`'s defaults. Edit those defaults directly
    in the script if you need a one-off different path.

Set `Config.show_plots = False` if you're batch-running many times and
don't want a window/inline plot appearing on every single run at all
(PNGs still get saved either way).

The script recurses through the data root looking for any CSV, so it doesn't
care how deeply nested the `Test Run .../DAQ...csv` folders are relative to
the `T###N##±1GFR###` condition folder.

## What it produces

- `run_summary.csv` — one row per test run: UTS, strain at UTS, fitted
  modulus, final recorded strain/stress, and data-cleaning stats.
- `all_curves_tidy.csv` — long-format stress/strain points for every run
  (decimated for plotting), tagged with nozzle size / gravity / run.
- `stress_strain_nozzle_<X>mm.png` — one figure per nozzle size, +1G vs -1G
  overlaid.
- `stress_strain_grid_all.png` — all four nozzle sizes in one grid.
- `summary_uts_modulus.png` — bar chart of UTS and modulus by condition.
- `data_quality_report.txt` — every automatic correction the script made
  (see below), so nothing is silently changed.

## The "local environment" (`Config` in `tensile_analysis.py`)

All tunable values used to parse and clean the data live in one
`@dataclass Config` block at the top of the script — nothing else in the
code hardcodes a path, threshold, or geometry number:

- `data_root`, `output_dir`
- `specimen_width_mm` / `specimen_thickness_mm` (custom dogbone: 20 mm x 4 mm
  gauge section, constant across all runs → area = 80 mm²)
- CSV header-detection marker and expected column names
- `crosshead_zero_threshold_mm` — defines "pre-motion" samples used to
  baseline-correct the load cell
- `max_extensometer_step_mm` — sensor-glitch trim threshold (see below)
- `modulus_strain_lo` / `modulus_strain_hi` — the engineering-strain window
  used to fit the elastic modulus by linear regression
- `outlier_strain_ratio` — how far a run's recorded strain range can exceed
  its replicates before it gets flagged (not dropped) for manual review

Edit these and rerun — no other code changes needed for a different
specimen geometry, sampling rig, or nomenclature.

## Data-quality issues found and how they're handled

1. **Load-cell zero offset.** 10 of 17 runs record small negative Load
   values for the first handful of samples (up to ~180, always <0.5% of the
   file) before the crosshead starts moving — load-cell noise around zero,
   not real load. The script computes the mean Load over all
   pre-motion samples (`Crosshead < crosshead_zero_threshold_mm`) and
   subtracts it from the whole run.
2. **NaN rows.** Two runs have a handful of NaN rows (extensometer sensor
   dropout): `Test Run 63` (272 rows, at the very end) and `Test Run 71` (1
   row, the very first sample). Dropped.
3. **Extensometer sensor glitch.** In `Test Run 63`, immediately before the
   sensor drops to NaN, it logs one physically-impossible instantaneous
   jump (Extensometer drops ~36 mm in a single 0.01 s sample — normal motion
   is ~0.01-0.05 mm/sample). The script detects any single-sample jump
   larger than `max_extensometer_step_mm` and trims that sample and
   everything after it.
4. **No true fracture signal.** None of the 17 runs show a classic brittle
   load collapse — every recording simply ends while load is still
   ~70-90% of its peak (the print material draws/yields rather than
   snapping within the recorded window). So "elongation at break" isn't
   meaningful here; the script reports **peak stress (UTS) and strain at
   peak**, plus the **final recorded strain/stress** labeled explicitly as
   "recorded", not "at break".
5. **`Test Run 73` (2.2 mm, +1G) is a genuine outlier in duration**, not a
   sensor artifact: it runs ~2x longer than its sibling replicate (80,750
   vs. 41,629 samples) and reaches ~2.5x the final strain, because the
   machine wasn't stopped where its replicate was. It's automatically
   flagged in `run_summary.csv` (`outlier_extension=True`) and drawn with
   reduced opacity in the per-nozzle plots — it's kept in the dataset (real
   data, not garbage) but should be sanity-checked manually before treating
   it as equivalent to its 2-3x shorter replicates.
6. **Gauge length is per-specimen, not a fixed nominal value.** The
   Extensometer channel reports an absolute position (~24.6-25.3 mm at the
   start of each run, not a round number), so strain is computed per run as
   `(Extensometer − L0) / L0`, with `L0` = that run's first valid
   Extensometer reading — not a single hardcoded gauge length.
7. Column headers in the raw CSVs have trailing spaces (`"Crosshead "`,
   etc.) and the files open with a BOM + an unterminated quoted preamble
   line — the parser strips both.
