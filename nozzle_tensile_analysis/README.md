# Nozzle Tensile Test Analysis

Stress-vs-strain analysis of tensile pulls across nozzle sizes and print
modes, combining multiple experiment campaigns printed at 105 C.

| Campaign | Nozzles | Print modes | Flow ratio | Runs |
|---|---|---|---|---|
| `Raw Data/20AUG2026` | 0.6, 1.4, 2.2, 2.6 mm | +1G, -1G | 0.85 | 17 |
| `08_Tensile_Tests` | 0.6, 1.4 mm | +1G, -1G, CLS | not encoded | 18 |

## Folder nomenclature

Two conventions are supported; the script detects which applies per folder.

**New** -- `T105N06-1GFR085`
- `T105` temperature 105 C; `N06` 0.6 mm nozzle (tenths implied)
- `-1G` inverted (`+1G` upright); `FR085` flow ratio 0.85

**Old** -- `T105N0_6_1G`, `T105N0_6(1G)`, `T105N0_6CLS`
- `N0_6` 0.6 mm nozzle (underscore is the decimal point)
- `_1G` -> +1G upright, `(1G)` -> -1G inverted, `CLS` -> rotating print mode
- No flow ratio encoded, so those runs carry `flow_ratio = NaN`

Note `(1G)` is tested before the bare `1G` when matching, since the inverted
marker contains the upright one as a substring.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run

```bash
python3 tensile_analysis.py \
  --data-root "/Users/siddharthchandel/Claude/Projects/Research/07_Experiments/Raw Data/20AUG2026" \
  --data-root "/Users/siddharthchandel/Claude/Projects/Research/08_Tensile_Tests"
```

`--data-root` is repeatable -- each campaign folder becomes a `dataset`
label (the folder's own name) so campaigns stay distinguishable in the
summary table and plots. With no arguments it uses `CONFIG.data_roots`,
which already lists both campaigns.

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

- `run_summary.csv` — one row per test run: campaign, condition, nozzle,
  print mode, UTS, strain at UTS, fitted modulus, final recorded
  strain/stress, cleaning stats, and QC flags.
- `all_curves_tidy.csv` — long-format stress/strain points for every run
  (decimated for plotting), tagged with campaign / nozzle / print mode / run.
- `stress_strain_nozzle_<X>mm.png` — one figure per nozzle size. Colour =
  print mode, line style = campaign, faded thin lines = QC-flagged runs.
- `stress_strain_grid_all.png` — every nozzle size in one grid.
- `summary_uts_modulus.png` — UTS and modulus bars by nozzle, grouped by
  print mode + campaign (hatching separates campaigns sharing a colour).
- `print_mode_comparison.png` — UTS vs. nozzle size per print mode, one
  panel per campaign; the clearest view of the print-mode effect.
- `data_quality_report.txt` — every automatic correction and every QC flag,
  so nothing is silently changed or dropped.

## The "local environment" (`Config` in `tensile_analysis.py`)

All tunable values used to parse and clean the data live in one
`@dataclass Config` block at the top of the script — nothing else in the
code hardcodes a path, threshold, or geometry number:

- `data_roots` (a tuple — one entry per campaign), `output_dir`
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
- `max_abs_load_offset_n` — baseline offset above which the load cell was so
  badly zeroed that every reading in the run is suspect
- `max_nan_fraction` — share of rows lost to sensor dropout before a run is
  flagged untrustworthy
- `exclude_runs` — run labels to drop entirely; empty by default so nothing
  disappears without an explicit, recorded decision
- `skip_path_substrings` — ignores the `__MACOSX/` resource-fork tree that
  macOS zip archives carry alongside real CSVs

Edit these and rerun — no other code changes needed for a different
specimen geometry, sampling rig, or nomenclature.

## Data-quality issues found and how they're handled

Across all **35 runs** (17 new + 18 old). Everything below is applied
automatically and logged to `data_quality_report.txt`.

1. **Load-cell zero offset.** Most runs in both campaigns record a nonzero
   Load for the first handful of samples before the crosshead starts moving
   — load-cell noise around zero, not real load. The script averages the
   pre-motion samples (`Crosshead < crosshead_zero_threshold_mm`) and
   subtracts that offset from the whole run. Typical corrections are a few
   newtons; anything beyond `max_abs_load_offset_n` is flagged instead of
   quietly trusted.
2. **NaN rows** from extensometer dropout are dropped: `Test Run 63` (272
   rows at the very end), `Test Run 71` (1 row, first sample), `Test Run 62`
   (1 row), and `Test Run 52` (18,015 rows — see below).
3. **Extensometer sensor glitch.** In `Test Run 63`, immediately before the
   sensor drops to NaN, it logs one physically-impossible instantaneous jump
   (~36 mm in a single 0.01 s sample; normal motion is ~0.01-0.05 mm/sample).
   Any single-sample jump beyond `max_extensometer_step_mm` triggers
   trimming of that sample and everything after it.
4. **No true fracture signal in any run, either campaign.** Every recording
   ends while load is still ~70-90% of its peak — the material draws/yields
   rather than snapping within the recorded window. "Elongation at break" is
   therefore not measurable here, so the script reports **UTS and strain at
   UTS**, plus **final recorded strain/stress** labelled explicitly as
   "recorded", never "at break".
5. **Gauge length is per-specimen, not nominal.** Extensometer reports an
   absolute position starting anywhere from ~23.0-25.9 mm across the two
   campaigns, so strain is `(Extensometer − L0) / L0` with `L0` = that run's
   own first valid reading.
6. **Parser tolerances**: headers carry trailing spaces (`"Crosshead "`),
   files open with a BOM plus an unterminated quoted preamble line, and
   macOS zips add a `__MACOSX/` tree of fake CSVs. All handled.

### Runs flagged for manual review

Flagged runs stay in the dataset and plots (drawn faded/thin) — they are
real data, not garbage — but should not be treated as equivalent to their
replicates without a look:

- **`Test Run 52`** (old campaign, `T105N0_6CLS`) — **the worst run in
  either campaign, and the one to check first.** Three independent problems:
  its load cell was zeroed at **−306 N** (so every load reading in the run
  is suspect, not merely the first samples); **12.8%** of its rows are lost
  to extensometer dropout; and it ran to **3.13 strain** vs ~1.5 for its
  replicates, with the crosshead reaching 117 mm against ~30-40 mm typical.
  Consider excluding it via `Config.exclude_runs` once you have decided.
- **`Test Run 73`** (new campaign, `T105N22+1GFR085`) — ran ~2× longer than
  its replicate (80,750 vs 41,629 samples) to ~2.2× the final strain,
  because the machine wasn't stopped where its replicate was. The curve
  itself looks healthy.

## Caveats worth knowing before publishing

- **Specimen geometry is assumed identical across campaigns** (20 mm × 4 mm
  → 80 mm²). If the older experiments used different dimensions, that whole
  campaign's stress axis is rescaled by a constant and the cross-campaign
  comparison is invalid. Split `Config` into a per-dataset geometry lookup
  if so.
- **Flow ratio is not encoded in the old campaign's folder names.** The new
  campaign is all FR 0.85. If the old runs used a different flow ratio, then
  campaign and flow ratio are confounded — a cross-campaign difference
  cannot be attributed to print mode alone.
- **The CLS rotation speed is not recorded anywhere in the data.** All CLS
  runs are pooled under one label regardless of RPM. If more than one RPM
  was used, they need separating before the CLS numbers mean anything.
- **Replicate counts are uneven** (n=3 for the old campaign, n=2-3 for the
  new). Small-n additive-manufacturing studies draw reviewer fire for
  p-value-only reporting — report effect sizes alongside any significance
  test.
