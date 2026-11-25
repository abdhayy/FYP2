# Change log and comparison: vibration replay and bearing fault tool

One living file. Every time the code changes, add a new entry at the top of "Change log" and a row in "Version history". Keep old entries.

Project: Fault Feature Extraction in Vibration Signal of Rotational Machinery Using Entropy Measures (Ammar Na'im bin Borhan, UiTM, supervisor Dr. Abdul Malek Abdul Wahab).
Hardware: Teensy 4.1 + ADXL326 (A0) + MATLAB for acquisition. This tool replays the saved Excel recordings (Fs = 2000 Hz, 60 groups x 5 s).

---

## Version history

| File | What it adds | Date |
|---|---|---|
| `plot_dataset_time_fft.py` | Replay of the Excel dataset: time domain + FFT with Mode 1/2/3 labels | not recorded |
| `plot_dataset_time_fft_entropy.py` | + second window every 10 s: VMD and 7 entropies per mode, multiscale SampEn, optional Excel log | not recorded |
| `..._15oct26*.py` | + `--bearing-sim`: simulated bearing signal with BPFO, BPFI, BSF, FTF | 15 Oct 2026 |
| `plot_dataset_time_fft_entropy_bearing_30oct26.py` | + `--bearing`: replay of YOUR dataset with the four fault components, their sum, envelope FFT, condition banner | 30 Oct 2026 |
| `plot_dataset_time_fft_entropy_bearing_10nov26.py` | + fault-evidence graphics and VMD of the combined signal (third window) | 10 Nov 2026 |
| `plot_dataset_time_fft_entropy_bearing_10nov26.py` (comments revision) | Same program, technique guide and comments added. Program logic identical (checked by comparing syntax trees) | 15 Nov 2026 |
| `plot_dataset_time_fft_entropy_bearing_23nov26.py` | Three-window redesign of `--bearing`: window 2 = new entropy window with healthy-range gauges and entropy verdict; window 3 = VMD without the per-mode entropy table | 23 Nov 2026 |
| `export_dashboard_data.py` + `dashboard-app/` | React web dashboard of the three windows; results exported to JSON by importing the newest analysis script. Analysis code unchanged | 23 Nov 2026 |
| `mafaulda_analysis.py` + MaFaulDa tab | Same method run on the public MaFaulDa dataset (cage, outer race, ball faults; healthy range from speed-matched normal files; settings tuned on one half, tested on the other) | 25 Nov 2026 |

---

## Change log

### 25 Nov 2026: MaFaulDa dataset (`mafaulda_analysis.py`, dashboard tab "MaFaulDa dataset")

**Why:** test the method on independent data with more fault types (the own rig has only a cage fault), as listed under "Known limits".

**What the script does:** every MaFaulDa CSV (50 kHz, 5 s) -> shaft speed from the tachometer -> fault lines from the MaFaulDa bearing factors (FTF 0.375, BSF 1.871, BPFO 2.998, BPFI 5.002 x fr) -> evidence in dB (`diagnose`) for 6 accelerometers x 5 bands, and the 7 entropies (`entropy_of_window`, accelerometer decimated to 2000 Hz, entropies at 500 Hz) -> z-scores against the 10 normal files with the closest shaft speed -> best sensor/band chosen on the "tune" half -> detection reported on the "test" half.

**First look (4 files near 24 Hz, one normal file, NOT a result):** the 600-900 Hz band of the own rig does not transfer: the fault-minus-normal margin of the true fault line was only 1 to 8 dB depending on channel and band (best: outer race +8.1 dB, column 6, 3-6 kHz). Column 6 shows a +24.5 dB FTF line in the NORMAL file too, so a single fault file can look faulty without a reference. Dispersion entropy changes with every fault, but the direction depends on the channel (down and up).

**Status:** script and dashboard tab tested with these 4 files (pipeline only). The real result needs the full run on the laptop (all normal, underhang and overhang files).

**Known limits:** column order (tachometer, underhang axial/radial/tangential, overhang axial/radial/tangential, microphone) still to be confirmed against the dataset description; the healthy range needs enough normal files near each speed (49 normal files over 12-61 Hz).

---

### 23 Nov 2026: React dashboard (`dashboard-app/`) + exporter (`export_dashboard_data.py`)

**Why:** the user asked for a dashboard built as a React app that shows the output of the analysis code, and that keeps working when the Python code changes later.

**Comparison with the matplotlib windows**

| Item | matplotlib (23nov26 `--bearing`) | React dashboard |
|---|---|---|
| Where it runs | Desktop Python, three windows | Web browser, one page with tabs 1 / 2 / 3 + Overview |
| Analysis | Live, in background threads | Same functions, run beforehand by `export_dashboard_data.py`, saved as JSON |
| Window 1 | Rolling replay, updated every 0.5 s | One snapshot per 10 s window (at the end of the window) |
| Window 2 trend | Windows played so far | Same (windows up to the selected one) |
| Recordings | One per run | All four in one page, switchable, plus an overview table |
| Playback | SPACE pause, R restart | Play / pause, restart, slider, x1 to x20, keys SPACE / arrows / R |
| Analysis code | | Unchanged. The exporter imports the newest `plot_dataset_time_fft_entropy_bearing_*.py` |

**Settings used for the export:** `--balls 8`, reference = `Vibration_Dataset_Healthy_2000Hz.xlsx` (fault threshold and entropy range calibrated on it, threshold 7.44 dB), VMD combined, K = 18, 10 s windows.

**Test results (export of all four recordings, 30 windows each)**

| Recording | Spectrum verdict FAULTY | Entropy verdict FAULTY |
|---|---|---|
| Healthy | 0 / 30 | 0 / 30 |
| Healthy2 | 0 / 30 | 0 / 30 |
| Faulty1 | 30 / 30 (FTF) | 30 / 30 |
| Faulty2 | 30 / 30 (FTF) | 29 / 30 |

Window 30 of Faulty1 and Healthy matches the screenshots of the matplotlib windows (FTF 11.5 dB, BPFI 10.7 dB; Dispersion z -9.2, Fuzzy z -4.3; VMD r 0.999 / 0.999 / 0.995 / 0.999).

**Known limits:** the Healthy file is both the reference and one of the displayed recordings, so its own result is not an independent test. The exported data (about 38 MB of JSON) must be regenerated after every change to the analysis code.

---

### 23 Nov 2026: three-window redesign (`--bearing`): window 2 entropy, window 3 VMD

**Why:** the user did not understand the 7 x 18 entropy heat-map and asked for a window that shows which entropy separates Healthy from Faulty, with window 1 (live monitor) and VMD kept as they were.

**Comparison with the 10nov26 file (in `--bearing` mode)**

| Item | 10nov26 | 23nov26 |
|---|---|---|
| Window 1 | Live monitor, evidence graphics | Unchanged |
| Window 2 | Old entropy window: VMD K = 7 modes x 7 entropies heat-map, multiscale entropy | New: one gauge per entropy against the HEALTHY range (mean +/- 3 std), z-score, bar chart, window-by-window trend, multiscale curve vs healthy band, ENTROPY VERDICT banner |
| Window 3 | VMD of the combined signal + 7 entropies per mode (heat-map) + feature table | VMD of the combined signal + table of the modes + fault-feature table (evidence, RMS, r, share). Per-mode entropy removed (a mode is nearly a sine, its entropy mostly reflects its frequency) |
| Entropy computed on | VMD modes and raw signal | The whole 10 s window at 500 Hz (decimated with an anti-alias filter) |
| Entropy verdict | None | FAULTY when Dispersion AND Fuzzy are both more than 3 std below the Healthy mean |
| Healthy range | None | From `--reference` (every 10 s window of the Healthy file); built-in values if `--reference` is not given |
| `--log` | Old window-2 table | Window 2: 7 entropies, their z-scores, verdict, MSE scales |
| `--bearing-log` | Window 3 incl. all entropies | Window 3 without the entropy columns |
| Background jobs | Entropy job (old) + VMD job | Entropy job (new) + VMD job, held as one `_JobPair`; the VMD job no longer computes 21 x 7 entropies (about 5 s faster) |
| Untouched | | Plain mode (no `--bearing`) and `--bearing-sim` behave as before (smoke-tested) |

**New:** `ENT_*` constants, `entropy_of_window`, `calibrate_entropy_reference`, `entropy_verdict_from_z`, `analyze_entropy_window`, `BearingVMDSimulator.show_entropy` / `_draw_entropy_figure` / new `_draw_combined_figure`.
**Changed:** `analyze_combined_vmd` and `combined_result_to_row` (no entropies), `BearingVMDSimulator` (`__init__`, jobs, drawing), `make_bearing_dataset_simulator` (entropy calibration), help texts, technique guide section 8.

**Test results (headless run of the real `main()` on the Excel files; cut-offs from the Healthy file)**

| Recording | Windows analysed | Entropy verdict |
|---|---|---|
| Faulty1 (with `--reference`) | 16 | 16 FAULTY (Dispersion z -7 to -14, Fuzzy z -3.5 to -6.3) |
| Healthy2 (built-in cut-offs) | 14 | 14 HEALTHY |
Earlier fixed test on all windows (30 per recording): Dispersion 60/60 faulty and 0/30 healthy flagged, Fuzzy 59/60 and 0/30; Sample 30/60, Approximate 23/60, Spectral 15/60 (11/30 healthy), Shannon 8/60 (4/30 healthy), Permutation 0/60. These counts are a fixed note in the code (`ENT_TEST`), not recomputed.

**Known limits:** healthy range from ONE recording; one fault type (cage); entropy cut-offs and the counts above come from the same four recordings, so this is not yet an independent test; entropy cannot name the fault type.

---

### 10 Nov 2026: fault evidence + VMD of the combined signal

**Why:** the user asked (1) how can we see in the graph that it is an FTF fault and what defines the fault, (2) to take the time-domain combined signal (all four faults) and separate it again with VMD, in line with the FYP title.

**Comparison with the 30oct26 file**

| Item | 30oct26 | 10nov26 |
|---|---|---|
| Envelope FFT | Spectrum + vertical lines for the four faults | + grey noise floor, red threshold curve, red shading where the spectrum rises above it, every line above the threshold marked and named (e.g. `2xFTF 19.0 Hz +15 dB`), optional Healthy recording in green |
| Why the verdict is what it is | Banner + one line of numbers | + bar chart of the evidence per fault, text box that states the rule and the numbers for the current window |
| VMD | Only on the plain signal (window 2, K = 7) | + new third window: VMD of the COMBINED signal, K = 18, modes named by fault line, fault rebuilt from its modes, alignment r, RMS, energy share |
| Entropy features | Window 2: raw signal + 7 IMFs | + 7 entropies for every mode of the combined signal, for each rebuilt fault and for the combined signal; fault-feature table |
| Excel log | `--log` (window 2 only) | + `--bearing-log` (window 3: evidence, alignment, modes, all entropies) |
| Background jobs | One entropy job | Two jobs (entropy + VMD) held as one `_JobPair`; one pair in flight at a time |
| New command-line options | | `--vmd-source combined/envelope`, `--comb-k`, `--bearing-log`; `--reference` now also draws the Healthy spectrum |
| Untouched | | 47 of 49 base functions and classes are byte-for-byte identical; only `main` and `make_bearing_dataset_simulator` changed. The old `BearingDatasetSimulator` class still exists |

**New functions:** `reference_envelope_spectrum`, `noise_floor_curve`, `strongest_lines`, `analyze_combined_vmd`, `verdict_from`, `combined_result_to_row`, `_JobPair`, `BearingVMDSimulator`.

**VMD settings and why**

| Setting | Value | Reason (measured on the real recordings) |
|---|---|---|
| Signal given to VMD | Combined = BPFO + BPFI + BSF + FTF (16 lines) | The user's request |
| K | 18 (16 lines + 2 spare) | Over 120 real 10 s windows: K=16 failed in 5 (worst r = 0.58), K=17 in 1, K=18 in 0 (worst r = 0.953), K=19 in 2, K=20 in 1 |
| Bandwidth | alpha = (fe / 3 Hz)^2, fe = 1000 Hz | 1.5 Hz: 12-16 % error. 6 Hz and 10 Hz: median r falls to 0.89-0.90 |
| tau, tol, max iterations | 0, 1e-7, 500 | Noise-tolerant default of the algorithm |
| Labelling | Largest peak of the mode's own spectrum, must be >= 25 % of its strongest peak | Same rule as the older `label_modes` |

**What happens with only 4 modes (K = 4)** (10 windows each, Healthy and Faulty1)

| | K = 4 | K = 18 |
|---|---|---|
| Signal lost (residual) | 47 to 67 % | 3 to 4 % |
| Faults that get a mode | 2.5 to 2.8 of 4 | 4 of 4 |
| Alignment r, Faulty1 (BPFO / BPFI / BSF / FTF) | 0.43 / 0.87 / 0.00 / 0.90 | about 1.0 each |

Reason: VMD separates by frequency, and each fault is a family of 4 lines, so 4 modes cannot hold 4 faults.

**How close is VMD to the directly filtered components** (Faulty1, 15 windows, K = 18)

| Fault | Size of VMD result vs filtered | Difference |
|---|---|---|
| FTF | 99 % | 2 % |
| BPFO | 98 % | 6 % |
| BSF | 97 % | 8 % |
| BPFI | 96 % | 8 % |

BSF and BPFI differ most because their lines lie only 3 to 4 Hz apart (for example 97.4 Hz and 101.0 Hz), so neighbouring modes overlap.

**Diagnosis results with the ball count set to 8** (30 windows of 10 s per recording, threshold 8 dB)

| Recording | Verdict | Mean evidence BPFO / BPFI / BSF / FTF (dB) |
|---|---|---|
| Healthy | HEALTHY 30 of 30 | 4.0 / 4.1 / 4.0 / 4.1 |
| Healthy2 | HEALTHY 30 of 30 | 4.0 / 4.2 / 4.0 / 4.1 |
| Faulty1 | FTF FAULTY 30 of 30 | 3.6 / 8.7 / 4.9 / 12.6 |
| Faulty2 | FTF FAULTY 29, BPFI FAULTY 1 of 30 | 3.4 / 8.4 / 4.7 / 12.6 |

**Timing:** one VMD job takes about 8 to 10 s per 10 s window (entropies of 21 series run on 4 threads). At real-time speed it just keeps up; at high `--speed` some windows are skipped. Drawing the main window takes about 0.05 s per frame.

**Known limits (state these in the report)**
1. The combined signal is built from the same fault lines, so VMD recovering the four parts shows the decomposition works; it does not prove a fault exists (r is about 0.99 on Healthy too). The proof is the evidence in dB against the noise floor.
2. FTF = 9.52 Hz was measured from the Faulty recordings, so the FTF result is consistent with a cage fault but is not independent proof. Confirm with the bearing model number and a physical inspection.
3. The 8 dB threshold comes from only two Healthy recordings (highest value 7.1 dB, +1 dB margin).
4. BPFI sits just above the threshold in the Faulty files (about 8.4 to 8.7 dB) and won the banner in 1 of 30 Faulty2 windows. Likely a coincidence, not shown to be one.
5. Unexplained envelope lines in the Faulty files: about 5.75, 24.8, 49.35 and 98.7 Hz. 24.8 Hz is the shaft speed, so imbalance or looseness may also be present.
6. A narrow-band mode is close to a sine wave, so its entropy is low and mostly reflects its frequency. Compare the same mode between Healthy and Faulty, not different modes with each other.

---

### Template for the next entry (copy, fill in, put above the previous entry)

```
### <date>: <short title>

**Why:** <what was asked or what went wrong>

**Comparison with the previous version**
| Item | Before | After |
|---|---|---|

**Settings changed (old -> new, and the measurement that justified it)**

**Test results (what was run, on which recordings, the numbers)**

**Known limits / open questions**
```

---

## Quick reference: current settings

| Setting | Value | Where |
|---|---|---|
| Sampling rate | 2000 Hz | dataset |
| Shaft speed | 24.8 Hz (estimated, 1488 rpm) | `estimate_shaft_speed` |
| Cage frequency FTF | 9.52 Hz (measured) | `DATA_FTF_HZ` |
| Balls | 8 (confirmed by the user, pass `--balls 8`) | `DATA_BALLS` |
| Resonance band | 600 to 900 Hz | `DATA_BAND` |
| Envelope rate | 1000 Hz | `DATA_ENV_FS` |
| Fault line half-width (filtered components) | 0.8 Hz | `DATA_LINE_HW_HZ` |
| Fault threshold | 8 dB (or Healthy maximum + 1 dB with `--reference`) | `DATA_FAULT_DB` |
| Entropy cut-off | Healthy mean -/+ 3 std (from `--reference`, else built-in) | `ENT_Z_CLEAR`, `ENT_DEFAULT_STATS` |
| Entropy verdict | FAULTY if Dispersion AND Fuzzy both below the cut-off | `ENT_VERDICT_KEYS` |
| Window for diagnosis | 10 s | `DATA_SPEC_SECONDS` |
| VMD (combined) | K = 18, half-bandwidth 3 Hz | `COMB_VMD_K`, `COMB_VMD_HALFBW` |
| VMD (plain mode only) | K = 7 | `VMD_K` |
| Entropy settings | m = 2, r = 0.2 x std, permutation order 5, dispersion c = 6, 32 histogram bins, MSE scales 1 to 10 | top of the file |

## Run command

```
python plot_dataset_time_fft_entropy_bearing_23nov26.py --file Vibration_Dataset_Faulty1_2000Hz.xlsx --bearing --balls 8 --label Faulty --reference Vibration_Dataset_Healthy_2000Hz.xlsx --log entropy_faulty.xlsx --bearing-log vmd_faulty.xlsx
```
