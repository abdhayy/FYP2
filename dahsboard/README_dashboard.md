# Bearing fault dashboard (React)

A web dashboard that shows the output of the Python analysis
(`plot_dataset_time_fft_entropy_bearing_23nov26.py`) as the same three windows:

| Tab | Same as | Shows |
|---|---|---|
| 1 · Live monitor | matplotlib window 1 | time domain, FFT with Mode 1/2/3, four fault components and their sum, envelope FFT with noise floor / threshold / named lines, evidence bars, verdict rule |
| 2 · Entropy | window 2 | 7 entropies vs the Healthy range (gauges), z-score bars, trend, multiscale entropy, ENTROPY VERDICT |
| 3 · VMD | window 3 | VMD of the combined signal, rebuilt faults, mode table, fault-feature table |
| Overview | – | every recording: how many windows each method calls FAULTY |

Play / pause / restart and the slider step through the 10 s windows like the Python replay
(keys: SPACE, ← →, R).

## How it fits together

```
Excel recordings ──► plot_dataset_time_fft_entropy_bearing_<date>.py   (YOUR analysis code)
                              │  imported by
                              ▼
                     export_dashboard_data.py  ──►  dashboard-app/public/data/*.json
                                                              │
                                                              ▼
                                                    dashboard-app (React)  ──►  browser
```

`export_dashboard_data.py` has no analysis of its own. It imports the newest
`plot_dataset_time_fft_entropy_bearing_*.py` in this folder (newest = date in the file name)
and calls its functions. **When you change the Python code, run the exporter again and the
dashboard shows the new results.** Only if you add a completely new plot or value do the
exporter and a React file need a new entry.

## Update the results after changing the Python code

```bash
cd dahsboard
pip install numpy scipy pandas openpyxl matplotlib
python export_dashboard_data.py
```

Options: `--script <file.py>` (use a specific version), `--files a.xlsx b.xlsx`,
`--reference Healthy.xlsx` (default: the Healthy file), `--balls 8`, `--max-windows 3` (quick test).
All four recordings take about 30 s.

## Run the dashboard

Needs Node.js 18 or newer.

```bash
cd dahsboard/dashboard-app
npm install
npm run dev
```

Open the address it prints (http://localhost:5173). `npm run build` makes a static site in
`dashboard-app/dist/` that can be opened from any web server (for example GitHub Pages).

## Files

| File | Purpose |
|---|---|
| `export_dashboard_data.py` | runs the analysis, writes the JSON |
| `dashboard-app/public/data/` | exported results: `index.json`, one folder per recording (`summary.json`, `w01.json` … `w30.json`) |
| `dashboard-app/src/App.jsx` | page layout, recording picker, playback |
| `dashboard-app/src/Window1.jsx` | tab 1 (live monitor) |
| `dashboard-app/src/Window2.jsx` | tab 2 (entropy) |
| `dashboard-app/src/Window3.jsx` | tab 3 (VMD) |
| `dashboard-app/src/Overview.jsx` | overview table |
| `dashboard-app/src/common.js` | colours, fault-line lists, helpers |

## MaFaulDa dataset tab (public data, UFRJ)

`mafaulda_analysis.py` runs the same method (envelope evidence + 7 entropies, imported from your
analysis script) on the MaFaulDa CSV files. Run it on the computer that has the dataset:

```bash
cd dahsboard
python mafaulda_analysis.py --root "C:\Users\User\Desktop\python fyp"
```

Quick test first: add `--limit 5` (5 files per class). It writes
`mafaulda_results/mafaulda_features.csv`, `mafaulda_results/mafaulda_summary.json` and
`dashboard-app/public/data/mafaulda.json` (read by the "MaFaulDa dataset" tab). Stopped runs continue
where they stopped. `--evaluate-only` redoes only the evaluation.

The raw CSV files must NOT be committed (they are ignored by `.gitignore`); only the small result files.
