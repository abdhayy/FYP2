"""
MaFaulDa analysis: test the FYP method (envelope evidence + entropy) on the public
Machinery Fault Database (UFRJ). Run this on YOUR laptop, where the dataset is.

It does NOT have its own copy of the signal processing. It imports your newest analysis
script (plot_dataset_time_fft_entropy_bearing_<date>.py) and uses its functions:
demodulate, welch_amplitude_spectrum, diagnose (evidence in dB) and entropy_of_window
(the 7 entropies). Change those in your script, run this again, and the results follow.

STEP 1  features   every CSV file -> shaft speed (tachometer), fault evidence (dB) for
                   every accelerometer x frequency band, and the 7 entropies per accelerometer.
                   Saved to mafaulda_results/mafaulda_features.csv. Already processed files
                   are skipped, so you can stop and continue later.
STEP 2  evaluate   healthy range from the NORMAL files of similar speed -> z-scores ->
                   choose the best sensor/band on one half of the files ("tune") ->
                   report detection on the OTHER half ("test"), so the result is not
                   measured on the same files that chose the settings.
                   Saved to mafaulda_results/mafaulda_summary.json and to the dashboard
                   (dashboard-app/public/data/mafaulda.json).

Usage (Windows PowerShell, from the dahsboard folder):
  python mafaulda_analysis.py --root "C:\\Users\\User\\Desktop\\python fyp"
  python mafaulda_analysis.py --root "..." --limit 20         # quick test: 20 files per class
  python mafaulda_analysis.py --evaluate-only                  # only redo step 2

Expected folders under --root (nested copies like underhang/underhang are fine):
  normal/*.csv
  underhang/<cage_fault|outer_race|ball_fault>/<0g|6g|20g|35g>/*.csv
  overhang/<cage_fault|outer_race|ball_fault>/<0g|6g|20g|35g>/*.csv
Other folders (imbalance, misalignment) are ignored.

CSV columns (MaFaulDa, 50 kHz, 5 s, no header) - CHECK against the dataset description:
  1 tachometer | 2-4 underhang accelerometer (axial, radial, tangential)
  5-7 overhang accelerometer (axial, radial, tangential) | 8 microphone
"""
import argparse
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from export_dashboard_data import find_script, load_module  # noqa: E402

FS = 50000
ACC_COLS = {2: "UH axial", 3: "UH radial", 4: "UH tangential", 5: "OH axial", 6: "OH radial", 7: "OH tangential"}
BANDS = [(600, 900), (1000, 3000), (3000, 6000), (6000, 12000), (12000, 20000)]
FAULT_FACTORS = {"FTF": 0.3750, "BSF": 1.8710, "BPFO": 2.9980, "BPFI": 5.0020}   # x shaft speed, MaFaulDa bearing
CLASS_LINE = {"cage_fault": "FTF", "outer_race": "BPFO", "ball_fault": "BSF"}   # the fault family each class should show
CLASSES = ["normal", "cage_fault", "outer_race", "ball_fault"]
ENT_FS = 2000          # accelerometer decimated to this rate before the entropy step (then to 500 Hz inside it)
Z_CLEAR = 3.0          # same rule as the FYP: |z| >= 3 = clear change
K_NEIGHBOURS = 10      # normal files of the closest speed used as the healthy reference
RESULTS = os.path.join(HERE, "mafaulda_results")
DASH_OUT = os.path.join(HERE, "dashboard-app", "public", "data", "mafaulda.json")

M = None


def _init(script):
    global M
    M = load_module(script)


# ------------------------------------------------------------------ file list
def describe(path, root):
    parts = [p.lower() for p in os.path.relpath(path, root).replace("\\", "/").split("/")]
    if "normal" in parts:
        return {"cls": "normal", "bearing": "-", "weight": "-"}
    bearing = "underhang" if "underhang" in parts else "overhang" if "overhang" in parts else None
    cls = next((c for c in CLASS_LINE if c in parts), None)
    weight = next((p for p in parts if p.endswith("g") and p[:-1].isdigit()), "-")
    if bearing and cls:
        return {"cls": cls, "bearing": bearing, "weight": weight}
    return None


def list_files(root, limit=None):
    out = []
    for dp, _, fns in os.walk(root):
        for fn in fns:
            if fn.lower().endswith(".csv"):
                p = os.path.join(dp, fn)
                d = describe(p, root)
                if d:
                    out.append((p, d))
    out.sort(key=lambda t: t[0])
    if limit:
        keep, count = [], {}
        for p, d in out:
            key = (d["cls"], d["bearing"])
            if count.get(key, 0) < limit:
                keep.append((p, d))
                count[key] = count.get(key, 0) + 1
        out = keep
    return out


# ------------------------------------------------------------------ step 1: features of one file
def name_speed(path):
    """MaFaulDa file names are the nominal rotation frequency in Hz (e.g. 24.3712.csv)."""
    try:
        return float(os.path.splitext(os.path.basename(path))[0])
    except ValueError:
        return float("nan")


def speed_ok(fr, path):
    hint = name_speed(path)
    return np.isfinite(fr) and 8.0 <= fr <= 70.0 and (not np.isfinite(hint) or abs(fr - hint) <= 0.15 * hint)


def shaft_speed(tach, path=""):
    """Shaft speed from the tachometer (one pulse per turn). A Schmitt trigger (rise above 75 %, re-arm below 25 %
    of the range) ignores noise spikes that would count as extra pulses. If the result is not plausible
    (outside 8-70 Hz or more than 15 % from the speed in the file name) the speed in the file name is used."""
    lo, hi = np.percentile(tach, [5, 95])
    a, b = lo + 0.25 * (hi - lo), lo + 0.75 * (hi - lo)
    edges, armed = [], True
    for i in np.flatnonzero((tach < a) | (tach > b)):
        if armed and tach[i] > b:
            edges.append(i)
            armed = False
        elif not armed and tach[i] < a:
            armed = True
    fr = FS / float(np.median(np.diff(edges))) if len(edges) >= 3 else float("nan")
    return fr if speed_ok(fr, path) else name_speed(path)


def features(path):
    x = pd.read_csv(path, header=None).to_numpy(dtype=float)
    fr = shaft_speed(x[:, 0], path)
    F = {"fr": fr, **{k: v * fr for k, v in FAULT_FACTORS.items()}}
    row = {"fr": round(fr, 4)}
    for c in ACC_COLS:
        s = x[:, c - 1]
        for bi, band in enumerate(BANDS):
            env, fe = M.demodulate(s, FS, band, fe_target=M.DATA_ENV_FS, lowpass=None)
            fq, ae = M.welch_amplitude_spectrum(env, fe, 3)
            for k, v in M.diagnose(fq, ae, F).items():
                row[f"c{c}_b{bi}_{k}"] = round(v, 3)
        d = M.decimate_stages(s - s.mean(), FS // ENT_FS)
        vals, _ = M.entropy_of_window(d, ENT_FS, M.ANALYSIS_FS)
        for k, v in vals.items():
            row[f"c{c}_{k}"] = round(float(v), 5)
    return row


def _job(item):
    path, d = item
    t0 = time.time()
    try:
        return {"path": path, **d, **features(path), "error": "", "sec": round(time.time() - t0, 1)}
    except Exception as exc:                     # a broken file must not stop the whole run
        return {"path": path, **d, "error": str(exc), "sec": round(time.time() - t0, 1)}


def run_features(root, script, workers, limit):
    os.makedirs(RESULTS, exist_ok=True)
    out_csv = os.path.join(RESULTS, "mafaulda_features.csv")
    files = list_files(root, limit)
    if not files:
        raise SystemExit(f"No MaFaulDa CSV files found under {root}")
    done = pd.read_csv(out_csv) if os.path.exists(out_csv) else pd.DataFrame()
    have = set(done["path"]) if len(done) else set()
    if len(done):                                   # files saved with an impossible shaft speed (old tachometer bug): redo them
        bad = ~np.array([speed_ok(f, p) for f, p in zip(done.fr, done.path)])
        if bad.any():
            print(f"{int(bad.sum())} saved file(s) had an impossible shaft speed and will be processed again")
            done = done[~bad].reset_index(drop=True)
            have = set(done["path"])
    todo = [f for f in files if f[0] not in have]
    print(f"{len(files)} files found ({', '.join(f'{c}: {sum(d['cls'] == c for _, d in files)}' for c in CLASSES)}), "
          f"{len(todo)} still to do, {workers} worker(s)")
    rows, t0 = [], time.time()
    with ProcessPoolExecutor(max_workers=workers, initializer=_init, initargs=(script,)) as ex:
        futs = [ex.submit(_job, f) for f in todo]
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            rows.append(r)
            if r["error"]:
                print(f"  ! {r['path']}: {r['error']}")
            if i % 10 == 0 or i == len(todo):
                eta = (time.time() - t0) / i * (len(todo) - i)
                print(f"  {i}/{len(todo)} files, about {eta / 60:.0f} min left")
                pd.concat([done, pd.DataFrame(rows)], ignore_index=True).to_csv(out_csv, index=False)
    df = pd.concat([done, pd.DataFrame(rows)], ignore_index=True)
    df.to_csv(out_csv, index=False)
    return df


# ------------------------------------------------------------------ step 2: evaluate
def zscores(df, cols):
    """z of every feature against the K normal files with the closest shaft speed (a normal file never uses itself)."""
    nrm = df[df.cls == "normal"]
    Z = np.full((len(df), len(cols)), np.nan)
    for i, (idx, r) in enumerate(df.iterrows()):
        ref = nrm[nrm.index != idx]
        ref = ref.iloc[np.argsort(np.abs(ref.fr.to_numpy() - r.fr))[:K_NEIGHBOURS]]
        mu = ref[cols].mean().to_numpy()
        sd = np.maximum(ref[cols].std(ddof=1).to_numpy(), 1e-9)
        Z[i] = (r[cols].to_numpy(dtype=float) - mu) / sd
    return pd.DataFrame(Z, index=df.index, columns=cols)


def split_of(path):
    """Fixed half/half split from the file name, so it is the same on every run."""
    return "tune" if int(hashlib.md5(os.path.basename(os.path.dirname(path) + os.path.basename(path)).encode())
                         .hexdigest(), 16) % 2 == 0 else "test"


def score(z, cls, target, positive_if):
    """Youden-style score: hit rate on files of `target` class - false alarms on normal files."""
    hit = positive_if(z[cls == target]).mean() if (cls == target).any() else 0
    fa = positive_if(z[cls == "normal"]).mean() if (cls == "normal").any() else 0
    return hit - fa, hit, fa


def evaluate(df, script_name):
    df = df[(df.get("error", "").fillna("") == "") & df.fr.notna()].copy()
    df = df.reset_index(drop=True)
    n_norm = (df.cls == "normal").sum()
    if n_norm < 4:
        raise SystemExit(f"Only {n_norm} normal files: need at least 4 (better all ~49) to build the healthy range.")
    df["split"] = df.path.map(split_of)
    feat_cols = [c for c in df.columns if c.startswith("c") and "_" in c and c[1].isdigit()]
    Z = zscores(df, feat_cols)
    tune, test = df.split == "tune", df.split == "test"
    cls = df.cls.to_numpy()

    # ---- envelope: for every fault class choose the accelerometer x band whose fault evidence rises most
    env_choice = {}
    for target, line in CLASS_LINE.items():
        best = None
        for c in ACC_COLS:
            for bi in range(len(BANDS)):
                col = f"c{c}_b{bi}_{line}"
                s, hit, fa = score(Z.loc[tune, col].to_numpy(), cls[tune], target, lambda v: v >= Z_CLEAR)
                if best is None or s > best[0]:
                    best = (s, col, c, bi, hit, fa)
        env_choice[target] = {"feature": best[1], "column": best[2], "sensor": ACC_COLS[best[2]],
                              "band": list(BANDS[best[3]]), "line": line,
                              "tune_hit": round(best[4], 3), "tune_false_alarm": round(best[5], 3)}

    # ---- entropy: choose the accelerometer x entropy that moves most (either direction) for ANY bearing fault
    ent_keys = list(M.ENT_KEYS)
    any_fault = np.where(cls == "normal", "normal", "fault")
    best = None
    for c in ACC_COLS:
        for k in ent_keys:
            col = f"c{c}_{k}"
            s, hit, fa = score(Z.loc[tune, col].to_numpy(), any_fault[tune], "fault", lambda v: np.abs(v) >= Z_CLEAR)
            if best is None or s > best[0]:
                best = (s, col, c, k, hit, fa)
    ent_choice = {"feature": best[1], "column": best[2], "sensor": ACC_COLS[best[2]], "entropy": best[3],
                  "name": M.ENT_NAMES[best[3]], "tune_hit": round(best[4], 3), "tune_false_alarm": round(best[5], 3)}

    # ---- verdicts for every file
    zf = {t: Z[ch["feature"]].to_numpy() for t, ch in env_choice.items()}
    stack = np.vstack([zf[t] for t in CLASS_LINE])
    env_fault = (stack >= Z_CLEAR).any(axis=0)
    env_type = np.array(list(CLASS_LINE))[np.argmax(stack, axis=0)]
    df["env_verdict"] = np.where(env_fault, env_type, "normal")
    df["ent_z"] = Z[ent_choice["feature"]].to_numpy()
    df["ent_verdict"] = np.where(np.abs(df.ent_z) >= Z_CLEAR, "fault", "normal")
    for t in CLASS_LINE:
        df[f"z_{t}"] = zf[t]

    def rates(sub):
        out = {}
        for c in CLASSES:
            s = sub[sub.cls == c]
            if not len(s):
                continue
            det = (s.env_verdict != "normal") if c != "normal" else (s.env_verdict == "normal")
            typ = (s.env_verdict == c)
            ent = (s.ent_verdict == "fault") if c != "normal" else (s.ent_verdict == "normal")
            out[c] = {"n": int(len(s)), "env_correct": round(float(det.mean()), 3),
                      "env_type_correct": round(float(typ.mean()), 3), "ent_correct": round(float(ent.mean()), 3)}
        return out

    T = df[test]
    confusion = {a: {b: int(((T.cls == a) & (T.env_verdict == b)).sum()) for b in CLASSES} for a in CLASSES}
    by_group = []
    for (b, w, c), s in T[T.cls != "normal"].groupby(["bearing", "weight", "cls"]):
        by_group.append({"bearing": b, "weight": w, "cls": c, "n": int(len(s)),
                         "env_detected": round(float((s.env_verdict != "normal").mean()), 3),
                         "env_type_correct": round(float((s.env_verdict == c).mean()), 3),
                         "ent_detected": round(float((s.ent_verdict == "fault").mean()), 3)})

    summary = {
        "generated": time.strftime("%Y-%m-%d %H:%M"), "script": script_name,
        "files": {c: int((df.cls == c).sum()) for c in CLASSES},
        "settings": {"fs": FS, "bands": BANDS, "fault_factors": FAULT_FACTORS, "z_clear": Z_CLEAR,
                     "k_neighbours": K_NEIGHBOURS, "entropy_rate": M.ANALYSIS_FS,
                     "columns": {str(k): v for k, v in ACC_COLS.items()}},
        "envelope_choice": env_choice, "entropy_choice": ent_choice,
        "test": rates(T), "tune": rates(df[tune]), "confusion_test": confusion, "by_group_test": by_group,
        "files_detail": [{"file": os.path.relpath(r.path, os.path.commonpath(list(df.path))).replace("\\", "/"),
                          "cls": r.cls, "bearing": r.bearing, "weight": r.weight, "fr": round(r.fr, 2),
                          "split": r.split, "env": r.env_verdict, "ent": r.ent_verdict,
                          "z": {t: (None if not np.isfinite(getattr(r, f"z_{t}")) else round(float(getattr(r, f"z_{t}")), 2)) for t in CLASS_LINE},
                          "entZ": None if not np.isfinite(r.ent_z) else round(float(r.ent_z), 2)}
                         for r in df.itertuples()],
    }
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "mafaulda_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    if DASH_OUT:
        os.makedirs(os.path.dirname(os.path.abspath(DASH_OUT)), exist_ok=True)
        with open(DASH_OUT, "w") as fh:
            json.dump(summary, fh, separators=(",", ":"))

    # ---- print a short report
    print("\n=== MaFaulDa result (TEST half: files not used to choose the settings) ===")
    for t, ch in env_choice.items():
        print(f"  envelope {t:10s}: {ch['sensor']:13s} band {ch['band'][0]}-{ch['band'][1]} Hz, line {ch['line']}")
    print(f"  entropy          : {ent_choice['sensor']}, {ent_choice['name']}")
    print(f"  {'class':11s} {'n':>4s}  {'env correct':>11s}  {'env type ok':>11s}  {'entropy ok':>10s}")
    for c, v in summary["test"].items():
        print(f"  {c:11s} {v['n']:4d}  {v['env_correct']:11.0%}  {v['env_type_correct']:11.0%}  {v['ent_correct']:10.0%}")
    print(f"\nSaved {os.path.join(RESULTS, 'mafaulda_summary.json')}" + (f" and {DASH_OUT}" if DASH_OUT else ""))
    return summary


# ------------------------------------------------------------------ step 3: per-file views (windows 1, 2, 3)
VIEW_DIR = os.path.join(HERE, "dashboard-app", "public", "data", "mafaulda_views")


def _acc(path, col):
    return pd.read_csv(path, header=None).to_numpy(dtype=float)[:, [0, col - 1]]


def _view_job(args):
    """Windows 1-3 for ONE MaFaulDa file, with the healthy reference built from its speed-matched normal files."""
    import export_dashboard_data as X
    path, info, env_col, band, ent_col, ref_paths = args
    m = M
    X.CTX["group_len"] = 5 * FS
    raw = _acc(path, env_col)
    fr = shaft_speed(raw[:, 0], path)
    F = {"fr": fr, **{k: v * fr for k, v in FAULT_FACTORS.items()}}
    known = {k: True for k in FAULT_FACTORS}
    x = raw[:, 1] - np.median(raw[:, 1])
    n = len(x)

    # healthy reference of THIS file: its normal neighbours, same sensor and band
    ref_ev, ref_spec, ref_ent = [], None, []
    for rp in ref_paths:
        r = _acc(rp, env_col)
        r_fr = shaft_speed(r[:, 0], rp)
        rF = {"fr": r_fr, **{k: v * r_fr for k, v in FAULT_FACTORS.items()}}
        env, fe = m.demodulate(r[:, 1], FS, band, fe_target=m.DATA_ENV_FS, lowpass=None)
        fq, ae = m.welch_amplitude_spectrum(env, fe, 3)
        ref_ev.append(max(m.diagnose(fq, ae, rF).values()))
        ref_spec = ae if ref_spec is None else ref_spec + ae
        e = _acc(rp, ent_col)[:, 1]
        d = m.decimate_stages(e - e.mean(), FS // ENT_FS)
        ref_ent.append(m.entropy_of_window(d, ENT_FS, m.ANALYSIS_FS))
    fault_db = float(max(ref_ev)) + 1.0                        # same rule as calibrate_threshold: healthy maximum + 1 dB
    keys = list(m.ENT_KEYS)
    mses = np.array([v[1] for v in ref_ent])
    stats = {"mean": {k: float(np.nanmean([v[0][k] for v in ref_ent])) for k in keys},
             "std": {k: max(float(np.nanstd([v[0][k] for v in ref_ent], ddof=1)), 1e-9) for k in keys},
             "mse_mean": np.nanmean(mses, axis=0), "mse_std": np.maximum(np.nanstd(mses, axis=0, ddof=1), 1e-9),
             "source": f"{len(ref_paths)} MaFaulDa normal files with the closest shaft speed"}

    w1 = X.window1(m, x, FS, n, F, known, band, fault_db, 300.0)
    d = m.decimate_stages(x - x.mean(), FS // ENT_FS)          # time display from the 2 kHz copy (the 50 kHz trace is too long)
    t = np.arange(1, len(d) + 1) / ENT_FS
    w1["time"] = {"t": X.arr(t[::4], 5), "raw": X.arr(d[::4], 3),
                  "smooth": X.arr(m.smooth_for_display(d, ENT_FS, m.SMOOTH_CUTOFF)[::4], 3)}
    e_raw = _acc(path, ent_col)[:, 1]
    de = m.decimate_stages(e_raw - e_raw.mean(), FS // ENT_FS)
    e = m.analyze_entropy_window(de, ENT_FS, stats, m.ANALYSIS_FS, 0.0)
    w2 = {"vals": {k: X.r(v, 5) for k, v in e["vals"].items()}, "z": {k: X.r(v, 4) for k, v in e["z"].items()},
          "mse": X.arr(e["mse"], 5), "verdict": e["verdict"], "elapsed": X.r(e["elapsed"], 3)}
    v = m.analyze_combined_vmd(x, FS, F, known, band, fault_db, "combined", m.COMB_VMD_K, m.COMB_VMD_HALFBW, 0.0)
    w3 = X.window3(v, m)
    rf, ra = fq[fq <= 450], (ref_spec / len(ref_paths))[fq <= 450]
    label = NAMES_SHORT[info["cls"]]
    s = {"name": os.path.basename(path), "slug": "", "label": label, "fs": FS, "groups": 1, "groupSeconds": 5,
         "durationS": 5, "windowS": 5, "bearing": {k: X.r(val, 5) for k, val in F.items()}, "known": known, "balls": 8,
         "faultLines": {k: [{"name": nm, "f": X.r(f, 5)} for nm, f in vv] for k, vv in m.fault_lines(F).items()},
         "band": list(band), "faultDb": X.r(fault_db, 5), "envFs": m.DATA_ENV_FS, "reference": "speed-matched normal files",
         "refSpectrum": {"f": X.arr(rf), "a": X.arr(ra)},
         "entropy": {"keys": keys, "names": m.ENT_NAMES, "zClear": m.ENT_Z_CLEAR, "verdictKeys": list(m.ENT_VERDICT_KEYS),
                     "reliable": sorted(m.ENT_RELIABLE), "test": None,
                     "mean": {k: X.r(vv, 6) for k, vv in stats["mean"].items()},
                     "std": {k: X.r(vv, 6) for k, vv in stats["std"].items()},
                     "mseMean": X.arr(stats["mse_mean"], 5), "mseStd": X.arr(stats["mse_std"], 5),
                     "source": stats["source"], "analysisFs": m.ANALYSIS_FS},
         "vmd": {"K": m.COMB_VMD_K, "source": "combined", "halfBw": m.COMB_VMD_HALFBW},
         "windows": [{"index": 1, "tStart": 0, "tEnd": 5, "verdict": w1["verdict"], "evidence": w1["evidence"],
                      "entVerdict": w2["verdict"], "z": w2["z"], "vals": w2["vals"]}],
         "mafaulda": {**info, "fr": X.r(fr, 5), "envSensor": ACC_COLS[env_col], "entSensor": ACC_COLS[ent_col]}}
    w = {"index": 1, "tStart": 0, "tEnd": 5, "group": 1, "w1": w1, "w2": w2, "w3": w3}
    return {"s": s, "w": w}


NAMES_SHORT = {"normal": "Normal", "cage_fault": "Cage fault", "outer_race": "Outer race fault", "ball_fault": "Ball fault"}


def export_views(df, summary, workers, speed=25.0, extra=()):
    """Pick a few files (0g, test half, closest to `speed` Hz: one per fault class and bearing, two normal) and
    export the three windows for each, using the sensor/band the method chose for the fault it found."""
    df = df[(df.get("error", "").fillna("") == "")].reset_index(drop=True).copy()
    df["split"] = df.path.map(split_of)
    picks = []
    for (cls, bearing), g in df[(df.weight.isin(["0g", "-"]))].groupby(["cls", "bearing"]):
        g = g.iloc[np.argsort(np.abs(g.fr.to_numpy() - speed))]
        picks += list(g.path[: 2 if cls == "normal" else 1])
    picks += [p for p in extra if p not in picks]
    det = {f["file"]: f for f in summary["files_detail"]}
    common = os.path.commonpath(list(df.path))
    nrm = df[df.cls == "normal"]
    ent_col = summary["entropy_choice"]["column"]
    jobs = []
    for p in picks:
        row = df[df.path == p].iloc[0]
        fd = det.get(os.path.relpath(p, common).replace("\\", "/"), {})
        zs = fd.get("z", {}) or {}
        found = max(CLASS_LINE, key=lambda c: (zs.get(c) if zs.get(c) is not None else -1e9))
        ch = summary["envelope_choice"][found]
        ref = nrm[nrm.path != p]
        ref = list(ref.path.iloc[np.argsort(np.abs(ref.fr.to_numpy() - row.fr))[:K_NEIGHBOURS]])
        info = {"cls": row.cls, "bearing": row.bearing, "weight": row.weight, "file": os.path.relpath(p, common).replace("\\", "/"),
                "detector": found, "envVerdict": fd.get("env"), "entVerdict": fd.get("ent"), "split": fd.get("split")}
        jobs.append((p, info, ch["column"], tuple(ch["band"]), ent_col, ref))
    os.makedirs(VIEW_DIR, exist_ok=True)
    for old in os.listdir(VIEW_DIR):
        if old.endswith(".json"):
            os.remove(os.path.join(VIEW_DIR, old))
    print(f"\nExporting windows 1-3 for {len(jobs)} files (each uses {K_NEIGHBOURS} normal files as its healthy reference) ...")
    index = []
    script = summary["_script_path"]
    with ProcessPoolExecutor(max_workers=min(workers, len(jobs)), initializer=_init, initargs=(script,)) as ex:
        for i, out in enumerate(ex.map(_view_job, jobs), 1):
            vid = f"v{i:02d}"
            with open(os.path.join(VIEW_DIR, vid + ".json"), "w") as fh:
                json.dump(out, fh, separators=(",", ":"))
            mf = out["s"]["mafaulda"]
            index.append({"id": vid, "file": mf["file"], "cls": mf["cls"], "bearing": mf["bearing"], "weight": mf["weight"],
                          "fr": mf["fr"], "verdict": out["w"]["w1"]["verdict"], "entVerdict": out["w"]["w2"]["verdict"]})
            print(f"  {mf['file']:55s} spectrum: {out['w']['w1']['verdict']:<16s} entropy: {out['w']['w2']['verdict']}")
    with open(os.path.join(VIEW_DIR, "index.json"), "w") as fh:
        json.dump(index, fh, indent=1)


def main():
    global M, RESULTS, DASH_OUT, VIEW_DIR
    ap = argparse.ArgumentParser(description="Run the FYP method on the MaFaulDa dataset.")
    ap.add_argument("--root", help="folder that contains normal/, underhang/, overhang/")
    ap.add_argument("--script", default=None, help="analysis script to import (default: newest *_bearing_*.py here)")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--limit", type=int, default=None, help="at most N files per class and bearing (quick test)")
    ap.add_argument("--evaluate-only", action="store_true", help="skip step 1, use the saved features")
    ap.add_argument("--no-views", action="store_true", help="do not export the windows 1-3 views of example files")
    ap.add_argument("--view-speed", type=float, default=25.0, help="example files closest to this shaft speed (Hz)")
    ap.add_argument("--view-file", action="append", default=[], help="also export this CSV file (repeatable)")
    ap.add_argument("--results", default=RESULTS, help="folder for the features CSV and summary JSON")
    ap.add_argument("--dashboard-out", default=DASH_OUT, help="JSON file read by the dashboard ('' = do not write)")
    args = ap.parse_args()

    RESULTS, DASH_OUT = os.path.abspath(args.results), args.dashboard_out
    if DASH_OUT:
        VIEW_DIR = os.path.join(os.path.dirname(os.path.abspath(DASH_OUT)), "mafaulda_views")
    script = os.path.abspath(args.script) if args.script else find_script()
    M = load_module(script)
    print(f"Analysis script: {os.path.basename(script)}")
    if args.evaluate_only:
        df = pd.read_csv(os.path.join(RESULTS, "mafaulda_features.csv"))
    else:
        if not args.root:
            raise SystemExit("Give --root (the folder that contains normal/, underhang/, overhang/).")
        df = run_features(os.path.abspath(args.root), script, args.workers, args.limit)
    summary = evaluate(df, os.path.basename(script))
    if not args.no_views and DASH_OUT:
        summary["_script_path"] = script
        export_views(df, summary, args.workers, args.view_speed, [os.path.abspath(p) for p in args.view_file])


if __name__ == "__main__":
    main()
