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
def shaft_speed(tach):
    th = (tach.max() + tach.min()) / 2
    up = np.where((tach[:-1] < th) & (tach[1:] >= th))[0]
    if len(up) < 3:
        return float("nan")
    return FS / float(np.median(np.diff(up)))


def features(path):
    x = pd.read_csv(path, header=None).to_numpy(dtype=float)
    fr = shaft_speed(x[:, 0])
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


def main():
    global M, RESULTS, DASH_OUT
    ap = argparse.ArgumentParser(description="Run the FYP method on the MaFaulDa dataset.")
    ap.add_argument("--root", help="folder that contains normal/, underhang/, overhang/")
    ap.add_argument("--script", default=None, help="analysis script to import (default: newest *_bearing_*.py here)")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--limit", type=int, default=None, help="at most N files per class and bearing (quick test)")
    ap.add_argument("--evaluate-only", action="store_true", help="skip step 1, use the saved features")
    ap.add_argument("--results", default=RESULTS, help="folder for the features CSV and summary JSON")
    ap.add_argument("--dashboard-out", default=DASH_OUT, help="JSON file read by the dashboard ('' = do not write)")
    args = ap.parse_args()

    RESULTS, DASH_OUT = os.path.abspath(args.results), args.dashboard_out
    script = os.path.abspath(args.script) if args.script else find_script()
    M = load_module(script)
    print(f"Analysis script: {os.path.basename(script)}")
    if args.evaluate_only:
        df = pd.read_csv(os.path.join(RESULTS, "mafaulda_features.csv"))
    else:
        if not args.root:
            raise SystemExit("Give --root (the folder that contains normal/, underhang/, overhang/).")
        df = run_features(os.path.abspath(args.root), script, args.workers, args.limit)
    evaluate(df, os.path.basename(script))


if __name__ == "__main__":
    main()
