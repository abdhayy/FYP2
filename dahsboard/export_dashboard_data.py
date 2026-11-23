"""
Export the results of the bearing analysis to JSON for the React dashboard (dashboard-app/).

This file does NOT contain its own copy of the analysis. It imports YOUR analysis script
(by default the newest  plot_dataset_time_fft_entropy_bearing_*.py  in this folder) and calls its
functions, so when you change the analysis code and run this exporter again, the dashboard shows
the new results.

For every recording and every 10 s window it saves what the three matplotlib windows show:
  window 1  live monitor      : time domain, FFT with Mode 1/2/3, four fault components, their sum,
                                envelope FFT with noise floor / threshold / named lines, evidence, verdict
  window 2  entropy           : 7 entropies, z-scores against the Healthy range, verdict, multiscale entropy
  window 3  VMD               : VMD of the combined signal, modes, rebuilt faults, fault-feature table

Usage (from this folder):
  python export_dashboard_data.py
  python export_dashboard_data.py --reference ../Vibration_Dataset_Healthy_2000Hz.xlsx --balls 8
  python export_dashboard_data.py --files ../Vibration_Dataset_Faulty1_2000Hz.xlsx --max-windows 5
  python export_dashboard_data.py --script plot_dataset_time_fft_entropy_bearing_23nov26.py

Output: dashboard-app/public/data/index.json + one folder per recording with summary.json and w01.json, w02.json ...
"""
import argparse
import glob
import importlib.util
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DEFAULT = os.path.join(HERE, "dashboard-app", "public", "data")
DEFAULT_FILES = ["Vibration_Dataset_Healthy_2000Hz.xlsx", "Vibration_Dataset_Healthy2_2000Hz.xlsx",
                 "Vibration_Dataset_Faulty1_2000Hz.xlsx", "Vibration_Dataset_Faulty2_2000Hz.xlsx"]
DEFAULT_REFERENCE = "Vibration_Dataset_Healthy_2000Hz.xlsx"


# ---------------------------------------------------------------- helpers
def find_script():
    """Newest analysis script in this folder (by the date in its name, then by modification time)."""
    c = [p for p in glob.glob(os.path.join(HERE, "plot_dataset_time_fft_entropy_bearing_*.py"))]
    if not c:
        raise SystemExit("No plot_dataset_time_fft_entropy_bearing_*.py found next to this exporter (use --script).")

    def key(p):
        name = os.path.basename(p)
        d = name.rsplit("_", 1)[-1].split(".")[0]           # e.g. 23nov26
        try:
            return (datetime.strptime(d, "%d%b%y"), os.path.getmtime(p))
        except ValueError:
            return (datetime.min, os.path.getmtime(p))
    return max(c, key=key)


def load_module(path):
    import matplotlib
    matplotlib.use("Agg")                                   # the script imports pyplot; no windows here
    spec = importlib.util.spec_from_file_location("analysis", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def find_file(name):
    for p in (name, os.path.join(HERE, name), os.path.join(HERE, "..", name)):
        if os.path.exists(p):
            return os.path.abspath(p)
    raise SystemExit(f"Cannot find {name} (looked here and in the parent folder).")


def r(x, sig=4):
    """Round to `sig` significant digits; NaN/inf -> None (valid JSON)."""
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(x):
        return None
    if x == 0:
        return 0.0
    return float(f"{x:.{sig}g}")


def arr(a, sig=4):
    return [r(v, sig) for v in np.asarray(a, dtype=float).ravel()]


def thin(a, n_max):
    """Keep at most n_max points (min/max pairs per bucket so peaks stay visible)."""
    a = np.asarray(a, dtype=float)
    if len(a) <= n_max:
        return a, np.arange(len(a))
    b = max(1, len(a) // (n_max // 2))
    m = len(a) // b * b
    blk = a[:m].reshape(-1, b)
    i_min, i_max = blk.argmin(1), blk.argmax(1)
    idx = np.sort(np.concatenate([np.arange(len(blk)) * b + np.minimum(i_min, i_max),
                                  np.arange(len(blk)) * b + np.maximum(i_min, i_max)]))
    return a[idx], idx


def series(t, y, n_max=1500, sig=4):
    y2, idx = thin(y, n_max)
    return {"t": arr(np.asarray(t)[idx], 5), "y": arr(y2, sig)}


# ---------------------------------------------------------------- one window
CTX = {}


def _init_worker(script, ctx):
    CTX["m"] = load_module(script)
    CTX.update(ctx)


def window1(m, stream, fs, n, F, known, band, fault_db, xmax):
    """Window 1 (live monitor) as it looks at sample n."""
    W = CTX["group_len"]
    lo = max(0, n - W)
    buf = np.zeros(W)
    buf[W - (n - lo):] = stream[lo:n]
    smooth = m.smooth_for_display(buf, fs, m.SMOOTH_CUTOFF)
    t_axis = np.arange(1, W + 1) / fs

    valid = stream[lo:n]
    freqs, p1 = m.welch_amplitude_spectrum(valid, fs, m.N_SEG_WELCH)
    modes = []
    f0 = m.estimate_fundamental_hps(freqs, p1)
    amp0 = p1[np.argmin(np.abs(freqs - f0))]
    detected = bool(amp0 > m.PEAK_RATIO_MIN * np.median(p1))
    if detected:
        for k in (1, 2, 3):
            fk = f0 * k
            if fk <= freqs[-1]:
                modes.append({"k": k, "f": r(fk), "a": r(p1[np.argmin(np.abs(freqs - fk))])})
    mf = freqs <= xmax

    have = min(n, int(m.DATA_SPEC_SECONDS * fs))
    seg = stream[n - have:n]
    env, fe = m.demodulate(seg, fs, band, fe_target=m.DATA_ENV_FS, lowpass=None)
    fq, ae = m.welch_amplitude_spectrum(env, fe, 3)
    floor = m.noise_floor_curve(fq, ae, fe)
    thr = floor * 10 ** (fault_db / 20.0)
    me = fq <= 450
    ev = m.diagnose(fq, ae, F)
    marks = []
    for fault in m.FAULT_ORDER:
        for name, f0_, pf, pa, snr in m.strongest_lines(fq, ae, F, fault, 2):
            marks.append({"fault": fault, "name": name, "f": r(pf), "a": r(pa), "db": r(snr, 3),
                          "above": bool(snr >= fault_db and pf <= 450)})
    verdict, top, topv = m.verdict_from(ev, known, fault_db)
    others = [k for k, v in ev.items() if known[k] and k != top and v >= fault_db] if verdict != "HEALTHY" else []

    comps = {}
    for fault, lines in m.fault_lines(F).items():
        s = np.zeros(len(env))
        for _, fl in lines:
            s += m.narrowband(env, fe, fl)
        comps[fault] = s
    total = sum(comps.values())
    cw = min(int(round(m.DATA_COMP_SECONDS * fe)), len(env))
    tc = np.arange(cw) / fe
    return {
        "time": {"t": arr(t_axis[::4], 5), "raw": arr(buf[::4], 3), "smooth": arr(smooth[::4], 3)},
        "fft": {"f": arr(freqs[mf]), "a": arr(p1[mf]), "modes": modes, "detected": detected},
        "env": {"f": arr(fq[me]), "a": arr(ae[me], 3), "floor": arr(floor[me], 3), "thr": arr(thr[me], 3), "marks": marks},
        "compT": arr(tc, 5),
        "components": {k: {"y": arr(comps[k][-cw:], 3), "peak": r(np.max(np.abs(comps[k])))} for k in m.FAULT_ORDER},
        "sum": {"env": arr(env[-cw:], 3), "sum": arr(total[-cw:], 3)},
        "evidence": {k: r(v, 4) for k, v in ev.items()},
        "verdict": verdict, "top": top, "topDb": r(topv, 4), "others": others,
        "topLines": [{"name": nm, "f": r(pf), "db": r(s, 3)} for nm, _, pf, _, s in m.strongest_lines(fq, ae, F, top, 2)],
    }


def _spectrum(f, a, n_max=1600):
    a2, idx = thin(a, n_max)
    return {"f": arr(f[idx]), "a": arr(a2, 3)}


def window3(res, m):
    fe = res["fe"]
    nt = int(m.COMB_TAIL_SECONDS * fe)
    t = np.arange(nt) / fe
    sig = res["sig"]
    modes_sum = res["modes"].sum(axis=0)
    sf = res["spec_f"] <= 330
    out = {
        "K": res["K"], "source": res["source"], "resid": r(res["resid"]), "nIter": int(res["n_iter"]),
        "elapsed": r(res["elapsed"], 3),
        "combined": {"t": arr(t[::2], 5), "sig": arr(sig[-nt:][::2], 3), "modesSum": arr(modes_sum[-nt:][::2], 3)},
        "faults": {},
        "spectrum": _spectrum(res["spec_f"][sf], res["spec_c"][sf] * 1000),
        "modes": [{"i": j + 1, "center": r(res["centers"][j]), "fault": res["labels"][j][0],
                   "line": res["labels"][j][1], "rms": r(res["rms"][j] * 1000)} for j in range(len(res["modes"]))],
    }
    for k in m.FAULT_ORDER:
        rec = res["recovered"][k]
        out["faults"][k] = {"orig": arr(res["comps"][k][-nt:][::2], 3), "rec": arr(rec[-nt:][::2], 3),
                            "imfs": [i + 1 for i in res["which"][k]], "r": r(res["align"][k]),
                            "rms": r(np.sqrt(np.mean(rec ** 2)) * 1000), "share": r(res["share"][k]),
                            "evidence": r(res["evidence"][k], 4)}
    out["verdict"] = m.verdict_from(res["evidence"], res["known"], res["fault_db"])[0]
    return out


def do_window(i):
    m = CTX["m"]
    stream, fs, L = CTX["stream"], CTX["fs"], CTX["L"]
    F, known, band, fault_db = CTX["F"], CTX["known"], CTX["band"], CTX["fault_db"]
    n = (i + 1) * L
    window = stream[n - L:n].copy()
    t_start = (n - L) / fs
    w1 = window1(m, stream, fs, n, F, known, band, fault_db, CTX["xmax"])
    e = m.analyze_entropy_window(window, fs, CTX["ent_stats"], CTX["analysis_fs"], t_start)
    w2 = {"vals": {k: r(v, 5) for k, v in e["vals"].items()}, "z": {k: r(v, 4) for k, v in e["z"].items()},
          "mse": arr(e["mse"], 5), "verdict": e["verdict"], "elapsed": r(e["elapsed"], 3)}
    v = m.analyze_combined_vmd(window, fs, F, known, band, fault_db, CTX["vmd_source"], CTX["comb_k"],
                               m.COMB_VMD_HALFBW, t_start)
    w3 = window3(v, m)
    return {"index": i + 1, "tStart": r(t_start, 6), "tEnd": r(n / fs, 6), "group": int((n - 1) // CTX["group_len"] + 1),
            "w1": w1, "w2": w2, "w3": w3}


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="Export the bearing analysis to JSON for the React dashboard.")
    ap.add_argument("--script", default=None, help="analysis script to import (default: newest *_bearing_*.py here)")
    ap.add_argument("--files", nargs="+", default=None, help="Excel recordings (default: the four in the repo)")
    ap.add_argument("--reference", default=DEFAULT_REFERENCE, help="Healthy recording for threshold + entropy range ('' = none)")
    ap.add_argument("--balls", type=int, default=8, help="ball count (default 8, as confirmed)")
    ap.add_argument("--fr", type=float, default=None)
    ap.add_argument("--ftf", type=float, default=None)
    ap.add_argument("--phi", type=float, default=0.0)
    ap.add_argument("--fault-db", type=float, default=None)
    ap.add_argument("--vmd-source", choices=["combined", "envelope"], default="combined")
    ap.add_argument("--comb-k", type=int, default=None)
    ap.add_argument("--max-windows", type=int, default=None, help="only the first N windows per recording (quick test)")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--out", default=OUT_DEFAULT)
    args = ap.parse_args()

    script = os.path.abspath(args.script) if args.script else find_script()
    print(f"Analysis script: {os.path.basename(script)}")
    m = load_module(script)
    files = [find_file(f) for f in (args.files or DEFAULT_FILES)]
    os.makedirs(args.out, exist_ok=True)

    ref = None
    if args.reference:
        _, fs_ref, _, ref = m.load_dataset(find_file(args.reference))
        print(f"Healthy reference: {os.path.basename(find_file(args.reference))}")

    recordings = []
    for path in files:
        name = os.path.splitext(os.path.basename(path))[0]
        slug = name.replace("Vibration_Dataset_", "").replace("_2000Hz", "")
        label = "Faulty" if "fault" in name.lower() else ("Healthy" if "health" in name.lower() else "Unlabelled")
        print(f"\n=== {name} ===")
        _, fs, _, data = m.load_dataset(path)
        stream_raw = data.T.reshape(-1)
        stream = stream_raw - np.median(stream_raw)
        fake = argparse.Namespace(fr=args.fr, ftf=args.ftf, balls=args.balls, ball_d=None, pitch_d=None, phi=args.phi)
        F, known, n_balls, how = m.resolve_bearing_params(fake, stream_raw, fs)
        band = tuple(m.DATA_BAND)
        fault_db = args.fault_db if args.fault_db is not None else m.DATA_FAULT_DB
        ref_spec, ent_stats = None, m.ENT_DEFAULT_STATS
        if ref is not None:
            fault_db = m.calibrate_threshold(ref, fs_ref, F, known, band) if args.fault_db is None else fault_db
            ref_spec = m.reference_envelope_spectrum(ref, fs_ref, band)
            ent_stats = m.calibrate_entropy_reference(ref, fs_ref, m.ANALYSIS_SECONDS, m.ANALYSIS_FS)
        comb_k = args.comb_k or (m.COMB_VMD_K if args.vmd_source == "combined" else m.ENV_VMD_K)
        L = int(round(m.ANALYSIS_SECONDS * fs))
        n_win = len(stream) // L if args.max_windows is None else min(args.max_windows, len(stream) // L)
        ctx = dict(stream=stream, fs=fs, L=L, F=F, known=known, band=band, fault_db=fault_db, ent_stats=ent_stats,
                   analysis_fs=m.ANALYSIS_FS, vmd_source=args.vmd_source, comb_k=comb_k, xmax=300.0,
                   group_len=data.shape[0])
        rec_dir = os.path.join(args.out, slug)
        os.makedirs(rec_dir, exist_ok=True)
        for old in glob.glob(os.path.join(rec_dir, "w*.json")):
            os.remove(old)

        t0 = time.time()
        summary_windows = []
        with ProcessPoolExecutor(max_workers=args.workers, initializer=_init_worker, initargs=(script, ctx)) as ex:
            for w in ex.map(do_window, range(n_win)):
                with open(os.path.join(rec_dir, f"w{w['index']:02d}.json"), "w") as fh:
                    json.dump(w, fh, separators=(",", ":"))
                summary_windows.append({"index": w["index"], "tStart": w["tStart"], "tEnd": w["tEnd"],
                                        "verdict": w["w1"]["verdict"], "evidence": w["w1"]["evidence"],
                                        "entVerdict": w["w2"]["verdict"], "z": w["w2"]["z"], "vals": w["w2"]["vals"]})
                print(f"  window {w['index']:2d}/{n_win}  {w['w1']['verdict']:<12}  entropy {w['w2']['verdict']}")
        print(f"  {n_win} windows in {time.time() - t0:.0f} s")

        summary = {
            "name": name, "slug": slug, "label": label, "fs": fs, "groups": int(data.shape[1]),
            "groupSeconds": data.shape[0] / fs, "durationS": len(stream) / fs, "windowS": m.ANALYSIS_SECONDS,
            "bearing": {k: r(v, 5) for k, v in F.items()}, "known": known, "balls": n_balls, "how": how,
            "faultLines": {k: [{"name": nm, "f": r(f, 5)} for nm, f in v] for k, v in m.fault_lines(F).items()},
            "band": list(band), "faultDb": r(fault_db, 5), "envFs": m.DATA_ENV_FS,
            "reference": os.path.basename(find_file(args.reference)) if args.reference else None,
            "refSpectrum": ({"f": arr(ref_spec[0][ref_spec[0] <= 450]), "a": arr(ref_spec[1][ref_spec[0] <= 450])}
                            if ref_spec is not None else None),
            "entropy": {"keys": list(m.ENT_KEYS), "names": m.ENT_NAMES, "zClear": m.ENT_Z_CLEAR,
                        "verdictKeys": list(m.ENT_VERDICT_KEYS), "reliable": sorted(m.ENT_RELIABLE),
                        "test": {k: list(v) for k, v in m.ENT_TEST.items()},
                        "mean": {k: r(v, 6) for k, v in ent_stats["mean"].items()},
                        "std": {k: r(v, 6) for k, v in ent_stats["std"].items()},
                        "mseMean": arr(ent_stats["mse_mean"], 5), "mseStd": arr(ent_stats["mse_std"], 5),
                        "source": ent_stats["source"], "analysisFs": m.ANALYSIS_FS},
            "vmd": {"K": comb_k, "source": args.vmd_source, "halfBw": m.COMB_VMD_HALFBW},
            "windows": summary_windows,
        }
        with open(os.path.join(rec_dir, "summary.json"), "w") as fh:
            json.dump(summary, fh, separators=(",", ":"))
        recordings.append({"slug": slug, "name": name, "label": label, "windows": n_win,
                           "faulty": sum(w["verdict"] != "HEALTHY" for w in summary_windows),
                           "entFaulty": sum(w["entVerdict"] == "FAULTY" for w in summary_windows)})

    index = {"generated": datetime.now().strftime("%Y-%m-%d %H:%M"), "script": os.path.basename(script),
             "recordings": recordings}
    with open(os.path.join(args.out, "index.json"), "w") as fh:
        json.dump(index, fh, indent=1)
    print(f"\nDone. Wrote {args.out}")


if __name__ == "__main__":
    sys.exit(main())
