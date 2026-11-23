"""
Vibration Dataset Playback Simulator  (Time Domain + FFT + VMD + Sample Entropy)
------------------------------------------------------------
Replays a saved dataset (the Excel file written by
collect_vibration_dataset.m: a "Time (s)" column + one "Group_N"
column per 5-second window) as if it were being acquired live.

The groups were cut from ONE continuous, gapless capture, so they
are joined back together in order and "streamed" at the original
sampling rate. A 5-minute dataset therefore plays for 5 minutes
(use --speed to go faster).

WINDOW 1 - live monitor (like the MATLAB real-time monitor)
  * rolling 5-second display, starts flat at zero and fills as data arrives
  * Time Domain: raw + smoothed (150 Hz low-pass) traces
  * FFT: Welch-averaged spectrum, fundamental + 2nd/3rd harmonics labelled
  * status text: VIBRATION DETECTED / IDLE

WINDOW 2 - VMD + Sample Entropy (NEW)
  Every 10 seconds of playback, the most recent 10 seconds of data are
  analysed in a background thread (so the live display never freezes):
    1. VMD splits the window into K modes (IMFs)       [default K = 7]
    2. Sample Entropy is computed for the raw signal and for every IMF
       (decimated to a fixed 500 Hz first, same as the MATLAB pipeline)
  and a pop-up figure shows: the IMF waveforms, the VMD energy distribution
  (RMS vs centre frequency), and the Sample Entropy bar chart with values.
  Results can also be appended to an Excel log (--log).

IMPORTANT: the VMD here is a numpy implementation of the standard
Dragomiretskiy & Zosso algorithm, NOT MATLAB's built-in vmd(). Its
initialisation/penalty defaults differ, so centre frequencies, RMS and
entropy values will be similar in spirit but not identical to your
MATLAB numbers. Do not mix Python and MATLAB feature values in one
Healthy-vs-Faulty comparison - analyse everything with the same tool.

Keyboard:  SPACE = pause / resume      R = restart

Requirements:  pip install numpy scipy pandas matplotlib openpyxl

Usage examples:
  python plot_dataset_time_fft.py                              # 5-minute playback
  python plot_dataset_time_fft.py --file Vibration_Dataset_Faulty_2000Hz.xlsx --label Faulty
  python plot_dataset_time_fft.py --log playback_features.xlsx  # also save results to Excel
  python plot_dataset_time_fft.py --speed 10                    # 10x faster (30 s)
  python plot_dataset_time_fft.py --vmd-k 7 --analysis-seconds 10
  python plot_dataset_time_fft.py --no-analysis                 # monitor only
  python plot_dataset_time_fft.py --static --group 5            # still plots instead
"""
import argparse
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.signal import butter, filtfilt, find_peaks, decimate
from scipy.signal.windows import hann
from scipy.spatial import cKDTree

# ---------------- USER-ADJUSTABLE DEFAULTS ----------------
DATASET_FILE   = "Vibration_Dataset_Faulty1_2000Hz.xlsx"  # <-- file to play back
SMOOTH_CUTOFF  = 150.0    # Hz, low-pass for the DISPLAYED smooth trace only
N_SEG_WELCH    = 6        # Welch segments averaged for the FFT
PEAK_RATIO_MIN = 5.0      # a fundamental must stand this many times above the median
                          # spectrum level to be labelled (avoids labelling noise as
                          # "harmonics" when the machine is idle)
MIN_SECONDS_FOR_FFT = 0.5 # need at least this much data before an FFT is drawn

ANALYSIS_SECONDS = 10.0   # run VMD + Sample Entropy every this many seconds of playback
VMD_K            = 7      # number of VMD modes (from your K-sweep result)
VMD_ALPHA        = 1000.0 # VMD bandwidth penalty
SAMPEN_M         = 2      # Sample Entropy embedding dimension
SAMPEN_FS        = 500    # Sample Entropy is always computed at this fixed rate (Hz)
# ----------------------------------------------------------


# ===================== signal processing =====================
def welch_amplitude_spectrum(x, fs, n_seg=6):
    """Hann-windowed, ~50%-overlap, averaged amplitude spectrum (units of x)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    seg_len = min(max(int(np.floor(2 * n / (n_seg + 1))), 16), n)
    step = max(seg_len // 2, 1)
    win = hann(seg_len)
    win_gain = win.mean()

    half = seg_len // 2 + 1
    accum = np.zeros(half)
    count = 0
    for s in range(0, n - seg_len + 1, step):
        seg = x[s:s + seg_len]
        seg = (seg - seg.mean()) * win
        spec = np.abs(np.fft.fft(seg) / seg_len) / win_gain
        p1 = spec[:half].copy()
        p1[1:-1] *= 2
        accum += p1
        count += 1
    freqs = fs * np.arange(half) / seg_len
    return freqs, accum / count


def estimate_fundamental_hps(freqs, p1, max_harm=3, min_freq=5.0):
    """Harmonic Product Spectrum: finds the true fundamental even when a
    harmonic is taller than the fundamental itself."""
    # Bin 0 is DC, so candidate bin j has frequency j*df and its h-th harmonic
    # sits at bin h*j. Every factor is truncated to the same length (min_len).
    min_len = len(p1) // max_harm
    hps = p1[:min_len].copy()
    for h in range(2, max_harm + 1):
        hps *= p1[::h][:min_len]                   # element j is p1[h*j]
    start = np.searchsorted(freqs[:min_len], min_freq)
    idx = start + int(np.argmax(hps[start:]))
    return freqs[idx]


def smooth_for_display(x, fs, cutoff):
    cutoff = min(cutoff, 0.45 * fs)
    b, a = butter(6, cutoff / (fs / 2), btype="low")
    return filtfilt(b, a, x)


# ===================== VMD + Sample Entropy =====================
def vmd_decompose(x, K, alpha=VMD_ALPHA, tau=0.0, tol=1e-7, max_iter=500, init="peaks"):
    """Variational Mode Decomposition (Dragomiretskiy & Zosso, 2014).

    Splits x into K band-limited modes by minimising the summed bandwidth of
    all modes subject to the modes adding back up to x (solved with ADMM:
    update each mode with a Wiener-like filter around its centre frequency,
    then move the centre frequency to the mode's spectral centre of mass).

    Returns (modes [K x N], centre_freqs [K] in cycles/sample, iterations).
    Multiply the centre frequencies by Fs to get Hz."""
    x = np.asarray(x, dtype=float)
    x = x[: len(x) // 2 * 2]                      # even length
    x = x - x.mean()
    T0 = len(x)
    half = T0 // 2
    f = np.concatenate([x[:half][::-1], x, x[-half:][::-1]])   # mirror extension
    L = len(f)
    nb = L // 2                                    # positive-frequency bins 0..L/2-1
    w = np.arange(nb) / L                          # cycles/sample
    F = np.fft.fft(f)[:nb]                         # one-sided (analytic) spectrum

    # initial centre frequencies: strongest, well-separated spectral peaks
    if init == "peaks":
        mag = np.abs(F)
        sm = np.convolve(mag, np.ones(21) / 21.0, mode="same")      # ignore noise ripple
        pk, _ = find_peaks(sm, distance=max(nb // (8 * K), 1))      # keep peaks apart
        pk = pk[np.argsort(sm[pk])[::-1]][:K]
        om = np.sort(w[pk])
        if len(om) < K:                            # not enough peaks -> pad uniformly
            extra = (0.5 / K) * np.arange(K - len(om)) + 0.5 / (2 * K)
            om = np.sort(np.concatenate([om, extra]))
    else:
        om = (0.5 / K) * np.arange(K)
    om = om.astype(float)

    u = np.zeros((K, nb), dtype=complex)
    lam = np.zeros(nb, dtype=complex)
    total = np.zeros(nb, dtype=complex)
    n_iter = 0
    for n_iter in range(1, max_iter + 1):
        diff = 0.0
        for k in range(K):
            others = total - u[k]
            new = (F - others - lam / 2) / (1.0 + alpha * (w - om[k]) ** 2)
            diff += np.sum(np.abs(new - u[k]) ** 2) / L
            total += new - u[k]
            u[k] = new
            p = np.abs(u[k]) ** 2
            s = p.sum()
            if s > 0:
                om[k] = np.dot(w, p) / s           # centre of mass of the mode's power
        if tau:
            lam = lam + tau * (total - F)
        if diff < tol:
            break

    order = np.argsort(om)
    om, u = om[order], u[order]
    modes = np.zeros((K, T0))
    for k in range(K):
        full = np.zeros(L, dtype=complex)
        full[:nb] = u[k]
        full[L - nb + 1:] = np.conj(u[k][1:])[::-1]     # conjugate-symmetric negative half
        modes[k] = np.real(np.fft.ifft(full))[half:half + T0]
    return modes, om, n_iter


def sample_entropy(x, m=2, r=None):
    """Sample Entropy (Richman & Moorman, 2000): -ln(A/B), where B counts pairs
    of length-m templates within tolerance r (Chebyshev distance, no
    self-matches) and A counts the same pairs at length m+1. Uses a KD-tree
    for the neighbour counting; verified identical to a brute-force version."""
    x = np.asarray(x, dtype=float)
    if r is None:
        r = 0.2 * np.std(x, ddof=1)
    n_vec = len(x) - m
    if n_vec < 2:
        return np.nan

    def count_pairs(dim):
        X = np.column_stack([x[k:k + n_vec] for k in range(dim)])
        tree = cKDTree(X)
        c = tree.count_neighbors(tree, r, p=np.inf)    # ordered pairs incl. self
        return (c - n_vec) / 2.0                       # unordered pairs, i < j

    B, A = count_pairs(m), count_pairs(m + 1)
    return np.nan if (A == 0 or B == 0) else float(-np.log(A / B))


def decimate_if_needed(x, fs, target_fs):
    if fs > target_fs:
        q = int(round(fs / target_fs))
        return decimate(x, q), fs / q
    return np.asarray(x, dtype=float), float(fs)


def analyze_window(x, fs, K, m, sampen_fs, t_start):
    """VMD + Sample Entropy of ONE captured window (same pipeline as the MATLAB
    combined script: VMD on the full-rate data, then Sample Entropy of the raw
    signal and of every IMF after decimating to a fixed rate)."""
    t0 = time.perf_counter()
    x = np.asarray(x, dtype=float)
    modes, om, n_iter = vmd_decompose(x, K)
    x = x[:modes.shape[1]] - x[:modes.shape[1]].mean()
    centers = om * fs
    mode_rms = np.sqrt(np.mean(modes ** 2, axis=1))

    xd, fs_dec = decimate_if_needed(x, fs, sampen_fs)
    sampen_raw = sample_entropy(xd, m, 0.2 * np.std(xd, ddof=1))
    sampen = []
    for k in range(K):
        md, _ = decimate_if_needed(modes[k], fs, sampen_fs)
        sampen.append(sample_entropy(md, m, 0.2 * np.std(md, ddof=1)))

    return dict(t_start=t_start, t_end=t_start + len(x) / fs, fs=fs, fs_dec=fs_dec, K=K, m=m,
                x=x, modes=modes, centers=centers, mode_rms=mode_rms,
                raw_rms=float(np.sqrt(np.mean(x ** 2))), sampen_raw=sampen_raw,
                sampen=np.array(sampen), n_iter=n_iter, elapsed=time.perf_counter() - t0)


def result_to_row(res, label):
    row = {"Timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "Label": label,
           "WindowStart_s": round(res["t_start"], 3), "WindowEnd_s": round(res["t_end"], 3),
           "AcqFs": res["fs"], "AnalysisFs": res["fs_dec"], "K": res["K"],
           "RawRMS": res["raw_rms"], "SampEn_Raw": res["sampen_raw"]}
    for k in range(res["K"]):
        row[f"IMF{k + 1}_CenterFreqHz"] = res["centers"][k]
        row[f"IMF{k + 1}_RMS"] = res["mode_rms"][k]
        row[f"IMF{k + 1}_SampEn"] = res["sampen"][k]
    return row


def load_dataset(path):
    df = pd.read_excel(path)
    time_col = df.iloc[:, 0].to_numpy(dtype=float)
    fs = int(round(1.0 / np.mean(np.diff(time_col))))
    group_names = list(df.columns[1:])
    data = df.iloc[:, 1:].to_numpy(dtype=float)   # samples x groups
    return time_col, fs, group_names, data


def fmt_clock(seconds):
    seconds = int(round(seconds))
    return f"{seconds // 60}:{seconds % 60:02d}"


# ===================== playback simulator =====================
class PlaybackSimulator:
    """Replays a continuous stream as a rolling window, like the live monitor,
    and runs VMD + Sample Entropy on every ANALYSIS_SECONDS of playback."""

    def __init__(self, data, fs, xmax=None, speed=1.0, analysis=True,
                 analysis_seconds=ANALYSIS_SECONDS, vmd_k=VMD_K, sampen_fs=SAMPEN_FS,
                 sampen_m=SAMPEN_M, label="Unlabelled", log_file=None):
        self.fs = fs
        self.group_len = data.shape[0]                 # samples per 5 s group
        self.num_groups = data.shape[1]
        # Groups are consecutive windows of one gapless capture -> join in order.
        stream = data.T.reshape(-1)
        # Remove the sensor's constant zero-g offset once (like the Teensy's
        # calibrated bias), so "no vibration" really sits at zero.
        self.stream = stream - np.median(stream)
        self.total = len(self.stream)
        self.W = self.group_len                        # rolling display = one group
        self.speed = speed
        self.xmax = xmax if xmax else fs / 2

        cutoff = min(SMOOTH_CUTOFF, 0.45 * fs)
        self.b, self.a = butter(6, cutoff / (fs / 2), btype="low")

        self.y_lim = 0.02            # smoothed auto-scale state (time axis)
        self.fft_top = 1e-3          # smoothed auto-scale state (FFT axis)
        self.labels = []
        self.ani = None
        self._tick = None

        # ---- VMD + Sample Entropy settings / state ----
        self.win = int(round(analysis_seconds * fs))
        self.analysis = bool(analysis) and self.win <= self.total
        if analysis and not self.analysis:
            print("Dataset is shorter than the analysis window - VMD/Sample Entropy disabled.")
        self.vmd_k, self.sampen_fs, self.sampen_m = vmd_k, sampen_fs, sampen_m
        self.label, self.log_file = label, log_file
        self.executor = ThreadPoolExecutor(max_workers=1) if self.analysis else None
        self.pending = None
        self.last_idx = 0
        self.rows = []
        self.n_results = 0
        self.fig2 = None

        t_axis = np.arange(1, self.W + 1) / fs
        self.fig, (self.ax_t, self.ax_f) = plt.subplots(1, 2, figsize=(14, 6.4))

        self.raw, = self.ax_t.plot(t_axis, np.zeros(self.W), color="lightsteelblue", lw=0.6, label="Raw")
        self.smooth, = self.ax_t.plot(t_axis, np.zeros(self.W), color="blue", lw=1.1,
                                      label=f"Smoothed ({cutoff:.0f} Hz LP)")
        self.ax_t.set_xlim(t_axis[0], t_axis[-1])
        self.ax_t.set_ylim(-self.y_lim, self.y_lim)
        self.ax_t.set_xlabel("Time (s)")
        self.ax_t.set_ylabel("Amplitude (g)")
        self.ax_t.set_title(f"Time Domain Waveform (Fs = {fs} Hz)")
        self.ax_t.grid(True, alpha=0.4)
        self.ax_t.legend(loc="upper right")
        self.status = self.ax_t.text(0.02, 0.95, "Starting...", transform=self.ax_t.transAxes,
                                     fontweight="bold", color="gray", va="top")

        self.fft_line, = self.ax_f.plot([0, fs / 2], [0, 0], color="red", lw=1.2)
        self.peaks, = self.ax_f.plot([], [], "kv", markerfacecolor="yellow", markersize=8)
        self.ax_f.set_xlim(0, self.xmax)
        self.ax_f.set_ylim(0, self.fft_top)
        self.ax_f.set_xlabel("Frequency (Hz)")
        self.ax_f.set_ylabel("Amplitude (g)")
        self.ax_f.set_title("FFT Spectrum (3 Harmonic Modes)")
        self.ax_f.grid(True, alpha=0.4)

        info = (f"VMD + Sample Entropy: first analysis after {analysis_seconds:g} s of data..."
                if self.analysis else "VMD + Sample Entropy: disabled")
        self.info = self.fig.text(0.5, 0.012, info, ha="center", fontsize=10, color="navy")
        self.fig.tight_layout(rect=(0, 0.04, 1, 0.94))
        self.render(0)

        if self.analysis:
            self.fig2 = plt.figure(figsize=(15, 8.5))
            self.fig2.text(0.5, 0.5, f"Waiting for the first {analysis_seconds:g} s of data...\n"
                           "(VMD + Sample Entropy results will appear here)",
                           ha="center", va="center", fontsize=14, color="gray")
            for fg in (self.fig, self.fig2):
                fg.canvas.mpl_connect("close_event", self._on_close)

    def _on_close(self, _event):
        if self.executor is not None:
            self.executor.shutdown(wait=False, cancel_futures=True)

    # ---- draw the state of the system after `n` samples have "arrived" ----
    def render(self, n):
        n = int(min(max(n, 0), self.total))
        lo = max(0, n - self.W)
        valid = self.stream[lo:n]
        nv = len(valid)

        buf = np.zeros(self.W)                 # no input -> flat zero
        if nv:
            buf[self.W - nv:] = valid          # newest sample at the right edge

        # --- time domain ---
        smooth = filtfilt(self.b, self.a, buf) if nv > 30 else buf
        self.raw.set_ydata(buf)
        self.smooth.set_ydata(smooth)
        target = max(1.15 * np.max(np.abs(buf)), 0.02)
        self.y_lim = 0.8 * self.y_lim + 0.2 * target
        self.ax_t.set_ylim(-self.y_lim, self.y_lim)

        # --- frequency domain ---
        if nv >= self.fs * MIN_SECONDS_FOR_FFT:
            freqs, p1 = welch_amplitude_spectrum(valid, self.fs, N_SEG_WELCH)
        else:
            freqs, p1 = np.array([0.0, self.fs / 2]), np.array([0.0, 0.0])
        self.fft_line.set_data(freqs, p1)
        top = max(1.25 * float(np.max(p1)), 1e-4)
        self.fft_top = 0.8 * self.fft_top + 0.2 * top
        self.ax_f.set_ylim(0, self.fft_top)

        # --- fundamental / harmonics, only if a clear fundamental exists ---
        for t in self.labels:
            t.remove()
        self.labels = []
        detected = False
        if nv >= self.fs and len(p1) > 12:
            f0 = estimate_fundamental_hps(freqs, p1)
            amp0 = p1[np.argmin(np.abs(freqs - f0))]
            detected = amp0 > PEAK_RATIO_MIN * np.median(p1)
        if detected:
            xs, ys = [], []
            for k in (1, 2, 3):
                fk = f0 * k
                if fk <= freqs[-1]:
                    ak = p1[np.argmin(np.abs(freqs - fk))]
                    xs.append(fk)
                    ys.append(ak)
                    self.labels.append(self.ax_f.annotate(
                        f"Mode {k}: {fk:.1f} Hz", (fk, ak), textcoords="offset points",
                        xytext=(8, 8 + 16 * (3 - k)), fontsize=9,
                        arrowprops=dict(arrowstyle="-", color="gray", lw=0.6)))
            self.peaks.set_data(xs, ys)
            self.status.set_text("VIBRATION DETECTED")
            self.status.set_color((0, 0.55, 0))
        else:
            self.peaks.set_data([], [])
            if nv < self.fs:
                self.status.set_text("Acquiring...")
                self.status.set_color("gray")
            else:
                self.status.set_text("IDLE - no significant vibration")
                self.status.set_color((0.6, 0.6, 0.6))

        # --- simulation clock ---
        group_now = min(self.num_groups, max(1, (n - 1) // self.group_len + 1)) if n else 0
        self.fig.suptitle(
            f"PLAYBACK SIMULATION   t = {fmt_clock(n / self.fs)} / {fmt_clock(self.total / self.fs)}"
            f"   |   Group {group_now}/{self.num_groups}   |   speed x{self.speed:g}"
            f"   |   [SPACE] pause  [R] restart",
            fontsize=11, fontweight="bold")

    # ---- VMD + Sample Entropy: background job, drawn from the main thread ----
    def _maybe_start_analysis(self, n, start_n):
        delivered = n - start_n
        idx = delivered // self.win
        if idx > self.last_idx and self.pending is None and n >= self.win:
            self.last_idx = idx
            window = self.stream[n - self.win:n].copy()
            t_start = (n - self.win) / self.fs
            self.pending = self.executor.submit(analyze_window, window, self.fs, self.vmd_k,
                                                self.sampen_m, self.sampen_fs, t_start)
            self.info.set_text(f"VMD + Sample Entropy: analysing {t_start:.0f}-{n / self.fs:.0f} s ...")

    def _poll_analysis(self):
        if self.pending is not None and self.pending.done():
            fut, self.pending = self.pending, None
            try:
                res = fut.result()
            except Exception as exc:                       # keep the monitor alive
                print(f"VMD / Sample Entropy analysis failed: {exc}")
                self.info.set_text(f"VMD + Sample Entropy: analysis failed ({exc})")
                return
            self.show_analysis(res)

    def show_analysis(self, res):
        self.n_results += 1
        K, fs = res["K"], res["fs"]
        print(f"[analysis #{self.n_results}] window {res['t_start']:.0f}-{res['t_end']:.0f} s | "
              f"Raw SampEn = {res['sampen_raw']:.3f} | IMF SampEn = "
              f"{np.round(res['sampen'], 3).tolist()} | {res['elapsed']:.1f} s")
        self.info.set_text(
            f"VMD + Sample Entropy #{self.n_results}: window {res['t_start']:.0f}-{res['t_end']:.0f} s"
            f"   |   Raw Sample Entropy = {res['sampen_raw']:.3f}   |   K = {K}")

        # ---- optional Excel log ----
        if self.log_file:
            self.rows.append(result_to_row(res, self.label))
            try:
                pd.DataFrame(self.rows).to_excel(self.log_file, index=False)
            except Exception as exc:
                print(f"Could not write {self.log_file}: {exc}")

        # ---- the pop-up figure ----
        fig = self.fig2
        fig.clf()
        gs = fig.add_gridspec(K + 1, 3, width_ratios=[1.35, 1, 1], hspace=0.95, wspace=0.28,
                              left=0.05, right=0.985, top=0.9, bottom=0.07)
        t = np.arange(len(res["x"])) / fs
        for k in range(K):
            ax = fig.add_subplot(gs[k, 0])
            ax.plot(t, res["modes"][k], color="blue", lw=0.6)
            ax.set_title(f"IMF {k + 1}  ({res['centers'][k]:.1f} Hz, RMS {res['mode_rms'][k]:.4f} g)",
                         fontsize=8, pad=2)
            ax.set_xlim(t[0], t[-1])
            ax.tick_params(labelsize=7)
            ax.grid(True, alpha=0.3)
            if k < K - 1:
                ax.set_xticklabels([])
        ax = fig.add_subplot(gs[K, 0])
        ax.plot(t, res["x"], color="black", lw=0.5)
        ax.set_title("Original signal (DC removed)", fontsize=8, pad=2)
        ax.set_xlim(t[0], t[-1])
        ax.set_xlabel("Time within window (s)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)

        ax_e = fig.add_subplot(gs[:, 1])
        ax_e.bar(res["centers"], res["mode_rms"], width=0.02 * fs / 2, color="tab:blue",
                 edgecolor="black", linewidth=0.5)
        ax_e.set_xlabel("Mode centre frequency (Hz)")
        ax_e.set_ylabel("Mode RMS energy (g)")
        ax_e.set_title("VMD Energy Distribution")
        ax_e.grid(True, alpha=0.4)

        ax_s = fig.add_subplot(gs[:, 2])
        vals = np.concatenate([[res["sampen_raw"]], res["sampen"]])
        names = ["Raw"] + [f"IMF{k + 1}" for k in range(K)]
        colors = ["dimgray"] + ["tab:blue"] * K
        bars = ax_s.bar(range(K + 1), np.nan_to_num(vals), color=colors, edgecolor="black", linewidth=0.5)
        for bar, v in zip(bars, vals):
            txt = "n/a" if np.isnan(v) else f"{v:.3f}"
            ax_s.annotate(txt, (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                          ha="center", va="bottom", fontsize=8, xytext=(0, 2), textcoords="offset points")
        ax_s.set_xticks(range(K + 1))
        ax_s.set_xticklabels(names, fontsize=8)
        ax_s.set_ylabel("Sample Entropy")
        ax_s.set_ylim(0, max(1.0, float(np.nanmax(vals)) * 1.15))
        ax_s.set_title(f"Sample Entropy Feature Vector (at {res['fs_dec']:g} Hz)")
        ax_s.grid(True, alpha=0.4, axis="y")

        fig.suptitle(
            f"VMD (K={K}) + Sample Entropy  -  analysis #{self.n_results}  -  "
            f"window {res['t_start']:.0f}-{res['t_end']:.0f} s of playback  -  {self.label}   "
            f"[computed in {res['elapsed']:.1f} s]", fontsize=12, fontweight="bold")
        fig.canvas.draw_idle()

    # ---- real-time driver (wall-clock based, so 5 min really is 5 min) ----
    def run(self, start_group=1, interval_ms=100):
        state = dict(t0=time.perf_counter(), paused=False, pause_t=0.0,
                     paused_total=0.0, done=False,
                     start_n=(max(1, start_group) - 1) * self.group_len)

        def tick(_frame):
            if self.analysis:
                self._poll_analysis()
            if state["done"]:
                if self.analysis and self.pending is None:
                    self._maybe_start_analysis(self.total, state["start_n"])   # final window
                if self.pending is None:                  # all results have been shown
                    self.ani.event_source.stop()
                return
            if state["paused"]:
                return
            elapsed = (time.perf_counter() - state["t0"] - state["paused_total"]) * self.speed
            n = min(self.total, state["start_n"] + int(elapsed * self.fs))
            self.render(n)
            if self.analysis:
                self._maybe_start_analysis(n, state["start_n"])
            if n >= self.total:
                state["done"] = True
                self.status.set_text("SIMULATION COMPLETE")
                self.status.set_color("navy")

        def on_key(event):
            if event.key == " ":
                now = time.perf_counter()
                if state["paused"]:
                    state["paused_total"] += now - state["pause_t"]
                    state["paused"] = False
                else:
                    state["pause_t"] = now
                    state["paused"] = True
            elif event.key in ("r", "R"):
                state.update(t0=time.perf_counter(), paused=False, paused_total=0.0, done=False)
                self.last_idx = 0
                self.n_results = 0
                self.rows = []
                self.ani.event_source.start()

        self._tick = tick
        self.fig.canvas.mpl_connect("key_press_event", on_key)
        self.ani = FuncAnimation(self.fig, tick, interval=interval_ms,
                                 cache_frame_data=False, blit=False)
        plt.show()


# ===================== static (still) plots, old behaviour =====================
def plot_single_group(time_col, fs, name, x, save=False, xmax=None):
    x = x - x.mean()
    smooth = smooth_for_display(x, fs, SMOOTH_CUTOFF)
    freqs, p1 = welch_amplitude_spectrum(x, fs, N_SEG_WELCH)

    fig, (ax_t, ax_f) = plt.subplots(1, 2, figsize=(14, 5.5))
    ax_t.plot(time_col, x, color="lightsteelblue", lw=0.6, label="Raw")
    ax_t.plot(time_col, smooth, color="blue", lw=1.0, label=f"Smoothed ({SMOOTH_CUTOFF:.0f} Hz LP)")
    ax_t.set_xlabel("Time (s)")
    ax_t.set_ylabel("Amplitude (g)")
    ax_t.set_title(f"Time Domain Waveform - {name} (Fs = {fs} Hz)")
    ax_t.set_xlim(time_col[0], time_col[-1])
    ax_t.grid(True, alpha=0.4)
    ax_t.legend(loc="upper right")

    ax_f.plot(freqs, p1, color="red", lw=1.2)
    ax_f.set_xlabel("Frequency (Hz)")
    ax_f.set_ylabel("Amplitude (g)")
    ax_f.set_title("FFT Spectrum (3 Harmonic Modes)")
    ax_f.set_xlim(0, xmax if xmax else fs / 2)
    ax_f.set_ylim(bottom=0)
    ax_f.grid(True, alpha=0.4)

    f0 = estimate_fundamental_hps(freqs, p1)
    if p1[np.argmin(np.abs(freqs - f0))] > PEAK_RATIO_MIN * np.median(p1):
        for k in (1, 2, 3):
            fk = f0 * k
            if fk <= freqs[-1]:
                ak = p1[np.argmin(np.abs(freqs - fk))]
                ax_f.plot(fk, ak, "kv", markerfacecolor="yellow", markersize=8)
                ax_f.annotate(f"Mode {k}: {fk:.1f} Hz", (fk, ak), textcoords="offset points",
                              xytext=(8, 8 + 16 * (3 - k)), fontsize=9,
                              arrowprops=dict(arrowstyle="-", color="gray", lw=0.6))
    else:
        ax_f.text(0.02, 0.95, "No clear fundamental (idle / noise only)",
                  transform=ax_f.transAxes, color="gray", fontweight="bold", va="top")
    fig.tight_layout()
    if save:
        fig.savefig(f"time_fft_{name}.png", dpi=150)
        print(f"Saved time_fft_{name}.png")
    return fig


def plot_average_spectrum(fs, data, save=False):
    spectra = []
    for g in range(data.shape[1]):
        freqs, p1 = welch_amplitude_spectrum(data[:, g], fs, N_SEG_WELCH)
        spectra.append(p1)
    spectra = np.array(spectra)
    mean_p1, std_p1 = spectra.mean(axis=0), spectra.std(axis=0)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.plot(freqs, mean_p1, color="darkred", lw=1.3, label=f"Mean of {data.shape[1]} groups")
    ax.fill_between(freqs, np.maximum(mean_p1 - std_p1, 0), mean_p1 + std_p1,
                    color="red", alpha=0.2, label="+/- 1 std across groups")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Amplitude (g)")
    ax.set_title(f"Average FFT Spectrum across all groups (Fs = {fs} Hz)")
    ax.set_xlim(0, fs / 2)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.4)
    ax.legend()
    fig.tight_layout()
    if save:
        fig.savefig("average_fft_all_groups.png", dpi=150)
        print("Saved average_fft_all_groups.png")
    return fig


# ===================== main =====================
def main():
    ap = argparse.ArgumentParser(description="Playback simulation of a saved vibration dataset.")
    ap.add_argument("--file", default=DATASET_FILE, help="Excel dataset to load")
    ap.add_argument("--speed", type=float, default=1.0, help="playback speed (1 = real time, 5 min)")
    ap.add_argument("--start-group", type=int, default=1, help="group to start playing from (1-based)")
    ap.add_argument("--xmax", type=float, default=None, help="zoom the FFT x-axis to this max frequency (Hz)")
    ap.add_argument("--interval", type=int, default=100, help="screen refresh interval in ms")
    ap.add_argument("--no-analysis", action="store_true", help="disable the VMD + Sample Entropy analysis")
    ap.add_argument("--analysis-seconds", type=float, default=ANALYSIS_SECONDS,
                    help="run VMD + Sample Entropy every this many seconds of playback")
    ap.add_argument("--vmd-k", type=int, default=VMD_K, help="number of VMD modes")
    ap.add_argument("--label", default="Unlabelled", help="condition label stored in the log (e.g. Healthy)")
    ap.add_argument("--log", default=None, help="Excel file to append VMD/Sample Entropy results to")
    ap.add_argument("--static", action="store_true", help="draw still plots instead of the simulation")
    ap.add_argument("--group", type=int, default=1, help="(static mode) group number to plot")
    ap.add_argument("--save", action="store_true", help="(static mode) save PNGs instead of opening windows")
    args = ap.parse_args()

    print(f"Loading {args.file} ... (large Excel files can take a little while)")
    time_col, fs, group_names, data = load_dataset(args.file)
    print(f"Loaded {data.shape[1]} groups x {data.shape[0]} samples, inferred Fs = {fs} Hz "
          f"-> {data.size / fs / 60:.1f} minutes of data")

    if args.static:
        if not 1 <= args.group <= data.shape[1]:
            raise SystemExit(f"--group must be between 1 and {data.shape[1]}")
        plot_single_group(time_col, fs, group_names[args.group - 1], data[:, args.group - 1],
                          args.save, args.xmax)
        plot_average_spectrum(fs, data, args.save)
        if not args.save:
            plt.show()
        return

    if not 1 <= args.start_group <= data.shape[1]:
        raise SystemExit(f"--start-group must be between 1 and {data.shape[1]}")
    sim = PlaybackSimulator(data, fs, xmax=args.xmax, speed=args.speed,
                            analysis=not args.no_analysis, analysis_seconds=args.analysis_seconds,
                            vmd_k=args.vmd_k, label=args.label, log_file=args.log)
    remaining = (data.shape[1] - args.start_group + 1) * data.shape[0] / fs / args.speed
    print(f"Starting playback - about {fmt_clock(remaining)} (min:sec) at x{args.speed:g} speed. "
          f"SPACE = pause/resume, R = restart.")
    if sim.analysis:
        print(f"VMD (K={args.vmd_k}) + Sample Entropy will run every {args.analysis_seconds:g} s of "
              f"playback and appear in a second window.")
    sim.run(start_group=args.start_group, interval_ms=args.interval)


if __name__ == "__main__":
    main()