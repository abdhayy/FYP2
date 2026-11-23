"""
Vibration Dataset Playback Simulator  -  Time Domain + FFT + ALL-ENTROPY analysis
------------------------------------------------------------
Replays a saved dataset (the Excel file written by
collect_vibration_dataset.m: a "Time (s)" column + one "Group_N"
column per 5-second window) as if it were being acquired live.
The groups are joined back into the original continuous capture and
"streamed" at the original sampling rate: a 5-minute dataset plays for
5 minutes (use --speed to go faster).

WINDOW 1 - live monitor  (identical to plot_dataset_time_fft.py)
  * rolling 5-second display, flat at zero until data "arrives"
  * Time Domain: raw + smoothed (150 Hz low-pass) traces
  * FFT: Welch-averaged spectrum, fundamental + 2nd/3rd harmonics labelled
  * status text: VIBRATION DETECTED / IDLE

WINDOW 2 - entropy analysis (NEW), refreshed every 10 s of playback
  The latest 10 s of data are analysed in a background thread (the live
  display never freezes) with SEVEN entropy measures:

    Sample Entropy        SampEn   regularity of the waveform (0 = perfectly regular)
    Approximate Entropy   ApEn     older relative of SampEn (includes self-matches)
    Fuzzy Entropy         FuzzyEn  SampEn with a smooth similarity, more stable on short data
    Permutation Entropy   PermEn   randomness of the up/down ordering of samples (0..1)
    Dispersion Entropy    DispEn   randomness of amplitude-class patterns (0..1)
    Spectral Entropy      SpecEn   how spread-out the power spectrum is (0..1)
    Shannon Entropy       ShanEn   flatness of the amplitude histogram (0..1) - NOTE: this
                                   ignores time order, it is not a regularity measure

  plus Multiscale Sample Entropy (scales 1-10) of the raw signal.
  Each measure is computed for the raw signal AND (unless --no-vmd) for every
  VMD mode (IMF), shown as a heat-map table with the actual values printed in
  every cell, next to the IMF waveforms, the MSE curve and the VMD band-RMS bars.
  Results can be appended to an Excel log (--log).

ANALYSIS BAND WARNING
  Entropy is computed after decimating to --analysis-fs (default 500 Hz), whose
  anti-alias filter removes everything above ~0.4 x 500 = 200 Hz. Therefore:
    * the "Raw" entropies describe only the 0-200 Hz part of the vibration, and
    * an IMF whose centre frequency is above ~200 Hz (e.g. a 550 Hz resonance
      band) is almost entirely filtered away, so its entropies are meaningless.
  Such IMFs are hatched grey and marked with * in the figure, and flagged
  (IMFk_InBand = False) in the Excel log. Use --analysis-fs 2000 to analyse the
  full band (slower).

VMD here is a numpy implementation of the standard Dragomiretskiy & Zosso
algorithm, NOT MATLAB's vmd(), so numbers are similar in spirit but not
identical to your MATLAB results - do not mix the two in one comparison.

Keyboard:  SPACE = pause / resume      R = restart

Requirements:  pip install numpy scipy pandas matplotlib openpyxl

Usage examples:
  python plot_dataset_time_fft_entropy.py --file Vibration_Dataset_Healthy_2000Hz.xlsx --label Healthy
  python plot_dataset_time_fft_entropy.py --log entropy_features.xlsx    # save results to Excel
  python plot_dataset_time_fft_entropy.py --no-vmd                       # entropies of the raw signal only
  python plot_dataset_time_fft_entropy.py --speed 10 --xmax 300
  python plot_dataset_time_fft_entropy.py --analysis-fs 2000             # full band (slower)
  python plot_dataset_time_fft_entropy.py --static --group 5            # still plots instead
"""
import argparse
import math
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
from scipy.signal import butter, filtfilt, find_peaks, decimate
from scipy.signal.windows import hann
from scipy.spatial import cKDTree
from scipy.stats import norm

# ---------------- USER-ADJUSTABLE DEFAULTS ----------------
DATASET_FILE   = "Vibration_Dataset_Healthy_2000Hz.xlsx"  # <-- file to play back
SMOOTH_CUTOFF  = 150.0    # Hz, low-pass for the DISPLAYED smooth trace only
N_SEG_WELCH    = 6        # Welch segments averaged for the FFT
PEAK_RATIO_MIN = 5.0      # a fundamental must stand this many times above the median
                          # spectrum level to be labelled (avoids labelling noise as
                          # "harmonics" when the machine is idle)
MIN_SECONDS_FOR_FFT = 0.5 # need at least this much data before an FFT is drawn

ANALYSIS_SECONDS = 10.0   # run the entropy analysis every this many seconds of playback
USE_VMD          = True   # also compute every entropy for each VMD mode
VMD_K            = 7      # number of VMD modes (from your K-sweep result)
VMD_ALPHA        = 1000.0 # VMD bandwidth penalty
ANALYSIS_FS      = 500    # entropies are computed at this fixed rate (Hz)

# entropy parameters (common literature defaults)
EN_M          = 2         # embedding dimension for SampEn / ApEn / FuzzyEn
EN_R          = 0.2       # tolerance = EN_R x standard deviation
FUZZY_N       = 2         # fuzzy membership exponent
PERM_ORDER    = 5         # permutation entropy order (5! = 120 patterns)
PERM_DELAY    = 1
DISP_M        = 2         # dispersion entropy embedding dimension
DISP_C        = 6         # dispersion entropy number of classes
SHANNON_BINS  = 32
MSE_MAX_SCALE = 10        # multiscale entropy scales 1..10
FUZZY_MAX_ROWS = 6000     # FuzzyEn compares every template with all others; for very long
                          # series (e.g. --analysis-fs 2000) only this many templates are used
# ----------------------------------------------------------


# ===================== time-domain / FFT processing (unchanged) =====================
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


# ===================== VMD =====================
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


# ===================== the entropy measures =====================
def _std(x):
    return np.std(x, ddof=1)


def sample_entropy(x, m=2, r=None):
    """Sample Entropy (Richman & Moorman, 2000): -ln(A/B). B = number of pairs of
    length-m templates closer than r (Chebyshev distance, no self-matches);
    A = the same pairs still matching at length m+1. Lower = more regular."""
    x = np.asarray(x, dtype=float)
    if r is None:
        r = EN_R * _std(x)
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


def approximate_entropy(x, m=2, r=None):
    """Approximate Entropy (Pincus, 1991): phi(m) - phi(m+1), where phi(k) is the mean
    log of the fraction of length-k templates within r of each template
    (INCLUDING the self-match, which biases ApEn low on short data)."""
    x = np.asarray(x, dtype=float)
    if r is None:
        r = EN_R * _std(x)
    N = len(x)
    if N - m < 2:
        return np.nan

    def phi(dim):
        n = N - dim + 1
        X = np.column_stack([x[k:k + n] for k in range(dim)])
        counts = cKDTree(X).query_ball_point(X, r, p=np.inf, return_length=True)
        return float(np.mean(np.log(counts / n)))

    return phi(m) - phi(m + 1)


def fuzzy_entropy(x, m=2, r=EN_R, n_exp=FUZZY_N, chunk=256, max_rows=FUZZY_MAX_ROWS):
    """Fuzzy Entropy (Chen et al., 2007): like SampEn, but two templates are
    'similar' to a degree exp(-(d^n)/r) instead of a hard yes/no, after
    standardising the series and removing each template's own mean."""
    x = np.asarray(x, dtype=float)
    z = (x - x.mean()) / _std(x)                  # standardised, so r is relative to the std
    n_vec = len(z) - m
    if n_vec < 2:
        return np.nan
    stride = max(1, int(math.ceil(n_vec / max_rows)))
    rows = np.arange(0, n_vec, stride)

    def phi(dim):
        X = np.column_stack([z[k:k + n_vec] for k in range(dim)])
        X = (X - X.mean(axis=1, keepdims=True)).astype(np.float32)   # remove template baseline
        total = 0.0
        for s in range(0, len(rows), chunk):
            sel = rows[s:s + chunk]
            blk = X[sel]
            d = np.abs(blk[:, 0:1] - X[:, 0][None, :])
            for k in range(1, dim):
                d = np.maximum(d, np.abs(blk[:, k:k + 1] - X[:, k][None, :]))   # Chebyshev
            sim = np.exp(-(d ** n_exp) / r)
            total += float(sim.sum()) - len(sel)                    # drop self-matches (=1 each)
        return total / (len(rows) * (n_vec - 1))

    pm, pm1 = phi(m), phi(m + 1)
    return np.nan if (pm <= 0 or pm1 <= 0) else float(np.log(pm) - np.log(pm1))


def permutation_entropy(x, order=PERM_ORDER, delay=PERM_DELAY, normalize=True):
    """Permutation Entropy (Bandt & Pompe, 2002): Shannon entropy of the ordinal
    (up/down ranking) patterns of 'order' consecutive samples. 0 = always the
    same pattern, 1 = all patterns equally likely (white noise)."""
    x = np.asarray(x, dtype=float)
    n = len(x) - (order - 1) * delay
    if n < 2:
        return np.nan
    idx = np.arange(order) * delay
    X = x[np.arange(n)[:, None] + idx[None, :]]
    patterns = np.argsort(X, axis=1, kind="stable")
    code = (patterns * (order ** np.arange(order))).sum(axis=1)
    p = np.unique(code, return_counts=True)[1] / n
    H = -np.sum(p * np.log(p))
    return float(H / np.log(math.factorial(order))) if normalize else float(H)


def dispersion_entropy(x, m=DISP_M, c=DISP_C, delay=1, normalize=True):
    """Dispersion Entropy (Rostaghi & Azami, 2016): map the signal to c amplitude
    classes via the normal CDF, then take the Shannon entropy of the class
    patterns of m consecutive samples. Normalised to 0..1."""
    x = np.asarray(x, dtype=float)
    y = norm.cdf(x, loc=x.mean(), scale=_std(x))
    z = np.clip(np.round(c * y + 0.5).astype(int), 1, c)
    n = len(z) - (m - 1) * delay
    if n < 2:
        return np.nan
    idx = np.arange(m) * delay
    Z = z[np.arange(n)[:, None] + idx[None, :]]
    code = (Z - 1) @ (c ** np.arange(m))
    p = np.unique(code, return_counts=True)[1] / n
    H = -np.sum(p * np.log(p))
    return float(H / np.log(c ** m)) if normalize else float(H)


def spectral_entropy(x, fs, n_seg=N_SEG_WELCH):
    """Shannon entropy of the normalised power spectrum, divided by ln(#bins).
    Few strong peaks -> low; flat (white-noise-like) spectrum -> close to 1."""
    _, amp = welch_amplitude_spectrum(x, fs, n_seg)
    p = amp ** 2
    p = p / p.sum()
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)) / np.log(len(amp)))


def shannon_entropy(x, bins=SHANNON_BINS):
    """Shannon (information) entropy of the amplitude histogram, normalised to 0..1.
    Measures how evenly the amplitudes are spread - it does NOT see time order."""
    counts, _ = np.histogram(np.asarray(x, dtype=float), bins=bins)
    p = counts[counts > 0] / counts.sum()
    return float(-np.sum(p * np.log(p)) / np.log(bins))


def multiscale_sample_entropy(x, max_scale=MSE_MAX_SCALE, m=EN_M):
    """Multiscale Sample Entropy (Costa et al., 2002): SampEn of the signal
    coarse-grained (averaged in non-overlapping blocks) at scales 1..max_scale,
    with the tolerance fixed from the original series."""
    x = np.asarray(x, dtype=float)
    r = EN_R * _std(x)
    out = []
    for tau in range(1, max_scale + 1):
        n = len(x) // tau
        y = x[:n * tau].reshape(n, tau).mean(axis=1)
        out.append(sample_entropy(y, m, r))
    return np.array(out)


# key, row label for the table, function(series, fs)
ENTROPY_DEFS = [
    ("SampEn",  "Sample En",           lambda s, fs: sample_entropy(s, EN_M, EN_R * _std(s))),
    ("ApEn",    "Approximate En",      lambda s, fs: approximate_entropy(s, EN_M, EN_R * _std(s))),
    ("FuzzyEn", "Fuzzy En",            lambda s, fs: fuzzy_entropy(s)),
    ("PermEn",  "Permutation En (0-1)", lambda s, fs: permutation_entropy(s)),
    ("DispEn",  "Dispersion En (0-1)",  lambda s, fs: dispersion_entropy(s)),
    ("SpecEn",  "Spectral En (0-1)",    lambda s, fs: spectral_entropy(s, fs)),
    ("ShanEn",  "Shannon En (0-1)",     lambda s, fs: shannon_entropy(s)),
]


def compute_all_entropies(series, fs):
    """All seven entropy measures for one series -> dict {key: value}."""
    series = np.asarray(series, dtype=float)
    if len(series) < 50 or not np.isfinite(series).all() or _std(series) < 1e-15:
        return {key: np.nan for key, _, _ in ENTROPY_DEFS}
    out = {}
    with np.errstate(all="ignore"):
        for key, _, fn in ENTROPY_DEFS:
            try:
                out[key] = float(fn(series, fs))
            except Exception:
                out[key] = np.nan
    return out


def decimate_if_needed(x, fs, target_fs):
    if fs > target_fs:
        q = int(round(fs / target_fs))
        return decimate(x, q), fs / q
    return np.asarray(x, dtype=float), float(fs)


def analyze_window(window, fs, use_vmd, K, analysis_fs, t_start):
    """Entropy analysis of ONE window: optional VMD at the full rate, then every
    entropy measure for the raw signal (and each IMF) at a fixed analysis rate."""
    t0 = time.perf_counter()
    x = np.asarray(window, dtype=float)
    x = x[: len(x) // 2 * 2]
    x = x - x.mean()

    names, series = ["Raw"], [x]
    centers = mode_rms = None
    n_iter = 0
    if use_vmd:
        modes, om, n_iter = vmd_decompose(x, K)
        centers = om * fs
        mode_rms = np.sqrt(np.mean(modes ** 2, axis=1))
        names += [f"IMF{k + 1}" for k in range(K)]
        series += list(modes)

    dec, fs_dec = [], float(fs)
    for s in series:
        d, fs_dec = decimate_if_needed(s, fs, analysis_fs)
        dec.append(d)

    matrix = np.full((len(ENTROPY_DEFS), len(dec)), np.nan)
    for j, s in enumerate(dec):
        vals = compute_all_entropies(s, fs_dec)
        for i, (key, _, _) in enumerate(ENTROPY_DEFS):
            matrix[i, j] = vals[key]

    mse = multiscale_sample_entropy(dec[0], MSE_MAX_SCALE, EN_M)

    # decimation's low-pass cuts off at ~0.8 x Nyquist = 0.4 x fs_dec
    in_band = np.ones(len(dec), dtype=bool)
    if use_vmd and fs_dec < fs:
        in_band[1:] = centers <= 0.4 * fs_dec

    return dict(t_start=t_start, t_end=t_start + len(x) / fs, fs=fs, fs_dec=fs_dec,
                use_vmd=use_vmd, K=K if use_vmd else 0, names=names, x=x,
                modes=(modes if use_vmd else None), centers=centers, mode_rms=mode_rms,
                raw_rms=float(np.sqrt(np.mean(x ** 2))), matrix=matrix, mse=mse,
                mse_ci=float(np.nansum(mse)), in_band=in_band, n_iter=n_iter,
                elapsed=time.perf_counter() - t0)


def result_to_row(res, label):
    row = {"Timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "Label": label,
           "WindowStart_s": round(res["t_start"], 3), "WindowEnd_s": round(res["t_end"], 3),
           "AcqFs": res["fs"], "AnalysisFs": res["fs_dec"], "K": res["K"],
           "RawRMS": res["raw_rms"]}
    for i, (key, _, _) in enumerate(ENTROPY_DEFS):
        for j, name in enumerate(res["names"]):
            row[f"{key}_{name}"] = res["matrix"][i, j]
    for s, v in enumerate(res["mse"], start=1):
        row[f"MSE_scale{s}"] = v
    row["MSE_ComplexityIndex"] = res["mse_ci"]
    if res["use_vmd"]:
        for k in range(res["K"]):
            row[f"IMF{k + 1}_CenterFreqHz"] = res["centers"][k]
            row[f"IMF{k + 1}_RMS"] = res["mode_rms"][k]
            row[f"IMF{k + 1}_InBand"] = bool(res["in_band"][k + 1])
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
    and runs the entropy analysis on every ANALYSIS_SECONDS of playback."""

    def __init__(self, data, fs, xmax=None, speed=1.0, analysis=True,
                 analysis_seconds=ANALYSIS_SECONDS, use_vmd=USE_VMD, vmd_k=VMD_K,
                 analysis_fs=ANALYSIS_FS, label="Unlabelled", log_file=None):
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

        # ---- entropy analysis settings / state ----
        self.win = int(round(analysis_seconds * fs))
        self.analysis = bool(analysis) and self.win <= self.total
        if analysis and not self.analysis:
            print("Dataset is shorter than the analysis window - entropy analysis disabled.")
        self.use_vmd, self.vmd_k, self.analysis_fs = use_vmd, vmd_k, analysis_fs
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

        info = (f"Entropy analysis: first result after {analysis_seconds:g} s of data..."
                if self.analysis else "Entropy analysis: disabled")
        self.info = self.fig.text(0.5, 0.012, info, ha="center", fontsize=10, color="navy")
        self.fig.tight_layout(rect=(0, 0.04, 1, 0.94))
        self.render(0)

        if self.analysis:
            self.fig2 = plt.figure(figsize=(16, 8.8))
            self.fig2.text(0.5, 0.5, f"Waiting for the first {analysis_seconds:g} s of data...\n"
                           "(entropy analysis results will appear here)",
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

    # ---- entropy analysis: background job, drawn from the main thread ----
    def _maybe_start_analysis(self, n, start_n):
        delivered = n - start_n
        idx = delivered // self.win
        if idx > self.last_idx and self.pending is None and n >= self.win:
            self.last_idx = idx
            window = self.stream[n - self.win:n].copy()
            t_start = (n - self.win) / self.fs
            self.pending = self.executor.submit(analyze_window, window, self.fs, self.use_vmd,
                                                self.vmd_k, self.analysis_fs, t_start)
            self.info.set_text(f"Entropy analysis: analysing {t_start:.0f}-{n / self.fs:.0f} s ...")

    def _poll_analysis(self):
        if self.pending is not None and self.pending.done():
            fut, self.pending = self.pending, None
            try:
                res = fut.result()
            except Exception as exc:                       # keep the monitor alive
                print(f"Entropy analysis failed: {exc}")
                self.info.set_text(f"Entropy analysis: failed ({exc})")
                return
            self.show_analysis(res)

    def show_analysis(self, res):
        self.n_results += 1
        raw_vals = {key: res["matrix"][i, 0] for i, (key, _, _) in enumerate(ENTROPY_DEFS)}
        summary = " | ".join(f"{k} {v:.3f}" for k, v in raw_vals.items())
        print(f"[analysis #{self.n_results}] window {res['t_start']:.0f}-{res['t_end']:.0f} s | "
              f"RAW: {summary} | MSE-CI {res['mse_ci']:.2f} | {res['elapsed']:.1f} s")
        self.info.set_text(f"Entropy #{self.n_results} (raw, {res['t_start']:.0f}-{res['t_end']:.0f} s):  "
                           + "   ".join(f"{k} {v:.3f}" for k, v in raw_vals.items()))

        if self.log_file:
            self.rows.append(result_to_row(res, self.label))
            try:
                pd.DataFrame(self.rows).to_excel(self.log_file, index=False)
            except Exception as exc:
                print(f"Could not write {self.log_file}: {exc}")

        self._draw_analysis_figure(res)

    def _draw_analysis_figure(self, res):
        fig = self.fig2
        fig.clf()
        fs, names, matrix = res["fs"], res["names"], res["matrix"]
        n_ent, n_ser = matrix.shape
        in_band = res["in_band"]
        ent_labels = [lab for _, lab, _ in ENTROPY_DEFS]

        if res["use_vmd"]:
            K = res["K"]
            gs = fig.add_gridspec(K + 1, 3, width_ratios=[1.2, 1.7, 1.0], hspace=1.0, wspace=0.32,
                                  left=0.045, right=0.985, top=0.84, bottom=0.075)
            t = np.arange(len(res["x"])) / fs
            for k in range(K):
                ax = fig.add_subplot(gs[k, 0])
                ax.plot(t, res["modes"][k], color="blue" if in_band[k + 1] else "gray", lw=0.6)
                star = "" if in_band[k + 1] else " *"
                ax.set_title(f"IMF {k + 1}  ({res['centers'][k]:.1f} Hz, RMS {res['mode_rms'][k]:.4f} g){star}",
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

            # --- heat-map table: rows = entropy measures, columns = Raw + IMFs ---
            ax_h = fig.add_subplot(gs[:, 1])
            norm_m = np.full_like(matrix, np.nan)
            for i in range(n_ent):
                row = matrix[i]
                ok = in_band & np.isfinite(row)
                if ok.sum() >= 1:
                    lo, hi = row[ok].min(), row[ok].max()
                    norm_m[i] = (row - lo) / (hi - lo) if hi > lo else 0.5
            ax_h.imshow(np.clip(np.nan_to_num(norm_m, nan=0.0), 0, 1), aspect="auto",
                        cmap="YlOrRd", vmin=0, vmax=1)
            for i in range(n_ent):
                for j in range(n_ser):
                    v = matrix[i, j]
                    ax_h.text(j, i, "n/a" if np.isnan(v) else f"{v:.3f}", ha="center", va="center",
                              fontsize=8, color="black" if in_band[j] else "dimgray")
            for j in range(n_ser):
                if not in_band[j]:
                    ax_h.add_patch(Rectangle((j - 0.5, -0.5), 1, n_ent, fill=True, facecolor="white",
                                             alpha=0.55, hatch="///", edgecolor="gray", linewidth=0))
            xl = [("Raw" if j == 0 else f"IMF{j}\n{res['centers'][j - 1]:.0f} Hz")
                  + ("" if in_band[j] else " *") for j in range(n_ser)]
            ax_h.set_xticks(range(n_ser))
            ax_h.set_xticklabels(xl, fontsize=8)
            ax_h.set_yticks(range(n_ent))
            ax_h.set_yticklabels(ent_labels, fontsize=9)
            ax_h.xaxis.tick_top()
            ax_h.set_title("Entropy of the raw signal and of every VMD mode\n"
                           "(colour = relative level within each row; numbers = actual values)",
                           fontsize=10, pad=30)
            fig.text(0.28, 0.02,
                     f"* hatched = mode centre frequency above the analysis band (~{0.4 * res['fs_dec']:.0f} Hz at "
                     f"{res['fs_dec']:g} Hz analysis rate): its entropies describe filtered-out content - ignore.",
                     fontsize=8, color="dimgray")

            gs_r = gs[:, 2].subgridspec(2, 1, hspace=0.6)
            ax_m = fig.add_subplot(gs_r[0])
            ax_b = fig.add_subplot(gs_r[1])
            bar_colors = ["tab:blue" if ok else "lightgray" for ok in in_band[1:]]
            ax_b.bar(res["centers"], res["mode_rms"], width=0.02 * fs / 2, color=bar_colors,
                     edgecolor="black", linewidth=0.5)
            ax_b.set_xlabel("Mode centre frequency (Hz)", fontsize=9)
            ax_b.set_ylabel("Mode RMS (g)", fontsize=9)
            ax_b.set_title("VMD band RMS", fontsize=10)
            ax_b.grid(True, alpha=0.4)
        else:
            gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1.0], wspace=0.3,
                                  left=0.16, right=0.97, top=0.88, bottom=0.12)
            ax_h = fig.add_subplot(gs[0, 0])
            vals = matrix[:, 0]
            ax_h.barh(range(n_ent), np.nan_to_num(vals), color="tab:blue", edgecolor="black", linewidth=0.5)
            for i, v in enumerate(vals):
                ax_h.text(max(v, 0) + 0.01, i, "n/a" if np.isnan(v) else f"{v:.3f}", va="center", fontsize=9)
            ax_h.set_yticks(range(n_ent))
            ax_h.set_yticklabels(ent_labels, fontsize=9)
            ax_h.invert_yaxis()
            ax_h.set_xlim(0, max(1.1, float(np.nanmax(vals)) * 1.2))
            ax_h.set_xlabel("Entropy value")
            ax_h.set_title("Entropy of the raw signal")
            ax_h.grid(True, alpha=0.4, axis="x")
            ax_m = fig.add_subplot(gs[0, 1])

        scales = np.arange(1, len(res["mse"]) + 1)
        ax_m.plot(scales, res["mse"], "o-", color="darkgreen", lw=1.5)
        ax_m.set_xlabel("Scale factor", fontsize=9)
        ax_m.set_ylabel("Sample Entropy", fontsize=9)
        ax_m.set_title(f"Multiscale Sample Entropy (raw)\ncomplexity index = {res['mse_ci']:.2f}", fontsize=10)
        ax_m.set_xticks(scales)
        ax_m.grid(True, alpha=0.4)

        mode_txt = f"VMD K={res['K']} + " if res["use_vmd"] else ""
        fig.suptitle(
            f"{mode_txt}All-entropy analysis  -  #{self.n_results}  -  window {res['t_start']:.0f}-"
            f"{res['t_end']:.0f} s of playback  -  {self.label}   [computed in {res['elapsed']:.1f} s, "
            f"entropies at {res['fs_dec']:g} Hz]", fontsize=12, fontweight="bold", y=0.985)
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


# ===================== static (still) plots =====================
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
    global PERM_ORDER
    ap = argparse.ArgumentParser(description="Playback simulation with time domain, FFT and all-entropy analysis.")
    ap.add_argument("--file", default=DATASET_FILE, help="Excel dataset to load")
    ap.add_argument("--speed", type=float, default=1.0, help="playback speed (1 = real time, 5 min)")
    ap.add_argument("--start-group", type=int, default=1, help="group to start playing from (1-based)")
    ap.add_argument("--xmax", type=float, default=None, help="zoom the FFT x-axis to this max frequency (Hz)")
    ap.add_argument("--interval", type=int, default=100, help="screen refresh interval in ms")
    ap.add_argument("--no-analysis", action="store_true", help="disable the entropy analysis entirely")
    ap.add_argument("--no-vmd", action="store_true", help="entropies of the raw signal only (no VMD modes)")
    ap.add_argument("--analysis-seconds", type=float, default=ANALYSIS_SECONDS,
                    help="run the entropy analysis every this many seconds of playback")
    ap.add_argument("--vmd-k", type=int, default=VMD_K, help="number of VMD modes")
    ap.add_argument("--analysis-fs", type=float, default=ANALYSIS_FS,
                    help="rate (Hz) the entropies are computed at (default 500; see ANALYSIS BAND WARNING)")
    ap.add_argument("--perm-order", type=int, default=PERM_ORDER, help="permutation entropy order")
    ap.add_argument("--label", default="Unlabelled", help="condition label stored in the log (e.g. Healthy)")
    ap.add_argument("--log", default=None, help="Excel file to write the entropy results to")
    ap.add_argument("--static", action="store_true", help="draw still plots instead of the simulation")
    ap.add_argument("--group", type=int, default=1, help="(static mode) group number to plot")
    ap.add_argument("--save", action="store_true", help="(static mode) save PNGs instead of opening windows")
    args = ap.parse_args()
    PERM_ORDER = args.perm_order

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
                            use_vmd=not args.no_vmd, vmd_k=args.vmd_k, analysis_fs=args.analysis_fs,
                            label=args.label, log_file=args.log)
    remaining = (data.shape[1] - args.start_group + 1) * data.shape[0] / fs / args.speed
    print(f"Starting playback - about {fmt_clock(remaining)} (min:sec) at x{args.speed:g} speed. "
          f"SPACE = pause/resume, R = restart.")
    if sim.analysis:
        what = f"VMD (K={args.vmd_k}) + 7 entropies" if not args.no_vmd else "7 entropies (raw signal)"
        print(f"{what} every {args.analysis_seconds:g} s of playback, shown in a second window "
              f"(entropies computed at {args.analysis_fs:g} Hz).")
    sim.run(start_group=args.start_group, interval_ms=args.interval)


if __name__ == "__main__":
    main()
