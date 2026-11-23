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

BEARING VIEW ON YOUR DATASET (--bearing)  [plays back YOUR recording like a replay; no random simulation]
  Your original replay (time domain, FFT with Mode labels) is kept. In addition:
    LIVE WINDOW
    * TIME DOMAIN: the four fault components - BPFO, BPFI, BSF, FTF - extracted from your real signal
      (envelope of the resonance band, narrow band-pass at each fault's lines) and their COMBINED sum
    * FFT: BPFO, BPFI, 2xBSF and FTF marked on your FFT; the FFT of the envelope shows the noise floor, the
      fault threshold, (with --reference) the Healthy spectrum, and the measured peaks of the winning fault
      with their dB value
    * WHAT DEFINES THE FAULT: a decision panel with the evidence of each fault family against the threshold,
      and a sentence saying why. A fault is declared when its strongest lines stand >= 8 dB above the local
      noise floor (8 dB = just above the strongest value seen in your two Healthy recordings, 7.1 dB)
    SECOND WINDOW (every 10 s of replay)
    * COMBINED fault signal -> VMD -> IMFs. Each IMF is aligned with a fault type through its correlation with
      the four fault components (heat-map + red box), then the seven entropies of the combined signal and of
      every IMF are computed = the fault feature vector (Data -> VMD -> entropy, as in the project flow)
  FTF and BSF follow from the cage frequency (measured from your Faulty recordings: 9.52 Hz). BPFO = n x FTF
  and BPFI = n x (fr - FTF) need the BALL COUNT n, which your data cannot tell: they are drawn for an assumed
  n (marked *) and excluded from the verdict until you pass --balls N.
  Your original raw-signal VMD + entropy window is off in this mode; --keep-base-analysis brings it back.
  Examples:  python <this file> --file Vibration_Dataset_Faulty1_2000Hz.xlsx --bearing --reference Vibration_Dataset_Healthy_2000Hz.xlsx
             python <this file> --file Vibration_Dataset_Healthy2_2000Hz.xlsx --bearing --speed 10
             python <this file> --file Faulty.xlsx --bearing --balls 9 --log results.xlsx

BEARING SIMULATION MODE (--bearing-sim)  [added; the dataset mode above is unchanged]
  Simulates a rolling-element bearing (SKF 6205-2RS JEM, 1800 rpm) whose signal contains all four
  defect signals at once and plays it back live:
    BPFO outer race (impacts at BPFO) | BPFI inner race (impacts at BPFI, modulated at shaft speed)
    BSF ball (impacts at 2xBSF, modulated at the cage speed) | FTF cage (FTF sines + impacts at FTF)
  Each impact rings a 3.5 kHz resonance; 1 % random slip, shaft harmonics and noise are added.
  Live window: the four fault signals and the COMBINED signal in the time domain; the raw FFT
  (the four frequencies carry no energy there - the energy sits at the resonance); and the FFT of
  the ENVELOPE with BPFO, BPFI, BSF/2xBSF and FTF marked - that is where they show up. A banner
  gives the diagnosed condition (default scenario: "BPFO FAULTY"). A second window (every 10 s)
  runs VMD on the demodulated signal, labels each mode with the fault line it carries and ranks
  the modes by fault relevance (IMF1 = strongest fault mode = BPFO), with the entropy table.
  Examples:  python <this file> --bearing-sim
             python <this file> --bearing-sim --dominant BPFI --speed 5
             python <this file> --bearing-sim --healthy
             python <this file> --bearing-sim --weights BPFO=1,BPFI=1,BSF=1,FTF=3 --components

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
from scipy.signal import butter, filtfilt, find_peaks, decimate, sosfiltfilt, hilbert, fftconvolve
from scipy.signal.windows import hann
from scipy.spatial import cKDTree
from scipy.stats import norm
from scipy.ndimage import median_filter

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


# ===================== BEARING FAULT SIMULATION (new: --bearing-sim) =====================
# Everything above this line is your original code, unchanged. This section simulates a bearing signal
# that contains BPFO, BPFI, BSF and FTF at once, and plays it back live (time domain, FFT, VMD, entropy).
# ---------------- BEARING SIMULATION SETTINGS (new) ----------------
FS_SIM        = 12800            # sampling rate of the simulated signal (Hz)
BEARING       = dict(n=9, d=7.94, D=39.04, phi=0.0)   # SKF 6205-2RS JEM: balls, ball dia, pitch dia (mm), contact angle (deg)
RPM           = 1800.0           # shaft speed
RESONANCE_HZ  = 3500.0           # structural resonance excited by the impacts
BAND          = (2500.0, 4500.0) # demodulation band around the resonance
DAMPING       = 0.05
SLIP          = 0.01             # 1 % random slip
NOISE_STD     = 0.05             # g
DURATION_S    = 60.0
# fault severity presets (impact-amplitude weights). FTF needs a larger weight than the others to be seen at all:
# a 12 Hz impact train puts far less energy into the envelope spectrum than a 107-162 Hz one.
PRESETS = {"BPFO": dict(BPFO=1.0, BPFI=0.3, BSF=0.5, FTF=1.5),
           "BPFI": dict(BPFO=0.2, BPFI=1.0, BSF=0.5, FTF=1.5),
           "BSF":  dict(BPFO=0.2, BPFI=0.3, BSF=1.5, FTF=1.5),
           "FTF":  dict(BPFO=0.3, BPFI=0.3, BSF=0.5, FTF=6.0)}
DEFAULT_WEIGHTS = PRESETS["BPFO"]

DISPLAY_SECONDS  = 0.5           # width of the rolling time-domain display
SPEC_SECONDS     = 5.0           # data used for the live FFT / envelope spectrum
SPEC_EVERY       = 0.5           # spectra refresh interval (s of playback)
BRG_ANALYSIS_SECONDS = 10.0          # VMD + entropy analysis interval
BRG_VMD_K        = 10            # number of VMD modes
VMD_HALFWIDTH_HZ = 5.0           # VMD mode bandwidth (alpha is derived from it)
ENV_FS_TARGET    = 800.0         # sampling rate of the demodulated envelope
ENV_LOWPASS_HZ   = 300.0         # envelope low-pass used before VMD (keeps the fault lines, drops the rest)
LIVE_ENV_FS      = 1600.0        # envelope sampling rate of the live spectra (no extra low-pass)
N_ENTROPY_MODES  = 6             # entropies are computed for the envelope + this many top modes
FAULT_DB         = 12.0          # line strength above the noise floor that declares a fault
SECONDARY_DB     = 9.0           # ... and that lists a secondary component
FAULT_COLORS = {"BPFO": "tab:red", "BPFI": "tab:orange", "BSF": "tab:purple", "FTF": "tab:green"}
# ----------------------------------------------------------

BRG_FUZZY_ROWS    = 1500          # FuzzyEn template rows in the bearing analysis (keeps a 10 s analysis fast)

# ===================== bearing physics =====================
def bearing_frequencies(n, d, D, phi_deg, rpm):
    """Characteristic frequencies (Hz): outer race passes, inner race passes, ball spin, cage."""
    fr = rpm / 60.0
    r = d / D * np.cos(np.radians(phi_deg))
    return {"fr": fr,
            "FTF": fr / 2 * (1 - r),
            "BPFO": n * fr / 2 * (1 - r),
            "BPFI": n * fr / 2 * (1 + r),
            "BSF": D / (2 * d) * fr * (1 - r ** 2)}


def impulse_response(fs, fn=RESONANCE_HZ, zeta=DAMPING, length=0.02):
    t = np.arange(0, length, 1 / fs)
    wn = 2 * np.pi * fn
    wd = wn * np.sqrt(1 - zeta ** 2)
    return np.exp(-zeta * wn * t) * np.sin(wd * t)


def simulate_bearing(duration, fs, F, weights, seed=42, noise_std=NOISE_STD,
                     resonance=RESONANCE_HZ, zeta=DAMPING, slip=SLIP):
    """Bearing vibration (g) containing BPFO, BPFI, BSF and FTF components.
    Returns (t, combined signal, dict of the individual components)."""
    rng = np.random.default_rng(seed)
    n = int(round(duration * fs))
    t = np.arange(n) / fs
    fr, ftf = F["fr"], F["FTF"]
    ir = impulse_response(fs, resonance, zeta)

    def train(fd, mod=None):
        x = np.zeros(n)
        T = 1.0 / fd
        tk = rng.uniform(0, T)
        while tk < (n - 1) / fs:
            x[int(round(tk * fs))] += 1.0 if mod is None else mod(tk)
            tk += T * (1 + slip * rng.standard_normal())       # random slip on the timing
        return fftconvolve(x, ir)[:n]                          # each impact rings the resonance

    shaft = (0.20 * np.sin(2 * np.pi * fr * t) + 0.08 * np.sin(2 * np.pi * 2 * fr * t + 0.5)
             + 0.03 * np.sin(2 * np.pi * 3 * fr * t + 1.1))
    comps = {"baseline": shaft + noise_std * rng.standard_normal(n)}
    w = weights
    zero = np.zeros(n)
    comps["BPFO"] = w.get("BPFO", 0) * 1.0 * train(F["BPFO"]) if w.get("BPFO", 0) else zero
    comps["BPFI"] = (w["BPFI"] * 0.8 * train(F["BPFI"], lambda tk: 0.5 * (1 + np.cos(2 * np.pi * fr * tk)) + 0.1)
                     if w.get("BPFI", 0) else zero)
    comps["BSF"] = (w["BSF"] * 0.6 * train(2 * F["BSF"], lambda tk: 0.5 * (1 + np.cos(2 * np.pi * ftf * tk)) + 0.1)
                    if w.get("BSF", 0) else zero)
    comps["FTF"] = (w["FTF"] * (0.10 * np.sin(2 * np.pi * ftf * t) + 0.05 * np.sin(2 * np.pi * 2 * ftf * t + 0.3)
                                + 0.4 * train(ftf)) if w.get("FTF", 0) else zero)
    x = sum(comps.values())
    return t, x, comps


# ===================== demodulation (envelope) =====================
def decimate_stages(x, q):
    """Decimate by q in stages of at most 8 (scipy recommends this for large q)."""
    while q > 1:
        for s in (8, 7, 6, 5, 4, 3, 2):
            if q % s == 0:
                x = decimate(x, s)
                q //= s
                break
        else:
            break
    return x


def demodulate(x, fs, band=BAND, fe_target=ENV_FS_TARGET, lowpass=ENV_LOWPASS_HZ):
    """band-pass around the resonance -> Hilbert envelope -> decimate -> low-pass.
    Returns (envelope, its sampling rate)."""
    x = np.asarray(x, dtype=float)
    hi = min(band[1], 0.45 * fs)
    sos = butter(4, [band[0] / (fs / 2), hi / (fs / 2)], btype="band", output="sos")
    env = np.abs(hilbert(sosfiltfilt(sos, x - x.mean())))
    q = max(1, int(round(fs / fe_target)))
    env = decimate_stages(env - env.mean(), q)
    fe = fs / q
    if lowpass and lowpass < 0.45 * fe:
        env = sosfiltfilt(butter(6, lowpass / (fe / 2), output="sos"), env)
    return env - env.mean(), fe


# ===================== fault lines + diagnosis =====================
def fault_lines(F, fmax=450.0):
    """Where each fault shows up in the ENVELOPE spectrum: {fault: [(name, Hz), ...]}."""
    fr, ftf, bpfo, bpfi, bsf = F["fr"], F["FTF"], F["BPFO"], F["BPFI"], F["BSF"]
    L = {"BPFO": [(f"{k}xBPFO" if k > 1 else "BPFO", k * bpfo) for k in range(1, 5)],
         "BPFI": [("BPFI", bpfi), ("BPFI-fr", bpfi - fr), ("BPFI+fr", bpfi + fr), ("2xBPFI", 2 * bpfi)],
         "BSF": [("2xBSF", 2 * bsf), ("2xBSF-FTF", 2 * bsf - ftf), ("2xBSF+FTF", 2 * bsf + ftf), ("4xBSF", 4 * bsf)],
         "FTF": [(f"{k}xFTF" if k > 1 else "FTF", k * ftf) for k in range(1, 5)]}
    return {k: [(n, f) for n, f in v if f < fmax] for k, v in L.items()}


def line_snr_db(f, a, f0, tol=0.6, guard=1.5, bw=15.0):
    """Peak at f0 (+/- tol) in dB above the median spectrum level around it."""
    w = (f > f0 - tol) & (f < f0 + tol)
    if not w.any():
        return 0.0
    nb = (f > f0 - bw) & (f < f0 + bw) & ~((f > f0 - guard) & (f < f0 + guard))
    noise = np.median(a[nb]) + 1e-15 if nb.any() else np.median(a) + 1e-15
    return float(20 * np.log10(max(a[w].max(), 1e-15) / noise))


def diagnose(f, a, F):
    """Evidence for each fault = mean of its two strongest lines (dB above the noise floor)."""
    ev = {}
    for fault, lines in fault_lines(F).items():
        snr = sorted([line_snr_db(f, a, fl) for _, fl in lines], reverse=True)
        ev[fault] = float(np.mean(snr[:2]))
    return ev


def condition_from_evidence(ev):
    """-> (verdict text, primary fault or None, list of secondary faults)."""
    ranked = sorted(ev.items(), key=lambda kv: -kv[1])
    top, topv = ranked[0]
    if topv < FAULT_DB:
        return "HEALTHY", None, []
    others = [k for k, v in ranked[1:] if v >= SECONDARY_DB]
    return f"{top} FAULTY", top, others


def label_modes(modes, fe, F, tol=0.6):
    """Label each VMD mode with the fault line that carries the largest peak of ITS spectrum
    (needs to be at least 25 % of the mode's strongest peak). -> list of (fault, line name)."""
    lines = [(fault, name, fl) for fault, ls in fault_lines(F).items() for name, fl in ls]
    lines += [("SHAFT", "fr", F["fr"]), ("SHAFT", "2xfr", 2 * F["fr"])]
    out = []
    for m in modes:
        sp = np.abs(np.fft.rfft(m * np.hanning(len(m)))) / len(m) * 2
        fq = np.fft.rfftfreq(len(m), 1 / fe)
        best, bestv = None, 0.0
        for fault, name, fl in lines:
            w = (fq > fl - tol) & (fq < fl + tol)
            if w.any() and sp[w].max() > bestv:
                best, bestv = (fault, name), sp[w].max()
        out.append(best if best and bestv >= 0.25 * sp[fq > 3].max() else (None, None))
    return out


# ===================== VMD + entropy analysis of one 10 s window =====================
def compute_entropies_bearing(series, fs):
    """The same seven entropies as compute_all_entropies(), but FuzzyEn uses fewer template rows so that
    a 10 s analysis finishes in a few seconds."""
    series = np.asarray(series, dtype=float)
    if len(series) < 50 or not np.isfinite(series).all() or _std(series) < 1e-15:
        return {key: np.nan for key, _, _ in ENTROPY_DEFS}
    out = {}
    with np.errstate(all="ignore"):
        for key, _, fn in ENTROPY_DEFS:
            try:
                out[key] = float(fuzzy_entropy(series, max_rows=BRG_FUZZY_ROWS) if key == "FuzzyEn" else fn(series, fs))
            except Exception:
                out[key] = np.nan
    return out


def analyze_window_bearing(window, fs, F, band, K, halfwidth_hz, n_modes, t_start):
    t0 = time.perf_counter()
    env, fe = demodulate(window, fs, band)
    f, a = welch_amplitude_spectrum(env, fe, 3)
    evidence = diagnose(f, a, F)

    modes, om, n_iter = vmd_decompose(env, K, alpha=(fe / halfwidth_hz) ** 2)
    centers = om * fe
    rms = np.sqrt(np.mean(modes ** 2, axis=1))
    labels = label_modes(modes, fe, F)
    # only keep a fault label if that fault is really detected in the envelope spectrum
    # (otherwise a healthy bearing would still get its biggest noise peak labelled as a "fault")
    labels = [lab if (lab[0] == "SHAFT" or (lab[0] and evidence[lab[0]] >= SECONDARY_DB)) else (None, None)
              for lab in labels]
    # rank: fault-labelled modes first, then by RMS (strongest first)
    order = sorted(range(K), key=lambda i: (labels[i][0] in (None, "SHAFT"), -rms[i]))
    modes, centers, rms = modes[order], centers[order], rms[order]
    labels = [labels[i] for i in order]

    names = ["Env"] + [f"IMF{k + 1}" for k in range(min(n_modes, K))]
    series = [env] + [modes[k] for k in range(min(n_modes, K))]
    matrix = np.full((len(ENTROPY_DEFS), len(series)), np.nan)
    for j, s in enumerate(series):
        vals = compute_entropies_bearing(s, fe)
        for i, (key, _, _) in enumerate(ENTROPY_DEFS):
            matrix[i, j] = vals[key]
    mse = multiscale_sample_entropy(env, MSE_MAX_SCALE, EN_M)

    return dict(t_start=t_start, t_end=t_start + len(window) / fs, fs=fs, fe=fe, K=K, env=env,
                spec_f=f, spec_a=a, evidence=evidence, modes=modes, centers=centers, rms=rms,
                labels=labels, names=names, matrix=matrix, mse=mse, mse_ci=float(np.nansum(mse)),
                n_iter=n_iter, elapsed=time.perf_counter() - t0)


def bearing_result_to_row(res, label, verdict):
    row = {"Timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "Label": label,
           "WindowStart_s": round(res["t_start"], 3), "WindowEnd_s": round(res["t_end"], 3),
           "Condition": verdict, "AcqFs": res["fs"], "EnvFs": res["fe"], "K": res["K"]}
    for k, v in res["evidence"].items():
        row[f"Evidence_{k}_dB"] = round(v, 2)
    for k in range(res["K"]):
        fault, line = res["labels"][k]
        row[f"IMF{k + 1}_CenterFreqHz"] = res["centers"][k]
        row[f"IMF{k + 1}_RMS"] = res["rms"][k]
        row[f"IMF{k + 1}_Fault"] = f"{fault}:{line}" if fault else ""
    for i, (key, _, _) in enumerate(ENTROPY_DEFS):
        for j, name in enumerate(res["names"]):
            row[f"{key}_{name}"] = res["matrix"][i, j]
    for s, v in enumerate(res["mse"], start=1):
        row[f"MSE_scale{s}"] = v
    row["MSE_ComplexityIndex"] = res["mse_ci"]
    return row


def mark_fault_lines(ax, F, xmax, ymax=None, label_size=8, show_text=True, note=None):
    """Draw the BPFO / BPFI / BSF / FTF markers on an envelope-spectrum axis."""
    note = note or {}
    fr = F["fr"]
    tr = ax.get_xaxis_transform()                    # x in data units, y in axes fraction
    spec = [("BPFO", F["BPFO"], "BPFO", "-"), ("BPFO", 2 * F["BPFO"], "2xBPFO", "--"),
            ("BPFO", 3 * F["BPFO"], "3xBPFO", "--"),
            ("BPFI", F["BPFI"], "BPFI", "-"), ("BPFI", F["BPFI"] - fr, "BPFI-fr", ":"),
            ("BPFI", F["BPFI"] + fr, "BPFI+fr", ":"),
            ("BSF", 2 * F["BSF"], "2xBSF", "-"), ("BSF", F["BSF"], "BSF", "--"),
            ("FTF", F["FTF"], "FTF", "-"), ("FTF", 2 * F["FTF"], "2xFTF", "--")]
    for fault, fq, name, ls in spec:
        if fq < xmax:
            ax.axvline(fq, color=FAULT_COLORS[fault], ls=ls, lw=1.2 if ls == "-" else 0.9, alpha=0.85)
            if show_text and ls != ":":
                ax.text(fq, 0.98 if ls == "-" else 0.62, f" {name} {fq:.1f}{note.get(fault, '')}", transform=tr,
                        color=FAULT_COLORS[fault], rotation=90, va="top", ha="right", fontsize=label_size)
    ax.axvline(fr, color="gray", ls="-.", lw=0.8)
    if show_text:
        ax.text(fr, 0.40, f" fr {fr:.1f}", transform=tr, color="dimgray", rotation=90, va="top", ha="right",
                fontsize=label_size)


def log_limits(a, f=None, lo_div=3.0, hi_mul=3.0, floor=1e-9):
    """y-limits for a log-scale spectrum: a little below the median, well above the peak."""
    a = np.asarray(a, float) if f is None else np.asarray(a, float)[f]
    a = a[a > 0]
    return max(np.median(a) / lo_div, floor), a.max() * hi_mul


def smooth_log(old, new, k=0.2):
    return 10 ** ((1 - k) * np.log10(old) + k * np.log10(new))


# ===================== live simulator =====================
class BearingSimulator:
    """Plays the simulated bearing signal back live: time domain, raw FFT, envelope spectrum with the
    four fault markers, live diagnosis banner and (every 10 s) the VMD + entropy analysis window."""

    def __init__(self, x, fs, F, comps=None, speed=1.0, band=BAND, display_seconds=DISPLAY_SECONDS,
                 spec_seconds=SPEC_SECONDS, spec_every=SPEC_EVERY, analysis=True,
                 analysis_seconds=BRG_ANALYSIS_SECONDS, vmd_k=BRG_VMD_K, label="Simulated bearing",
                 log_file=None, xmax_env=450.0):
        self.x, self.fs, self.F, self.comps = np.asarray(x, float), fs, F, comps
        self.total = len(self.x)
        self.speed, self.band, self.label, self.log_file = speed, band, label, log_file
        self.dw = int(round(display_seconds * fs))
        self.sw = int(round(spec_seconds * fs))
        self.spec_step = max(1, int(round(spec_every * fs)))
        self.last_spec_idx = -1
        self.xmax_env = xmax_env
        self.verdict_now = "ACQUIRING"
        self.ani = None
        self._tick = None
        self.vmd_k = vmd_k
        self.win = int(round(analysis_seconds * fs))
        self.analysis = bool(analysis) and self.win <= self.total
        self.executor = ThreadPoolExecutor(max_workers=1) if self.analysis else None
        self.pending, self.last_idx, self.rows, self.n_results, self.fig2 = None, 0, [], 0, None
        self.y_lim, self.raw_lim, self.env_lim = 0.5, (1e-4, 1.0), (1e-4, 0.5)
        self.fault_keys = ["BPFO", "BPFI", "BSF", "FTF"]
        self.comp_present, self.comp_scale = {}, {}
        for k in self.fault_keys:
            peak = float(np.max(np.abs(comps[k]))) if comps is not None else 0.0
            self.comp_present[k] = peak > 1e-9
            self.comp_scale[k] = peak if peak > 1e-9 else 1.0

        self.fig = plt.figure(figsize=(15, 10.2))
        gs = self.fig.add_gridspec(3, 2, height_ratios=[1.3, 1.0, 1.35], left=0.065, right=0.985, top=0.80,
                                   bottom=0.065, hspace=0.50, wspace=0.17)
        self.ax_c = self.fig.add_subplot(gs[0, 0])
        self.ax_t = self.fig.add_subplot(gs[1, 0])
        self.ax_r = self.fig.add_subplot(gs[0:2, 1])
        self.ax_e = self.fig.add_subplot(gs[2, :])
        t_axis = np.arange(1, self.dw + 1) / fs

        # the four fault signals, one trace each (every trace scaled to its own peak so all are readable)
        self.lines_c = {}
        for i, k in enumerate(self.fault_keys):
            self.lines_c[k], = self.ax_c.plot(t_axis, np.full(self.dw, 3.0 - i), color=FAULT_COLORS[k], lw=0.7)
        self.ax_c.set_xlim(t_axis[0], t_axis[-1])
        self.ax_c.set_ylim(-0.6, 3.6)
        self.ax_c.set_yticks([3, 2, 1, 0])
        self.ax_c.set_yticklabels([f"{k}\n" + (f"peak {self.comp_scale[k]:.2f} g" if self.comp_present[k] else "not present")
                                   for k in self.fault_keys], fontsize=8)
        self.ax_c.set_xlabel("Time (s)")
        self.ax_c.set_title("TIME DOMAIN - the four fault signals (each drawn to its own peak)", fontsize=10)
        self.ax_c.grid(True, alpha=0.4)

        # the combined signal
        self.line_t, = self.ax_t.plot(t_axis, np.zeros(self.dw), color="tab:blue", lw=0.6)
        self.ax_t.set_xlim(t_axis[0], t_axis[-1])
        self.ax_t.set_ylim(-self.y_lim, self.y_lim)
        self.ax_t.set_xlabel("Time (s)")
        self.ax_t.set_ylabel("Acceleration (g)")
        self.ax_t.set_title("TIME DOMAIN - COMBINED = baseline + BPFO + BPFI + BSF + FTF", fontsize=10, fontweight="bold")
        self.ax_t.grid(True, alpha=0.4)

        # raw FFT (log-log): shows the resonance and where the four bearing frequencies WOULD be
        self.line_r, = self.ax_r.plot([5, fs / 2], [1e-6, 1e-6], color="red", lw=0.8)
        self.ax_r.axvspan(band[0], min(band[1], fs / 2), color="gold", alpha=0.25,
                          label=f"demodulation band {band[0] / 1000:.1f}-{band[1] / 1000:.1f} kHz")
        self.ax_r.set_xscale("log")
        self.ax_r.set_yscale("log")
        self.ax_r.set_xlim(5, fs / 2)
        self.ax_r.set_ylim(*self.raw_lim)
        tr = self.ax_r.get_xaxis_transform()
        for key, fq, name in [("BPFO", F["BPFO"], "BPFO"), ("BPFI", F["BPFI"], "BPFI"),
                              ("BSF", 2 * F["BSF"], "2xBSF"), ("FTF", F["FTF"], "FTF")]:
            self.ax_r.axvline(fq, color=FAULT_COLORS[key], ls=":", lw=1.3)
            self.ax_r.text(fq, 0.97, f" {name} {fq:.0f}", transform=tr, color=FAULT_COLORS[key], rotation=90,
                           va="top", ha="right", fontsize=8)
        self.ax_r.set_xlabel("Frequency (Hz, log)")
        self.ax_r.set_ylabel("Amplitude (g, log)")
        self.ax_r.set_title("FFT of the raw signal: the four fault frequencies (dotted) carry no energy here,\n"
                            "the impacts ring the 3.5 kHz resonance (see the envelope FFT below)", fontsize=10)
        self.ax_r.legend(loc="lower right", fontsize=8)
        self.ax_r.grid(True, alpha=0.4, which="both")

        # envelope FFT: where the four fault frequencies actually show up
        self.line_e, = self.ax_e.plot([0, xmax_env], [1e-6, 1e-6], color="black", lw=0.9)
        self.ax_e.set_yscale("log")
        self.ax_e.set_xlim(0, xmax_env)
        self.ax_e.set_ylim(*self.env_lim)
        mark_fault_lines(self.ax_e, F, xmax_env)
        self.ax_e.set_xlabel("Frequency (Hz)")
        self.ax_e.set_ylabel("Envelope amplitude (g, log scale)")
        self.ax_e.set_title("FFT OF THE ENVELOPE (demodulated signal) - where BPFO, BPFI, BSF and FTF appear",
                            fontsize=10, fontweight="bold")
        self.ax_e.grid(True, alpha=0.4, which="both")

        self.banner = self.fig.text(0.5, 0.945, "CONDITION: ACQUIRING...", ha="center", fontsize=24,
                                    fontweight="bold", color="gray")
        self.sub = self.fig.text(0.5, 0.895, "", ha="center", fontsize=11, color="dimgray")
        self.clock = self.fig.text(0.5, 0.848, "", ha="center", fontsize=10, color="navy")
        self.info = self.fig.text(0.5, 0.012, "VMD + entropy analysis: first result after "
                                  f"{analysis_seconds:g} s of data..." if self.analysis else "", ha="center",
                                  fontsize=10, color="navy")
        self.render(0)

        if self.analysis:
            self.fig2 = plt.figure(figsize=(16, 9))
            self.fig2.text(0.5, 0.5, f"Waiting for the first {analysis_seconds:g} s of data...\n"
                           "(VMD + entropy analysis will appear here)", ha="center", va="center",
                           fontsize=14, color="gray")
            for fg in (self.fig, self.fig2):
                fg.canvas.mpl_connect("close_event", self._on_close)

    def _on_close(self, _event):
        if self.executor is not None:
            self.executor.shutdown(wait=False, cancel_futures=True)

    # ---- state of the system after n samples have "arrived" ----
    def render(self, n):
        n = int(min(max(n, 0), self.total))
        lo = max(0, n - self.dw)
        buf = np.zeros(self.dw)
        seg = self.x[lo:n]
        if len(seg):
            buf[self.dw - len(seg):] = seg
        self.line_t.set_ydata(buf)
        target = max(1.15 * float(np.max(np.abs(buf))), 0.3)
        self.y_lim = 0.8 * self.y_lim + 0.2 * target
        self.ax_t.set_ylim(-self.y_lim, self.y_lim)
        if self.comps is not None:
            for i, k in enumerate(self.fault_keys):
                cb = np.zeros(self.dw)
                cs = self.comps[k][lo:n]
                if len(cs):
                    cb[self.dw - len(cs):] = cs
                self.lines_c[k].set_ydata((3 - i) + 0.42 * cb / self.comp_scale[k])

        idx = n // self.spec_step
        if idx != self.last_spec_idx or n >= self.total:
            self.last_spec_idx = idx
            self._update_spectra(n)

        total_s = self.total / self.fs
        self.clock.set_text(f"PLAYBACK  t = {fmt_clock(n / self.fs)} / {fmt_clock(total_s)}   |   speed x{self.speed:g}"
                            "   |   [SPACE] pause  [R] restart")

    def _update_spectra(self, n):
        have = min(n, self.sw)
        if have < 2 * self.fs:                                  # need >= 2 s of data
            self.line_r.set_data([0, self.fs / 2], [1e-6, 1e-6])
            self.line_e.set_data([0, self.xmax_env], [1e-6, 1e-6])
            self.banner.set_text("CONDITION: ACQUIRING...")
            self.banner.set_color("gray")
            self.sub.set_text("")
            return
        seg = self.x[n - have:n]
        fr_, ar_ = welch_amplitude_spectrum(seg, self.fs, 6)
        self.line_r.set_data(fr_, np.maximum(ar_, 1e-9))
        lo, hi = log_limits(ar_[fr_ > 5], lo_div=2.0, hi_mul=2.5)
        self.raw_lim = (smooth_log(self.raw_lim[0], lo), smooth_log(self.raw_lim[1], hi))
        self.ax_r.set_ylim(*self.raw_lim)

        env, fe = demodulate(seg, self.fs, self.band, fe_target=LIVE_ENV_FS, lowpass=None)   # wider envelope: shows 3x/4x BPFO too
        fq, ae = welch_amplitude_spectrum(env, fe, 3)
        m = (fq <= self.xmax_env)
        self.line_e.set_data(fq[m], np.maximum(ae[m], 1e-9))
        lo, hi = log_limits(ae[m & (fq > 3)], lo_div=4.0, hi_mul=3.0)
        self.env_lim = (smooth_log(self.env_lim[0], lo), smooth_log(self.env_lim[1], hi))
        self.ax_e.set_ylim(*self.env_lim)
        for t in list(self.ax_e.texts):
            t.remove()
        for ln in list(self.ax_e.lines)[1:]:
            ln.remove()
        mark_fault_lines(self.ax_e, self.F, self.xmax_env)

        ev = diagnose(fq, ae, self.F)
        verdict, primary, others = condition_from_evidence(ev)
        self.verdict_now = verdict
        self.banner.set_text(f"CONDITION: {verdict}")
        self.banner.set_color("tab:red" if primary else (0, 0.55, 0))
        extra = f"   |   also present: {', '.join(others)}" if others else ""
        self.sub.set_text("evidence (dB above noise floor):  " +
                          "   ".join(f"{k} {v:.0f}" for k, v in ev.items()) + extra)

    # ---- VMD + entropy analysis: background job, drawn from the main thread ----
    def _maybe_start_analysis(self, n, start_n=0):
        idx = (n - start_n) // self.win
        if idx > self.last_idx and self.pending is None and n >= self.win:
            self.last_idx = idx
            window = self.x[n - self.win:n].copy()
            t_start = (n - self.win) / self.fs
            self.pending = self.executor.submit(analyze_window_bearing, window, self.fs, self.F, self.band,
                                                self.vmd_k, VMD_HALFWIDTH_HZ, N_ENTROPY_MODES, t_start)
            self.info.set_text(f"VMD + entropy analysis: analysing {t_start:.0f}-{n / self.fs:.0f} s ...")

    def _poll_analysis(self):
        if self.pending is not None and self.pending.done():
            fut, self.pending = self.pending, None
            try:
                res = fut.result()
            except Exception as exc:
                print(f"Analysis failed: {exc}")
                self.info.set_text(f"VMD + entropy analysis failed: {exc}")
                return
            self.show_analysis(res)

    def show_analysis(self, res):
        self.n_results += 1
        verdict, primary, others = condition_from_evidence(res["evidence"])
        k1 = res["labels"][0]
        imf1 = f"{k1[0]} ({res['centers'][0]:.1f} Hz)" if k1[0] else "none"
        print(f"[analysis #{self.n_results}] {res['t_start']:.0f}-{res['t_end']:.0f} s | {verdict} | IMF1 = {imf1} | "
              + " ".join(f"{k} {v:.0f}dB" for k, v in res["evidence"].items()) + f" | {res['elapsed']:.1f} s")
        self.info.set_text(f"VMD + entropy #{self.n_results} ({res['t_start']:.0f}-{res['t_end']:.0f} s):  "
                           f"{verdict}   |   IMF1 = {imf1}")
        if self.log_file:
            self.rows.append(bearing_result_to_row(res, self.label, verdict))
            try:
                pd.DataFrame(self.rows).to_excel(self.log_file, index=False)
            except Exception as exc:
                print(f"Could not write {self.log_file}: {exc}")
        self._draw_analysis_figure(res, verdict, primary, others)

    def _draw_analysis_figure(self, res, verdict, primary, others):
        fig = self.fig2
        fig.clf()
        fe, K = res["fe"], res["K"]
        n_show = min(8, K)
        gs = fig.add_gridspec(n_show, 3, width_ratios=[1.15, 1.6, 1.25], hspace=1.0, wspace=0.3,
                              left=0.045, right=0.985, top=0.74, bottom=0.07)
        zoom = int(0.4 * fe)                                    # first 0.4 s of every mode
        t = np.arange(zoom) / fe
        for k in range(n_show):
            ax = fig.add_subplot(gs[k, 0])
            fault, line = res["labels"][k]
            col = FAULT_COLORS.get(fault, "gray")
            ax.plot(t, res["modes"][k][:zoom], color=col, lw=0.9)
            tag = f"{fault}  [{line}]" if fault and fault != "SHAFT" else ("shaft" if fault == "SHAFT" else "no fault line")
            ax.set_title(f"IMF{k + 1} = {tag}   {res['centers'][k]:.1f} Hz   RMS {res['rms'][k]:.4f}",
                         fontsize=8, pad=2, color=col if fault else "dimgray",
                         fontweight="bold" if k == 0 else "normal")
            ax.set_xlim(0, t[-1])
            ax.tick_params(labelsize=7)
            ax.grid(True, alpha=0.3)
            if k < n_show - 1:
                ax.set_xticklabels([])
        ax.set_xlabel("Time (s) - first 0.4 s of each mode", fontsize=8)

        # ---- entropy heat-map ----
        ax_h = fig.add_subplot(gs[:, 1])
        matrix, names = res["matrix"], res["names"]
        n_ent, n_ser = matrix.shape
        norm_m = np.full_like(matrix, np.nan)
        for i in range(n_ent):
            row = matrix[i]
            ok = np.isfinite(row)
            if ok.sum():
                lo, hi = row[ok].min(), row[ok].max()
                norm_m[i] = (row - lo) / (hi - lo) if hi > lo else 0.5
        ax_h.imshow(np.clip(np.nan_to_num(norm_m, nan=0.0), 0, 1), aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
        for i in range(n_ent):
            for j in range(n_ser):
                v = matrix[i, j]
                ax_h.text(j, i, "n/a" if np.isnan(v) else f"{v:.3f}", ha="center", va="center", fontsize=8)
        xl = ["Envelope"]
        for j in range(1, n_ser):
            fault, line = res["labels"][j - 1]
            xl.append(f"IMF{j}\n{fault if fault else '-'}\n{res['centers'][j - 1]:.0f} Hz")
        ax_h.set_xticks(range(n_ser))
        ax_h.set_xticklabels(xl, fontsize=8)
        ax_h.set_yticks(range(n_ent))
        ax_h.set_yticklabels([lab for _, lab, _ in ENTROPY_DEFS], fontsize=9)
        ax_h.xaxis.tick_top()
        ax_h.set_title("Entropy of the envelope and of the top VMD modes\n(colour = relative level within a row; "
                       "numbers = actual values)", fontsize=10, pad=44)

        # ---- right column: envelope spectrum with the VMD mode centres, and MSE ----
        gs_r = gs[:, 2].subgridspec(2, 1, hspace=0.55)
        ax_s = fig.add_subplot(gs_r[0])
        f, a = res["spec_f"], res["spec_a"]
        m = f <= 450
        ax_s.plot(f[m], np.maximum(a[m], 1e-9), color="black", lw=0.8)
        ax_s.set_yscale("log")
        lo, hi = log_limits(a[m & (f > 3)], lo_div=4.0, hi_mul=6.0)
        ax_s.set_ylim(lo, hi)
        ax_s.set_xlim(0, 450)
        mark_fault_lines(ax_s, self.F, 450, label_size=6, show_text=False)
        for k in range(min(n_show, K)):
            fault = res["labels"][k][0]
            ax_s.plot(res["centers"][k], hi * 0.55, "v", color=FAULT_COLORS.get(fault, "gray"), ms=6)
            ax_s.text(res["centers"][k], hi * 0.85, f"{k + 1}", ha="center", fontsize=7)
        ax_s.set_xlabel("Frequency (Hz)", fontsize=9)
        ax_s.set_ylabel("Envelope amplitude (g, log)", fontsize=9)
        ax_s.set_title("Envelope spectrum; numbers = IMF centre frequencies", fontsize=9)
        ax_s.grid(True, alpha=0.4, which="both")
        ax_m = fig.add_subplot(gs_r[1])
        scales = np.arange(1, len(res["mse"]) + 1)
        ax_m.plot(scales, res["mse"], "o-", color="darkgreen")
        ax_m.set_xlabel("Scale factor", fontsize=9)
        ax_m.set_ylabel("Sample Entropy", fontsize=9)
        ax_m.set_title(f"Multiscale Sample Entropy (envelope), complexity index {res['mse_ci']:.2f}", fontsize=9)
        ax_m.set_xticks(scales)
        ax_m.grid(True, alpha=0.4)

        fig.suptitle(f"CONDITION: {verdict}", fontsize=22, fontweight="bold",
                     color="tab:red" if primary else (0, 0.55, 0), y=0.985)
        extra = f"   |   also present: {', '.join(others)}" if others else ""
        fig.text(0.5, 0.945, "evidence (dB above noise floor):  " +
                 "   ".join(f"{k} {v:.0f}" for k, v in res["evidence"].items()) + extra +
                 f"   |   window {res['t_start']:.0f}-{res['t_end']:.0f} s   |   analysis #{self.n_results}"
                 f"   |   computed in {res['elapsed']:.1f} s", ha="center", fontsize=10, color="dimgray")
        fig.text(0.5, 0.915, "IMFs are ranked by fault relevance: IMF1 is the strongest fault mode "
                 f"({res['labels'][0][0] or 'none'})", ha="center", fontsize=10, color="navy")
        fig.canvas.draw_idle()

    # ---- real-time driver (wall-clock based) ----
    def run(self, interval_ms=100):
        state = dict(t0=time.perf_counter(), paused=False, pause_t=0.0, paused_total=0.0, done=False)

        def tick(_frame):
            if self.analysis:
                self._poll_analysis()
            if state["done"]:
                if self.analysis and self.pending is None:
                    self._maybe_start_analysis(self.total)
                if self.pending is None:
                    self.ani.event_source.stop()
                return
            if state["paused"]:
                return
            elapsed = (time.perf_counter() - state["t0"] - state["paused_total"]) * self.speed
            n = min(self.total, int(elapsed * self.fs))
            self.render(n)
            if self.analysis:
                self._maybe_start_analysis(n)
            if n >= self.total:
                state["done"] = True
                self.clock.set_text(self.clock.get_text() + "   |   SIMULATION COMPLETE")

        def on_key(event):
            if event.key == " ":
                now = time.perf_counter()
                if state["paused"]:
                    state["paused_total"] += now - state["pause_t"]
                    state["paused"] = False
                else:
                    state["pause_t"], state["paused"] = now, True
            elif event.key in ("r", "R"):
                state.update(t0=time.perf_counter(), paused=False, paused_total=0.0, done=False)
                self.last_idx, self.n_results, self.rows = 0, 0, []
                self.ani.event_source.start()

        self._tick = tick
        self.fig.canvas.mpl_connect("key_press_event", on_key)
        self.ani = FuncAnimation(self.fig, tick, interval=interval_ms, cache_frame_data=False, blit=False)
        plt.show()


# ===================== components figure =====================
def plot_components(comps, fs, F, band=BAND, seconds=0.25, save=None):
    """The four fault signals, then their sum: time domain (left) and envelope spectrum (right)."""
    rows = [("BPFO", "outer race (BPFO)"), ("BPFI", "inner race (BPFI)"), ("BSF", "rolling element (2xBSF)"),
            ("FTF", "cage (FTF)")]
    total = sum(comps.values())
    n0 = int(seconds * fs)
    fig, ax = plt.subplots(len(rows) + 1, 2, figsize=(15, 11))
    for i, (key, title) in enumerate(rows + [("ALL", "COMBINED: baseline + BPFO + BPFI + BSF + FTF")]):
        sig = comps["baseline"] + comps[key] if key != "ALL" else total
        ax[i, 0].plot(np.arange(n0) / fs, sig[:n0], lw=0.5, color=FAULT_COLORS.get(key, "tab:blue"))
        ax[i, 0].set_title(title + " - time domain", fontsize=9)
        ax[i, 0].set_ylabel("g")
        env, fe = demodulate(sig[: int(10 * fs)], fs, band, fe_target=LIVE_ENV_FS, lowpass=None)
        f, a = welch_amplitude_spectrum(env, fe, 3)
        m = f <= 450
        ax[i, 1].plot(f[m], np.maximum(a[m], 1e-9), color="black", lw=0.8)
        ax[i, 1].set_yscale("log")
        ax[i, 1].set_ylim(*log_limits(a[m & (f > 3)], lo_div=4.0, hi_mul=6.0))
        ax[i, 1].set_xlim(0, 450)
        mark_fault_lines(ax[i, 1], F, 450, label_size=7, show_text=(key == "ALL"))
        ax[i, 1].set_title("envelope spectrum", fontsize=9)
    ax[-1, 0].set_xlabel("Time (s)")
    ax[-1, 1].set_xlabel("Frequency (Hz)")
    fig.suptitle("Fault components and their combination", fontweight="bold")
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=130)
    return fig


def parse_weights(text):
    w = {}
    for part in text.split(","):
        k, v = part.split("=")
        w[k.strip().upper()] = float(v)
    return w


def run_bearing_simulation(args):
    """--bearing-sim: simulate a bearing with BPFO + BPFI + BSF + FTF faults and play it back live."""
    F = bearing_frequencies(BEARING["n"], BEARING["d"], BEARING["D"], BEARING["phi"], args.rpm)
    if args.healthy:
        weights = {k: 0.0 for k in ("BPFO", "BPFI", "BSF", "FTF")}
    elif args.weights:
        weights = parse_weights(args.weights)
    else:
        weights = dict(PRESETS[args.dominant])
    print("Bearing: SKF 6205-2RS JEM (9 balls, d = 7.94 mm, D = 39.04 mm, phi = 0 deg)")
    print(f"Shaft {args.rpm:.0f} rpm (fr = {F['fr']:.2f} Hz):  FTF {F['FTF']:.2f}  BSF {F['BSF']:.2f} "
          f"(ball defect at 2xBSF = {2 * F['BSF']:.2f})  BPFO {F['BPFO']:.2f}  BPFI {F['BPFI']:.2f} Hz")
    print("Fault weights:", weights)
    print(f"Simulating {args.duration:g} s at {args.fs_sim} Hz ...")
    t, x, comps = simulate_bearing(args.duration, args.fs_sim, F, weights, seed=args.seed)
    if args.components:
        plot_components(comps, args.fs_sim, F)
    label = args.label if args.label != "Unlabelled" else (
        "healthy" if args.healthy else "weights " + ",".join(f"{k}={v:g}" for k, v in weights.items()))
    sim = BearingSimulator(x, args.fs_sim, F, comps, speed=args.speed, analysis=not args.no_analysis,
                           analysis_seconds=args.analysis_seconds, vmd_k=args.bearing_vmd_k,
                           label=label, log_file=args.log)
    print(f"Starting playback - about {fmt_clock(args.duration / args.speed)} (min:sec). SPACE = pause, R = restart.")
    sim.run(interval_ms=args.interval)




# ===================== DATASET MODE (new: --bearing) =====================
# Plays back YOUR recording (the base PlaybackSimulator, unchanged) and adds a bearing view on top of it:
# the four fault families (BPFO, BPFI, BSF, FTF) are extracted from the real signal, combined in the
# time domain, and their frequencies are marked in the FFT.
DATA_BAND         = (600.0, 900.0)   # resonance band of your rig (impacts ring it in your recordings)
DATA_FTF_HZ       = 9.52             # cage frequency measured from your two Faulty recordings (comb scan)
DATA_BALLS        = 8                # ASSUMED ball count - your data cannot tell it; pass --balls once counted
DATA_FAULT_DB     = 8.0              # threshold calibrated on your two Healthy recordings (10 s windows: max 7.1 dB)
DATA_SPEC_SECONDS = 10.0             # data used for the live envelope spectrum / diagnosis
DATA_COMP_SECONDS = 1.0              # width of the time-domain component display
DATA_LINE_HW_HZ   = 0.8              # half-width of the narrow band-pass around each fault line
DATA_ENV_FS       = 1000.0           # envelope sampling rate


def estimate_shaft_speed(stream, fs, max_seconds=60.0):
    """Shaft speed (Hz): the 15-40 Hz peak whose 2x and 3x are also strong."""
    x = np.asarray(stream[: int(max_seconds * fs)], dtype=float)
    f, a = welch_amplitude_spectrum(x, fs, 5)
    m = (f > 15) & (f < 40)
    cand = f[m][np.argsort(a[m])[::-1][:3]]
    best, score = cand[0], -1.0
    for c in cand:
        s = sum(a[np.argmin(np.abs(f - k * c))] for k in (1, 2, 3))
        if s > score:
            best, score = c, s
    return float(best)


def resolve_bearing_params(args, stream, fs):
    """Fault frequencies for the dataset mode. FTF and BSF follow from the cage frequency (or the geometry);
    BPFO = n x FTF and BPFI = n x (fr - FTF) additionally need the ball count n, which the data cannot give."""
    fr = args.fr if args.fr else estimate_shaft_speed(stream, fs)
    n_known = args.balls is not None
    n = args.balls if n_known else DATA_BALLS
    if args.ball_d and args.pitch_d:
        F = bearing_frequencies(n, args.ball_d, args.pitch_d, args.phi, fr * 60.0)
        how = f"geometry given: d = {args.ball_d} mm, D = {args.pitch_d} mm, phi = {args.phi:g} deg"
    else:
        ftf = args.ftf if args.ftf else DATA_FTF_HZ
        c = 1.0 - 2.0 * ftf / fr                                  # c = (d/D) cos(phi)
        F = {"fr": fr, "FTF": ftf, "BPFO": n * ftf, "BPFI": n * (fr - ftf),
             "BSF": fr * np.cos(np.radians(args.phi)) / (2 * c) * (1 - c ** 2)}
        how = (f"cage frequency FTF = {ftf:g} Hz {'(given)' if args.ftf else '(measured from your Faulty recordings)'}"
               f"  ->  d/D = {c / np.cos(np.radians(args.phi)):.3f} for phi = {args.phi:g} deg")
    known = {"FTF": True, "BSF": True, "BPFO": n_known, "BPFI": n_known}
    return F, known, n, how


def narrowband(x, fs, f0, hw=DATA_LINE_HW_HZ):
    """Zero-phase narrow band-pass around f0 (the part of x that oscillates at f0)."""
    lo, hi = max(f0 - hw, 0.3), f0 + hw
    if hi >= 0.45 * fs:
        return np.zeros_like(x)
    sos = butter(2, [lo / (fs / 2), hi / (fs / 2)], btype="band", output="sos")
    return sosfiltfilt(sos, x)


DATA_HEALTHY_MAX_DB = 7.1            # strongest evidence seen in your two Healthy recordings (10 s windows)
DATA_FAULT_K      = 10               # VMD modes used to separate the COMBINED fault signal
DATA_FAULT_HW_HZ  = 5.0              # VMD mode half-width (Hz)
DATA_ENT_FS       = 500.0            # entropies of the VMD modes are computed at this rate (fast; modes above 200 Hz are flagged)
FAULT_KEYS        = ("BPFO", "BPFI", "BSF", "FTF")
FAULT_WHY = {"FTF": "FTF = cage frequency: damage on the cage / balls modulates the resonance once per cage turn,",
             "BSF": "BSF = ball spin: a ball defect hits both races twice per spin (lines at 2xBSF +/- FTF),",
             "BPFO": "BPFO = outer-race pass frequency: an outer-race defect is hit once per ball passage,",
             "BPFI": "BPFI = inner-race pass frequency: lines at BPFI +/- shaft speed,"}


def fault_components(env, fe, F):
    """The four fault components of an envelope signal (narrow band-pass at each fault's lines) and their SUM."""
    comps = {fault: sum(narrowband(env, fe, fl) for _, fl in lines) for fault, lines in fault_lines(F).items()}
    return comps, sum(comps.values())


def analyze_reference(ref_data, fs, F, known, band):
    """Healthy reference: strongest evidence of any (known) fault, and the median envelope spectrum."""
    stream = ref_data.T.reshape(-1)
    L = int(DATA_SPEC_SECONDS * fs)
    mx, specs, fq = 0.0, [], None
    for k in range(len(stream) // L):
        env, fe = demodulate(stream[k * L:(k + 1) * L], fs, band, fe_target=DATA_ENV_FS, lowpass=None)
        fq, ae = welch_amplitude_spectrum(env, fe, 3)
        specs.append(ae)
        ev = diagnose(fq, ae, F)
        mx = max(mx, max(v for kk, v in ev.items() if known[kk]))
    return mx, fq, np.median(np.array(specs), axis=0)


def label_by_correlation(C, rms):
    """Which fault does each VMD mode belong to? Largest |correlation| with the four known components;
    two close competitors -> 'A+B' (mixed); weak -> 'mixed'; tiny energy -> 'residual'."""
    labels = []
    for i in range(C.shape[0]):
        a = np.abs(C[i])
        o = np.argsort(a)[::-1]
        if rms[i] < 0.05 * rms.max():
            labels.append("residual")
        elif a[o[0]] < 0.25:
            labels.append("mixed")
        elif a[o[1]] >= 0.35 * a[o[0]]:
            labels.append(f"{FAULT_KEYS[o[0]]}+{FAULT_KEYS[o[1]]}")
        else:
            labels.append(FAULT_KEYS[o[0]])
    return labels


def analyze_combined_fault(window, fs, F, known, band, K, halfwidth_hz, t_start):
    """COMBINED FAULT SIGNAL -> VMD -> IMFs aligned with fault types -> entropy features per IMF.
    1) envelope of the resonance band of the raw window   2) the four fault components and their sum
    3) VMD of the SUM into K modes   4) each mode matched to a fault (correlation with the components)
    5) the seven entropies for the sum and for every mode (the fault feature vector)."""
    t0 = time.perf_counter()
    env, fe = demodulate(window, fs, band, fe_target=DATA_ENV_FS, lowpass=None)
    comps, combined = fault_components(env, fe, F)
    fq, ae = welch_amplitude_spectrum(env, fe, 3)
    evidence = diagnose(fq, ae, F)

    modes, om, n_iter = vmd_decompose(combined, K, alpha=(fe / halfwidth_hz) ** 2)
    L = modes.shape[1]
    comps = {k: v[:L] for k, v in comps.items()}
    combined, env = combined[:L], env[:L]
    centers = om * fe
    rms = np.sqrt(np.mean(modes ** 2, axis=1))
    C = np.array([[np.corrcoef(modes[i], comps[k])[0, 1] if rms[i] > 1e-12 else 0.0 for k in FAULT_KEYS]
                  for i in range(K)])
    labels = label_by_correlation(C, rms)
    order = np.argsort(rms)[::-1]                              # IMF1 = strongest mode
    modes, centers, rms, C = modes[order], centers[order], rms[order], C[order]
    labels = [labels[i] for i in order]
    recon_err = float(np.linalg.norm(modes.sum(axis=0) - combined) / np.linalg.norm(combined))

    q = max(1, int(round(fe / DATA_ENT_FS)))
    fe_ent = fe / q
    names = ["Combined"] + [f"IMF{k + 1}" for k in range(K)]
    series = [combined] + list(modes)
    matrix = np.full((len(ENTROPY_DEFS), len(series)), np.nan)
    for j, s in enumerate(series):
        s = decimate(s, q) if q > 1 else s
        vals = compute_entropies_bearing(s, fe_ent)
        for i, (key, _, _) in enumerate(ENTROPY_DEFS):
            matrix[i, j] = vals[key]
    in_band = np.ones(len(series), dtype=bool)
    if q > 1:
        in_band[1:] = centers <= 0.4 * fe_ent
    mse = multiscale_sample_entropy(decimate(combined, q) if q > 1 else combined, MSE_MAX_SCALE, EN_M)
    return dict(t_start=t_start, t_end=t_start + len(window) / fs, fs=fs, fe=fe, fe_ent=fe_ent, K=K, F=F,
                known=known, env=env, comps=comps, combined=combined, modes=modes, centers=centers, rms=rms,
                corr=C, labels=labels, recon_err=recon_err, names=names, matrix=matrix, in_band=in_band,
                mse=mse, mse_ci=float(np.nansum(mse)), evidence=evidence, n_iter=n_iter,
                elapsed=time.perf_counter() - t0)


def fault_vmd_result_to_row(res, label, verdict):
    row = {"Timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "Label": label,
           "WindowStart_s": round(res["t_start"], 3), "WindowEnd_s": round(res["t_end"], 3),
           "Condition": verdict, "AcqFs": res["fs"], "EnvFs": res["fe"], "EntropyFs": res["fe_ent"], "K": res["K"],
           "ReconstructionError": res["recon_err"]}
    for k, v in res["evidence"].items():
        row[f"Evidence_{k}_dB"] = round(v, 2)
    for k in range(res["K"]):
        row[f"IMF{k + 1}_Fault"] = res["labels"][k]
        row[f"IMF{k + 1}_CenterFreqHz"] = res["centers"][k]
        row[f"IMF{k + 1}_RMS"] = res["rms"][k]
        for j, fk in enumerate(FAULT_KEYS):
            row[f"IMF{k + 1}_corr_{fk}"] = res["corr"][k, j]
    for i, (key, _, _) in enumerate(ENTROPY_DEFS):
        for j, name in enumerate(res["names"]):
            row[f"{key}_{name}"] = res["matrix"][i, j]
    for s, v in enumerate(res["mse"], start=1):
        row[f"MSE_scale{s}"] = v
    row["MSE_ComplexityIndex"] = res["mse_ci"]
    return row


class BearingDatasetSimulator(PlaybackSimulator):
    """The base PlaybackSimulator (your dataset replayed like a recording: time domain, FFT with Mode labels)
    plus the bearing view:
      * four fault components extracted from the real signal and their COMBINED sum in the time domain
      * the four fault frequencies marked on your FFT and on an envelope FFT that also shows the noise floor,
        the detection threshold, (optionally) the Healthy reference and the measured peaks of the winning fault
      * a decision panel: what defines the fault (evidence per fault vs the threshold) and why
      * a second window: COMBINED fault signal -> VMD -> IMFs aligned with the fault types -> entropy features"""

    def __init__(self, data, fs, F, known, n_assumed, fault_db=DATA_FAULT_DB, healthy_max_db=DATA_HEALTHY_MAX_DB,
                 band=DATA_BAND, ref_spec=None, fault_k=DATA_FAULT_K, fault_analysis=True, fault_log=None, **kw):
        self._brg_ready = False
        super().__init__(data, fs, **kw)
        self.F, self.known, self.n_assumed, self.fault_db, self.band = F, known, n_assumed, fault_db, band
        self.healthy_max_db, self.ref_spec, self.fault_k, self.fault_log = healthy_max_db, ref_spec, fault_k, fault_log
        self.keys = list(FAULT_KEYS)
        self.q = max(1, int(round(fs / DATA_ENV_FS)))
        self.fe = fs / self.q
        self.spec_sw = int(DATA_SPEC_SECONDS * fs)
        self.spec_step = max(1, int(round(0.5 * fs)))
        self.last_spec_idx = -1
        self.cw = int(round(DATA_COMP_SECONDS * self.fe))
        self.env_lim = (1e-5, 1e-2)
        self.lines_c = {}
        self._cur_n, self._start_n = 0, 0
        self.verdict, self.primary = "ACQUIRING", None
        self.fault_analysis = bool(fault_analysis) and self.win <= self.total
        self.executor2 = ThreadPoolExecutor(max_workers=1) if self.fault_analysis else None
        self.pending2, self.last_idx2, self.rows2, self.n_results2, self.fig3, self._timer = None, 0, [], 0, None, None

        # ---- re-arrange the base figure: your two axes stay, four new ones are added ----
        self.fig.set_size_inches(15, 12.8)
        gs = self.fig.add_gridspec(3, 2, left=0.065, right=0.985, top=0.845, bottom=0.145, hspace=0.62,
                                   wspace=0.17, height_ratios=[1.0, 1.15, 1.0])
        self.ax_t.set_subplotspec(gs[0, 0])
        self.ax_f.set_subplotspec(gs[0, 1])
        self.ax_c = self.fig.add_subplot(gs[1, 0])
        self.ax_s = self.fig.add_subplot(gs[2, 0])
        self.ax_e = self.fig.add_subplot(gs[1, 1])
        self.ax_d = self.fig.add_subplot(gs[2, 1])

        # four fault components (envelope domain), one trace each, every trace scaled to its own peak
        t_c = np.arange(self.cw) / self.fe
        for i, k in enumerate(self.keys):
            self.lines_c[k], = self.ax_c.plot(t_c, np.full(self.cw, 3.0 - i), color=FAULT_COLORS[k], lw=0.9)
        self.ax_c.set_xlim(0, t_c[-1])
        self.ax_c.set_ylim(-0.6, 3.6)
        self.ax_c.set_yticks([3, 2, 1, 0])
        self.ax_c.set_yticklabels(self.keys, fontsize=8)
        self.ax_c.set_xlabel("Time (s) - last 1 s")
        self.ax_c.set_title("TIME DOMAIN - the four fault components extracted from YOUR signal\n"
                            "(envelope domain, each drawn to its own peak)", fontsize=10)
        self.ax_c.grid(True, alpha=0.4)

        # their sum vs the measured envelope
        self.line_env, = self.ax_s.plot(t_c, np.zeros(self.cw), color="lightgray", lw=1.4, label="measured envelope")
        self.line_sum, = self.ax_s.plot(t_c, np.zeros(self.cw), color="black", lw=1.0,
                                        label="COMBINED = BPFO + BPFI + BSF + FTF")
        self.ax_s.set_xlim(0, t_c[-1])
        self.ax_s.set_xlabel("Time (s) - last 1 s")
        self.ax_s.set_ylabel("Envelope amplitude (g)")
        self.ax_s.set_title("TIME DOMAIN - the four components COMBINED (this sum is what VMD separates again)",
                            fontsize=10, fontweight="bold")
        self.ax_s.legend(loc="upper right", fontsize=8)
        self.ax_s.grid(True, alpha=0.4)

        # envelope FFT: spectrum, noise floor, threshold, optional Healthy reference, fault markers
        self.line_e, = self.ax_e.plot([0, 450], [1e-9, 1e-9], color="black", lw=0.9, label="envelope spectrum (this recording)")
        self.line_ref, = self.ax_e.plot([0, 450], [1e-9, 1e-9], color="tab:blue", lw=1.1, alpha=0.9,
                                        label="Healthy reference (median)" if ref_spec is not None else "_nolegend_")
        self.line_nf, = self.ax_e.plot([0, 450], [1e-9, 1e-9], color="gray", ls="--", lw=1.1, label="noise floor")
        self.line_thr, = self.ax_e.plot([0, 450], [1e-9, 1e-9], color="darkred", ls=":", lw=1.6,
                                        label=f"fault threshold (noise floor + {fault_db:g} dB)")
        self.persist = {self.line_e, self.line_ref, self.line_nf, self.line_thr}
        self.ax_e.set_yscale("log")
        self.ax_e.set_xlim(0, 450)
        self.ax_e.set_ylim(*self.env_lim)
        self.note = {} if all(known[k] for k in ("BPFO", "BPFI")) else {
            "BPFO": f" (n={n_assumed} assumed)", "BPFI": f" (n={n_assumed} assumed)"}
        mark_fault_lines(self.ax_e, F, 450, note=self.note)
        self.ax_e.set_xlabel("Frequency (Hz)")
        self.ax_e.set_ylabel("Envelope amplitude (g, log scale)")
        self.ax_e.set_title(f"FFT OF THE ENVELOPE ({band[0]:.0f}-{band[1]:.0f} Hz band): a fault = its lines rise above the threshold",
                            fontsize=10, fontweight="bold")
        self.ax_e.legend(loc="upper right", fontsize=7)
        self.ax_e.grid(True, alpha=0.4, which="both")

        # the same four frequencies on your original FFT axis
        tr = self.ax_f.get_xaxis_transform()
        for key, fq, name in [("BPFO", F["BPFO"], "BPFO" + ("" if known["BPFO"] else "*")),
                              ("BPFI", F["BPFI"], "BPFI" + ("" if known["BPFI"] else "*")),
                              ("BSF", 2 * F["BSF"], "2xBSF"), ("FTF", F["FTF"], "FTF")]:
            self.ax_f.axvline(fq, color=FAULT_COLORS[key], ls=":", lw=1.3)
            self.ax_f.text(fq, 0.02, f" {name} {fq:.1f}", transform=tr, color=FAULT_COLORS[key], rotation=90,
                           va="bottom", ha="right", fontsize=8)
        self.ax_f.set_title("FFT Spectrum (3 Harmonic Modes) - dotted lines = BPFO, BPFI, 2xBSF, FTF", fontsize=10)

        # decision panel: what DEFINES the fault
        self.bars = self.ax_d.barh(range(4), [0, 0, 0, 0], color=[FAULT_COLORS[k] for k in self.keys],
                                   edgecolor="black", linewidth=0.6)
        self.bar_txt = [self.ax_d.text(0, i, "", va="center", fontsize=9, fontweight="bold") for i in range(4)]
        self.ax_d.axvspan(0, healthy_max_db, color="green", alpha=0.13,
                          label=f"Healthy range (your Healthy recordings: up to {healthy_max_db:.1f} dB)")
        self.ax_d.axvline(fault_db, color="red", ls="--", lw=1.8, label=f"fault threshold {fault_db:g} dB")
        self.ax_d.set_yticks(range(4))
        self.ax_d.set_yticklabels([k + ("" if known[k] else "*\n(n unknown)") for k in self.keys], fontsize=9)
        self.ax_d.invert_yaxis()
        self.ax_d.set_xlim(0, 30)
        self.ax_d.set_xlabel("evidence = how far the fault's strongest lines stand above the noise floor (dB)")
        self.ax_d.set_title("WHAT DEFINES THE FAULT: evidence per fault family vs the threshold", fontsize=10, fontweight="bold")
        self.ax_d.legend(loc="lower right", fontsize=7)
        self.ax_d.grid(True, alpha=0.4, axis="x")

        self.banner = self.fig.text(0.5, 0.935, "CONDITION: ACQUIRING...", ha="center", fontsize=22,
                                    fontweight="bold", color="gray")
        self.sub = self.fig.text(0.5, 0.897, "", ha="center", fontsize=10, color="dimgray")
        self.why = self.fig.text(0.515, 0.115, "", ha="left", va="top", fontsize=9, color="black",
                                 bbox=dict(boxstyle="round", facecolor="#fff8e1", edgecolor="gray"))
        self.dyn = []
        self._brg_ready = True
        self._render_bearing(0)

        if self.fault_analysis:
            self.fig3 = plt.figure(figsize=(17.5, 10.5))
            self.fig3.text(0.5, 0.5, f"Waiting for the first {self.win / fs:g} s of data...\n"
                           "(COMBINED fault signal -> VMD -> IMFs aligned with the faults -> entropy features)",
                           ha="center", va="center", fontsize=14, color="gray")
            self.info.set_text(f"Combined-fault VMD + entropy: first result after {self.win / fs:g} s of data...")
        for fg in (self.fig, self.fig3):
            if fg is not None:
                fg.canvas.mpl_connect("close_event", self._on_close)

    # ------------------------------------------------------------------ housekeeping
    def _on_close(self, event):
        super()._on_close(event)
        if self.executor2 is not None:
            self.executor2.shutdown(wait=False, cancel_futures=True)
        if self._timer is not None:
            self._timer.stop()

    def run(self, start_group=1, interval_ms=100):
        self._start_n = (max(1, start_group) - 1) * self.group_len
        if self.fault_analysis:
            self._timer = self.fig.canvas.new_timer(interval=300)      # polls the VMD analysis, also after playback ends
            self._timer.add_callback(self._brg_timer)
            self._timer.start()
        super().run(start_group=start_group, interval_ms=interval_ms)

    def _brg_timer(self):
        self._poll_fault_analysis()
        self._maybe_start_fault_analysis(self._cur_n)

    def render(self, n):
        super().render(n)                                   # your original time domain + FFT, untouched
        self._cur_n = int(n)
        if getattr(self, "_brg_ready", False):
            self._render_bearing(n)

    # ------------------------------------------------------------------ live bearing view
    def _bearing_idle(self):
        flat = np.zeros(self.cw)
        for i, k in enumerate(self.keys):
            self.lines_c[k].set_ydata(np.full(self.cw, 3.0 - i))
        self.line_env.set_ydata(flat)
        self.line_sum.set_ydata(flat)
        self.line_e.set_data([0, 450], [1e-9, 1e-9])
        for ln in (self.line_ref, self.line_nf, self.line_thr):
            ln.set_data([0, 450], [1e-9, 1e-9])
        for b_, t_ in zip(self.bars, self.bar_txt):
            b_.set_width(0)
            t_.set_text("")
        self.banner.set_text("CONDITION: ACQUIRING...")
        self.banner.set_color("gray")
        self.sub.set_text("")
        self.why.set_text("")

    def _render_bearing(self, n):
        n = int(min(max(n, 0), self.total))
        idx = n // self.spec_step
        if idx == self.last_spec_idx and n < self.total:
            return
        self.last_spec_idx = idx
        have = min(n, self.spec_sw)
        if have < 3 * self.fs:                              # need >= 3 s of data
            self._bearing_idle()
            return
        seg = self.stream[n - have:n]
        env, fe = demodulate(seg, self.fs, self.band, fe_target=DATA_ENV_FS, lowpass=None)

        # ---- envelope FFT: spectrum, noise floor, threshold, Healthy reference ----
        fq, ae = welch_amplitude_spectrum(env, fe, 3)
        m = fq <= 450
        df = fq[1] - fq[0]
        nf = median_filter(ae, size=max(3, int(30.0 / df)), mode="nearest")      # running median = local noise floor
        thr = nf * 10 ** (self.fault_db / 20.0)
        self.line_e.set_data(fq[m], np.maximum(ae[m], 1e-12))
        self.line_nf.set_data(fq[m], np.maximum(nf[m], 1e-12))
        self.line_thr.set_data(fq[m], np.maximum(thr[m], 1e-12))
        if self.ref_spec is not None:
            self.line_ref.set_data(fq[m], np.maximum(np.interp(fq[m], self.ref_spec[0], self.ref_spec[1]), 1e-12))
        lo, hi = log_limits(ae[m & (fq > 3)], lo_div=4.0, hi_mul=3.0)
        self.env_lim = (smooth_log(self.env_lim[0], lo), smooth_log(self.env_lim[1], hi))
        self.ax_e.set_ylim(*self.env_lim)
        for t in list(self.ax_e.texts):
            t.remove()
        for ln in list(self.ax_e.lines):
            if ln not in self.persist:
                ln.remove()
        mark_fault_lines(self.ax_e, self.F, 450, note=self.note)

        # ---- diagnosis ----
        ev = diagnose(fq, ae, self.F)
        used = {k: v for k, v in ev.items() if self.known[k]}
        top, topv = max(used.items(), key=lambda kv: kv[1])
        faulty = topv >= self.fault_db
        self.primary = top if faulty else None
        self.verdict = f"{top} FAULTY" if faulty else "HEALTHY"
        others = [k for k, v in used.items() if k != top and v >= self.fault_db] if faulty else []
        self.banner.set_text(f"CONDITION: {self.verdict}")
        self.banner.set_color("tab:red" if faulty else (0, 0.55, 0))
        self.sub.set_text(f"fault threshold {self.fault_db:g} dB above the noise floor   |   "
                          + "   ".join(f"{k} {ev[k]:.0f}{'' if self.known[k] else '*'}" for k in self.keys)
                          + (f"   |   also present: {', '.join(others)}" if others else ""))

        # ---- measured peaks of the winning fault (or of the strongest family when healthy), with their dB ----
        show = top
        lines = fault_lines(self.F)[show]
        meas = []
        for name, fl in lines:
            w = (fq > fl - 0.6) & (fq < fl + 0.6)
            if not w.any():
                continue
            j = np.argmax(np.where(w, ae, 0))
            snr = line_snr_db(fq, ae, fl)
            meas.append((name, fl, snr))
            strong = snr >= self.fault_db
            self.ax_e.plot(fq[j], ae[j], "o", ms=8, markerfacecolor=FAULT_COLORS[show] if strong else "white",
                           markeredgecolor=FAULT_COLORS[show], markeredgewidth=1.6, zorder=5)
            self.ax_e.annotate(f"{name}\n{snr:+.0f} dB", (fq[j], ae[j]), textcoords="offset points", xytext=(6, 10),
                               fontsize=8, fontweight="bold" if strong else "normal", color=FAULT_COLORS[show],
                               arrowprops=dict(arrowstyle="-", color=FAULT_COLORS[show], lw=0.7))

        # ---- decision panel ----
        xmax = max(30.0, 1.25 * max(ev.values()))
        self.ax_d.set_xlim(0, xmax)
        for i, k in enumerate(self.keys):
            v = ev[k]
            self.bars[i].set_width(v)
            self.bars[i].set_alpha(1.0 if self.known[k] else 0.35)
            self.bars[i].set_hatch("" if self.known[k] else "//")
            self.bar_txt[i].set_position((v + 0.4, i))
            self.bar_txt[i].set_text(f"{v:.1f} dB" + ("  FAULT" if (self.known[k] and v >= self.fault_db) else ""))
        top2 = sorted(meas, key=lambda t: -t[2])[:2]
        if faulty:
            msg = (f"WHY {top}: its lines stand above the noise floor - "
                   + " and ".join(f"{nm} ({fl:.1f} Hz) {s:+.0f} dB" for nm, fl, s in top2) + ".\n"
                   f"Evidence {topv:.1f} dB >= threshold {self.fault_db:g} dB; your Healthy recordings never exceed "
                   f"{self.healthy_max_db:.1f} dB.\n{FAULT_WHY[top]} so peaks appear at its harmonics.")
        else:
            msg = (f"NO FAULT: the strongest evidence ({top} {topv:.1f} dB) stays below the {self.fault_db:g} dB threshold,\n"
                   f"inside the Healthy range (up to {self.healthy_max_db:.1f} dB): no line stands out of the noise floor.\n"
                   "A fault is declared when a family's lines rise clearly above the noise floor.")
        self.why.set_text(msg)

        # ---- the four fault components and their sum (time domain) ----
        comps, total = fault_components(env, fe, self.F)
        cw = min(self.cw, len(env))
        for i, k in enumerate(self.keys):
            tail = comps[k][-cw:]
            peak = float(np.max(np.abs(comps[k]))) or 1.0
            y = np.full(self.cw, 3.0 - i)
            y[-cw:] = (3.0 - i) + 0.42 * tail / peak
            self.lines_c[k].set_ydata(y)
        self.ax_c.set_yticklabels([f"{k}{'' if self.known[k] else '*'}\npeak {float(np.max(np.abs(comps[k]))) * 1000:.1f} mg"
                                   for k in self.keys], fontsize=8)
        e_tail, s_tail = np.zeros(self.cw), np.zeros(self.cw)
        e_tail[-cw:], s_tail[-cw:] = env[-cw:], total[-cw:]
        self.line_env.set_ydata(e_tail)
        self.line_sum.set_ydata(s_tail)
        lim = 1.15 * max(float(np.max(np.abs(e_tail))), float(np.max(np.abs(s_tail))), 1e-6)
        self.ax_s.set_ylim(-lim, lim)
        if not all(self.known[k] for k in ("BPFO", "BPFI")):
            self.ax_c.set_xlabel("Time (s) - last 1 s        * position depends on the ball count (assumed n="
                                 f"{self.n_assumed}, not used for the verdict)")

    # ------------------------------------------------------------------ combined fault -> VMD -> entropy window
    def _maybe_start_fault_analysis(self, n):
        if not self.fault_analysis:
            return
        idx = (n - self._start_n) // self.win
        if idx < self.last_idx2:                           # playback was restarted
            self.last_idx2, self.n_results2, self.rows2 = 0, 0, []
        if idx > self.last_idx2 and self.pending2 is None and n >= self.win:
            self.last_idx2 = idx
            window = self.stream[n - self.win:n].copy()
            t_start = (n - self.win) / self.fs
            self.pending2 = self.executor2.submit(analyze_combined_fault, window, self.fs, self.F, self.known,
                                                  self.band, self.fault_k, DATA_FAULT_HW_HZ, t_start)
            self.info.set_text(f"Combined-fault VMD + entropy: analysing {t_start:.0f}-{n / self.fs:.0f} s ...")

    def _poll_fault_analysis(self):
        if self.pending2 is not None and self.pending2.done():
            fut, self.pending2 = self.pending2, None
            try:
                res = fut.result()
            except Exception as exc:
                print(f"Combined-fault analysis failed: {exc}")
                self.info.set_text(f"Combined-fault VMD + entropy failed: {exc}")
                return
            self._show_fault_analysis(res)

    def _show_fault_analysis(self, res):
        self.n_results2 += 1
        ev = res["evidence"]
        used = {k: v for k, v in ev.items() if self.known[k]}
        top, topv = max(used.items(), key=lambda kv: kv[1])
        verdict = f"{top} FAULTY" if topv >= self.fault_db else "HEALTHY"
        lab1 = res["labels"][0]
        print(f"[combined-fault analysis #{self.n_results2}] {res['t_start']:.0f}-{res['t_end']:.0f} s | {verdict} | "
              f"IMF1 = {lab1} ({res['centers'][0]:.1f} Hz) | reconstruction error {res['recon_err']:.2f} | "
              f"{res['elapsed']:.1f} s")
        self.info.set_text(f"Combined-fault VMD + entropy #{self.n_results2} ({res['t_start']:.0f}-{res['t_end']:.0f} s):  "
                           f"{verdict}   |   IMF1 = {lab1} {res['centers'][0]:.1f} Hz   |   entropies at {res['fe_ent']:g} Hz")
        if self.fault_log:
            self.rows2.append(fault_vmd_result_to_row(res, self.label, verdict))
            try:
                pd.DataFrame(self.rows2).to_excel(self.fault_log, index=False)
            except Exception as exc:
                print(f"Could not write {self.fault_log}: {exc}")
        self._draw_fault_vmd_figure(res, verdict, topv)

    def _draw_fault_vmd_figure(self, res, verdict, topv):
        fig = self.fig3
        fig.clf()
        K, fe, known = res["K"], res["fe"], res["known"]
        gs = fig.add_gridspec(K, 3, width_ratios=[1.0, 1.15, 1.7], hspace=0.95, wspace=0.28,
                              left=0.045, right=0.985, top=0.80, bottom=0.07)
        zoom = int(1.0 * fe)
        t = np.arange(zoom) / fe

        def col(lab):
            return FAULT_COLORS.get(lab.split("+")[0], "gray") if "+" not in lab and lab in FAULT_COLORS else "gray"

        # ---- left: the K IMFs of the COMBINED fault signal ----
        for k in range(K):
            ax = fig.add_subplot(gs[k, 0])
            lab = res["labels"][k]
            ax.plot(t, res["modes"][k][:zoom], color=col(lab), lw=0.9)
            star = "*" if (lab in FAULT_KEYS and not known[lab]) else ""
            ax.set_title(f"IMF{k + 1} = {lab}{star}   {res['centers'][k]:.1f} Hz   RMS {res['rms'][k] * 1000:.2f} mg",
                         fontsize=8, pad=2, color=col(lab), fontweight="bold" if k == 0 else "normal")
            ax.set_xlim(0, t[-1])
            ax.tick_params(labelsize=7)
            ax.grid(True, alpha=0.3)
            if k < K - 1:
                ax.set_xticklabels([])
        ax.set_xlabel("Time (s) - first 1 s of each mode", fontsize=8)

        # ---- middle top: spectrum of the combined signal, IMF centre frequencies, fault lines ----
        half = K // 2
        ax_s = fig.add_subplot(gs[:half, 1])
        f, a = welch_amplitude_spectrum(res["combined"], fe, 3)
        m = f <= 300
        ax_s.plot(f[m], np.maximum(a[m], 1e-12), color="black", lw=0.9)
        ax_s.set_yscale("log")
        lo, hi = log_limits(a[m & (f > 3)], lo_div=50.0, hi_mul=8.0)
        ax_s.set_ylim(lo, hi)
        ax_s.set_xlim(0, 300)
        mark_fault_lines(ax_s, res["F"], 300, label_size=6, show_text=False)
        for k in range(K):
            ax_s.plot(res["centers"][k], hi * 0.45, "v", color=col(res["labels"][k]), ms=6)
            ax_s.text(res["centers"][k], hi * 0.8, f"{k + 1}", ha="center", fontsize=7)
        ax_s.set_xlabel("Frequency (Hz)", fontsize=8)
        ax_s.set_ylabel("amplitude (g, log)", fontsize=8)
        ax_s.set_title("Combined fault signal: spectrum, fault lines and the IMF centres (numbers)", fontsize=9)
        ax_s.grid(True, alpha=0.4, which="both")

        # ---- middle bottom: alignment of every IMF with the four known fault components ----
        ax_c = fig.add_subplot(gs[half:, 1])
        Cabs = np.abs(res["corr"])
        ax_c.imshow(Cabs, aspect="auto", cmap="Greens", vmin=0, vmax=1)
        for i in range(K):
            for j in range(4):
                ax_c.text(j, i, f"{Cabs[i, j]:.2f}", ha="center", va="center", fontsize=7)
            lab = res["labels"][i]
            if lab in FAULT_KEYS:
                j = FAULT_KEYS.index(lab)
                ax_c.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="red", linewidth=2))
        ax_c.set_xticks(range(4))
        ax_c.set_xticklabels([k + ("" if known[k] else "*") for k in FAULT_KEYS], fontsize=8)
        ax_c.set_yticks(range(K))
        ax_c.set_yticklabels([f"IMF{i + 1}  {res['centers'][i]:.0f} Hz" for i in range(K)], fontsize=7)
        ax_c.xaxis.tick_top()
        ax_c.set_title("ALIGNMENT: |correlation| of each IMF with the four\ncombined faults (red box = assigned fault)", fontsize=9, pad=24)

        # ---- right: entropy features of the combined signal and of every IMF ----
        ax_h = fig.add_subplot(gs[:, 2])
        M, names = res["matrix"], res["names"]
        n_ent, n_ser = M.shape
        norm_m = np.full_like(M, np.nan)
        for i in range(n_ent):
            row = M[i]
            ok = np.isfinite(row) & res["in_band"]
            if ok.sum():
                lo_, hi_ = row[ok].min(), row[ok].max()
                norm_m[i] = (row - lo_) / (hi_ - lo_) if hi_ > lo_ else 0.5
        ax_h.imshow(np.clip(np.nan_to_num(norm_m, nan=0.0), 0, 1), aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
        for i in range(n_ent):
            for j in range(n_ser):
                v = M[i, j]
                ax_h.text(j, i, "n/a" if np.isnan(v) else f"{v:.2f}", ha="center", va="center", fontsize=7,
                          color="black" if res["in_band"][j] else "dimgray")
        for j in range(n_ser):
            if not res["in_band"][j]:
                ax_h.add_patch(Rectangle((j - 0.5, -0.5), 1, n_ent, fill=True, facecolor="white", alpha=0.55,
                                         hatch="///", edgecolor="gray", linewidth=0))
        xl = ["Combined"] + [f"IMF{j}\n{res['labels'][j - 1].replace('+', '+' + chr(10))}\n{res['centers'][j - 1]:.0f}Hz"
                             for j in range(1, n_ser)]
        ax_h.set_xticks(range(n_ser))
        ax_h.set_xticklabels(xl, fontsize=6.5)
        ax_h.set_yticks(range(n_ent))
        ax_h.set_yticklabels([lab for _, lab, _ in ENTROPY_DEFS], fontsize=9)
        ax_h.xaxis.tick_top()
        ax_h.set_title("ENTROPY FEATURES per IMF (the fault feature vector)\n"
                       "colour = relative level in a row; numbers = values", fontsize=9, pad=40)

        col_t = "tab:red" if "FAULTY" in verdict else (0, 0.55, 0)
        fig.suptitle("COMBINED FAULT SIGNAL  ->  VMD  ->  IMFs ALIGNED WITH THE FAULTS  ->  ENTROPY FEATURES",
                     fontsize=15, fontweight="bold", y=0.985)
        fig.text(0.5, 0.925, f"CONDITION: {verdict}", ha="center", fontsize=18, fontweight="bold", color=col_t)
        fig.text(0.5, 0.885, f"window {res['t_start']:.0f}-{res['t_end']:.0f} s   |   analysis #{self.n_results2}   |   "
                 f"VMD K={K}, half-width {DATA_FAULT_HW_HZ:g} Hz, reconstruction error {res['recon_err']:.2f}   |   "
                 f"entropies computed at {res['fe_ent']:g} Hz (hatched = mode above {0.4 * res['fe_ent']:.0f} Hz, ignore)   |   "
                 f"computed in {res['elapsed']:.1f} s", ha="center", fontsize=9, color="dimgray")
        fig.text(0.5, 0.855, "IMF1 = strongest mode. A mode is 'A+B' (mixed) when two fault lines lie too close for VMD to split them; "
                 "* = position depends on the assumed ball count.", ha="center", fontsize=9, color="navy")
        fig.canvas.draw_idle()


def make_bearing_dataset_simulator(args, data, fs):
    stream = data.T.reshape(-1)
    F, known, n, how = resolve_bearing_params(args, stream, fs)
    band = tuple(args.band)
    fault_db, healthy_max, ref_spec = DATA_FAULT_DB, DATA_HEALTHY_MAX_DB, None
    print(f"\nBearing view for your dataset  (resonance band {band[0]:.0f}-{band[1]:.0f} Hz)")
    print(f"  shaft speed fr = {F['fr']:.2f} Hz ({F['fr'] * 60:.0f} rpm)")
    print(f"  {how}")
    print(f"  FTF  {F['FTF']:.2f} Hz   BSF {F['BSF']:.2f} Hz (ball defect at 2xBSF = {2 * F['BSF']:.2f} Hz)")
    if known["BPFO"]:
        print(f"  n = {n} balls:  BPFO {F['BPFO']:.2f} Hz   BPFI {F['BPFI']:.2f} Hz")
    else:
        print(f"  BPFO {F['BPFO']:.1f} Hz and BPFI {F['BPFI']:.1f} Hz are drawn for an ASSUMED n = {n} (marked *) and are NOT used "
              f"for the verdict - count the balls and pass --balls N.")
    if args.reference:
        print(f"  Healthy reference {args.reference}: calibrating the threshold and drawing its spectrum ...")
        _, fs_ref, _, ref = load_dataset(args.reference)
        healthy_max, rf, ra = analyze_reference(ref, fs_ref, F, known, band)
        ref_spec = (rf, ra)
        fault_db = healthy_max + 1.0
    if args.fault_db is not None:
        fault_db = args.fault_db
    print(f"  fault threshold = {fault_db:.1f} dB above the noise floor "
          f"({'calibrated on the Healthy reference' if args.reference and args.fault_db is None else 'given' if args.fault_db is not None else 'calibrated on your two Healthy recordings'})")
    keep_base = bool(args.keep_base_analysis) and not args.no_analysis
    fault_log = args.log
    if args.log and keep_base:
        fault_log = args.log.replace(".xlsx", "") + "_fault_vmd.xlsx"
    print(f"  analysis windows: combined-fault VMD + entropy{' + your original raw-signal VMD' if keep_base else ''}\n")
    return BearingDatasetSimulator(data, fs, F, known, n, fault_db=fault_db, healthy_max_db=healthy_max, band=band,
                                   ref_spec=ref_spec, fault_k=args.fault_k, fault_analysis=not args.no_analysis,
                                   fault_log=fault_log, xmax=args.xmax if args.xmax else 300.0, speed=args.speed,
                                   analysis=keep_base, analysis_seconds=args.analysis_seconds, use_vmd=not args.no_vmd,
                                   vmd_k=args.vmd_k, analysis_fs=args.analysis_fs, label=args.label,
                                   log_file=args.log if keep_base else None)



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
    # ---- bearing simulation (new) ----
    ap.add_argument("--bearing-sim", action="store_true",
                    help="simulate a bearing signal containing BPFO + BPFI + BSF + FTF instead of loading an Excel file")
    ap.add_argument("--dominant", choices=["BPFO", "BPFI", "BSF", "FTF"], default="BPFO",
                    help="(--bearing-sim) scenario preset in which this fault dominates (default BPFO)")
    ap.add_argument("--weights", default=None, help="(--bearing-sim) explicit weights, e.g. BPFO=1,BPFI=0.3,BSF=0.5,FTF=1.5")
    ap.add_argument("--healthy", action="store_true", help="(--bearing-sim) no faults, baseline only")
    ap.add_argument("--duration", type=float, default=DURATION_S, help="(--bearing-sim) simulated seconds")
    ap.add_argument("--rpm", type=float, default=RPM, help="(--bearing-sim) shaft speed")
    ap.add_argument("--fs-sim", type=int, default=FS_SIM, help="(--bearing-sim) sampling rate of the simulated signal")
    ap.add_argument("--seed", type=int, default=42, help="(--bearing-sim) random seed")
    ap.add_argument("--bearing-vmd-k", type=int, default=BRG_VMD_K, help="(--bearing-sim) number of VMD modes")
    ap.add_argument("--components", action="store_true", help="(--bearing-sim) also open the components figure")
    # ---- bearing view on YOUR dataset (new) ----
    ap.add_argument("--bearing", action="store_true",
                    help="dataset playback + bearing view: BPFO/BPFI/BSF/FTF components combined in the time domain and marked in the FFT")
    ap.add_argument("--fr", type=float, default=None, help="(--bearing) shaft speed in Hz (default: estimated from the data)")
    ap.add_argument("--ftf", type=float, default=None, help="(--bearing) cage frequency in Hz (default 9.52, measured from your Faulty recordings)")
    ap.add_argument("--balls", type=int, default=None, help="(--bearing) number of balls; BPFO/BPFI are only used for the verdict when this is given")
    ap.add_argument("--ball-d", type=float, default=None, help="(--bearing) ball diameter in mm (with --pitch-d: exact geometry)")
    ap.add_argument("--pitch-d", type=float, default=None, help="(--bearing) pitch diameter in mm")
    ap.add_argument("--phi", type=float, default=0.0, help="(--bearing) contact angle in degrees")
    ap.add_argument("--band", nargs=2, type=float, default=list(DATA_BAND), help="(--bearing) resonance band in Hz")
    ap.add_argument("--fault-db", type=float, default=None, help="(--bearing) fault threshold in dB (default 8, calibrated on your Healthy files)")
    ap.add_argument("--reference", default=None, help="(--bearing) a Healthy Excel dataset: calibrates the threshold and is drawn as the Healthy spectrum")
    ap.add_argument("--fault-k", type=int, default=DATA_FAULT_K, help="(--bearing) number of VMD modes for the COMBINED fault signal (default 10)")
    ap.add_argument("--keep-base-analysis", action="store_true", help="(--bearing) also run your original raw-signal VMD + entropy window")
    args = ap.parse_args()
    PERM_ORDER = args.perm_order

    if args.bearing_sim:
        run_bearing_simulation(args)
        return

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
    if args.bearing:
        sim = make_bearing_dataset_simulator(args, data, fs)
    else:
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
