"""
Bearing Fault Simulator  -  BPFO + BPFI + BSF + FTF combined  (Time Domain / FFT / VMD / Entropy)
------------------------------------------------------------
Built on plot_dataset_time_fft_entropy.py (its VMD and entropy code is reused
unchanged). Instead of replaying an Excel dataset, this tool SIMULATES a
rolling-element bearing vibration signal that contains all four characteristic
defect signals at once and plays it back live:

    BPFO  outer-race defect      impacts at BPFO, constant amplitude
    BPFI  inner-race defect      impacts at BPFI, amplitude modulated at shaft speed fr
    BSF   rolling-element defect impacts at 2 x BSF, modulated at the cage speed FTF
    FTF   cage defect            FTF sine components + impacts at FTF

Signal model = the reference generator you gave me (SKF 6205-2RS JEM, 1800 rpm):
each impact rings a structural resonance (3.5 kHz, 5 % damping), 1 % random slip
on the impact timing, shaft harmonics 1x/2x/3x and Gaussian noise are added.
Default scenario: BPFO is the dominant fault (weight 1.0), the other three are
present but weaker (0.3).

WINDOWS
  1  Live monitor
       Time domain  : the combined signal (rolling 0.5 s so the impacts are visible)
       FFT (raw)    : spectrum of the acceleration - the energy sits at the resonance;
                      the bearing frequencies are NOT visible here
       FFT (envelope): spectrum of the demodulated signal (band-pass around the
                      resonance -> Hilbert envelope -> FFT). THIS is where BPFO,
                      BPFI, BSF and FTF show up; each one is marked on the graph.
       Banner       : the diagnosed condition, e.g.  "CONDITION: BPFO FAULTY",
                      computed live from the line strength (dB above the local noise
                      floor) at each fault's frequencies - nothing is hard-coded.
  2  VMD + entropy analysis, refreshed every 10 s of playback
       The demodulated signal is split into K modes with VMD. Each mode is labelled
       with the fault line that carries its largest spectral peak, then the modes
       are RANKED by fault relevance (fault-labelled modes first, strongest RMS
       first). So "IMF1" is the strongest fault mode - with the default scenario
       that is the BPFO mode (centre frequency ~107.5 Hz). All seven entropies
       are shown for the envelope and for the top modes (heat-map table).
  3  (--components) the four fault components and their sum, stacked.

NOTE ON IMF ORDER: here IMFs are ranked by fault relevance, not sorted by centre
frequency as in the earlier scripts. VMD itself is unchanged.
NOTE: VMD is the numpy implementation from the attached file, not MATLAB's vmd().

Keyboard: SPACE = pause / resume      R = restart

Requirements:  pip install numpy scipy pandas matplotlib openpyxl

Examples:
  python bearing_fault_simulator.py                         # BPFO-dominant, 60 s, real time
  python bearing_fault_simulator.py --speed 5               # 5x faster
  python bearing_fault_simulator.py --dominant BPFI         # inner-race dominant
  python bearing_fault_simulator.py --weights BPFO=1,BPFI=1,BSF=1,FTF=1
  python bearing_fault_simulator.py --healthy               # no faults -> "HEALTHY"
  python bearing_fault_simulator.py --components --log sim_results.xlsx
"""
import argparse
import ast
import math
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.signal import butter, sosfiltfilt, hilbert, decimate, find_peaks, fftconvolve
from scipy.signal.windows import hann
from scipy.spatial import cKDTree
from scipy.stats import norm

# ---------------- USER-ADJUSTABLE DEFAULTS ----------------
FS_SIM        = 12800            # sampling rate of the simulated signal (Hz)
BEARING       = dict(n=9, d=7.94, D=39.04, phi=0.0)   # SKF 6205-2RS JEM: balls, ball dia, pitch dia (mm), contact angle (deg)
RPM           = 1800.0           # shaft speed
RESONANCE_HZ  = 3500.0           # structural resonance excited by the impacts
BAND          = (2500.0, 4500.0) # demodulation band around the resonance
DAMPING       = 0.05
SLIP          = 0.01             # 1 % random slip
NOISE_STD     = 0.05             # g
DURATION_S    = 60.0
DEFAULT_WEIGHTS = dict(BPFO=1.0, BPFI=0.3, BSF=0.3, FTF=0.3)

DISPLAY_SECONDS  = 0.5           # width of the rolling time-domain display
SPEC_SECONDS     = 5.0           # data used for the live FFT / envelope spectrum
SPEC_EVERY       = 0.5           # spectra refresh interval (s of playback)
ANALYSIS_SECONDS = 10.0          # VMD + entropy analysis interval
VMD_K            = 10            # number of VMD modes
VMD_HALFWIDTH_HZ = 5.0           # VMD mode bandwidth (alpha is derived from it)
ENV_FS_TARGET    = 800.0         # sampling rate of the demodulated envelope
ENV_LOWPASS_HZ   = 300.0         # envelope low-pass (keeps the fault lines, drops the rest)
N_ENTROPY_MODES  = 6             # entropies are computed for the envelope + this many top modes
FAULT_DB         = 12.0          # line strength above the noise floor that declares a fault
SECONDARY_DB     = 9.0           # ... and that lists a secondary component
FAULT_COLORS = {"BPFO": "tab:red", "BPFI": "tab:orange", "BSF": "tab:purple", "FTF": "tab:green"}
# ----------------------------------------------------------

# ===== reused UNCHANGED from plot_dataset_time_fft_entropy.py (VMD + entropy measures) =====
N_SEG_WELCH    = 6
VMD_ALPHA        = 1000.0
EN_M          = 2
EN_R          = 0.2
FUZZY_N       = 2
PERM_ORDER    = 5
PERM_DELAY    = 1
DISP_M        = 2
DISP_C        = 6
SHANNON_BINS  = 32
MSE_MAX_SCALE = 10
FUZZY_MAX_ROWS = 1500


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
def analyze_window_bearing(window, fs, F, band, K, halfwidth_hz, n_modes, t_start):
    t0 = time.perf_counter()
    env, fe = demodulate(window, fs, band)
    f, a = welch_amplitude_spectrum(env, fe, 3)
    evidence = diagnose(f, a, F)

    modes, om, n_iter = vmd_decompose(env, K, alpha=(fe / halfwidth_hz) ** 2)
    centers = om * fe
    rms = np.sqrt(np.mean(modes ** 2, axis=1))
    labels = label_modes(modes, fe, F)
    # rank: fault-labelled modes first, then by RMS (strongest first)
    order = sorted(range(K), key=lambda i: (labels[i][0] in (None, "SHAFT"), -rms[i]))
    modes, centers, rms = modes[order], centers[order], rms[order]
    labels = [labels[i] for i in order]

    names = ["Env"] + [f"IMF{k + 1}" for k in range(min(n_modes, K))]
    series = [env] + [modes[k] for k in range(min(n_modes, K))]
    matrix = np.full((len(ENTROPY_DEFS), len(series)), np.nan)
    for j, s in enumerate(series):
        vals = compute_all_entropies(s, fe)
        for i, (key, _, _) in enumerate(ENTROPY_DEFS):
            matrix[i, j] = vals[key]
    mse = multiscale_sample_entropy(env, MSE_MAX_SCALE, EN_M)

    return dict(t_start=t_start, t_end=t_start + len(window) / fs, fs=fs, fe=fe, K=K, env=env,
                spec_f=f, spec_a=a, evidence=evidence, modes=modes, centers=centers, rms=rms,
                labels=labels, names=names, matrix=matrix, mse=mse, mse_ci=float(np.nansum(mse)),
                n_iter=n_iter, elapsed=time.perf_counter() - t0)


def result_to_row(res, label, verdict):
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


def fmt_clock(seconds):
    seconds = int(round(seconds))
    return f"{seconds // 60}:{seconds % 60:02d}"


def mark_fault_lines(ax, F, xmax, ymax, label_size=8, show_text=True):
    """Draw the BPFO / BPFI / BSF / FTF markers on an envelope-spectrum axis."""
    fr = F["fr"]
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
                ax.text(fq, ymax * (0.98 if ls == "-" else 0.80), f" {name} {fq:.1f}", color=FAULT_COLORS[fault],
                        rotation=90, va="top", ha="right", fontsize=label_size)
    ax.axvline(fr, color="gray", ls="-.", lw=0.8)
    if show_text:
        ax.text(fr, ymax * 0.55, f" fr {fr:.1f}", color="dimgray", rotation=90, va="top", ha="right", fontsize=label_size)


# ===================== live simulator =====================
class BearingSimulator:
    """Plays the simulated bearing signal back live: time domain, raw FFT, envelope spectrum with the
    four fault markers, live diagnosis banner and (every 10 s) the VMD + entropy analysis window."""

    def __init__(self, x, fs, F, comps=None, speed=1.0, band=BAND, display_seconds=DISPLAY_SECONDS,
                 spec_seconds=SPEC_SECONDS, spec_every=SPEC_EVERY, analysis=True,
                 analysis_seconds=ANALYSIS_SECONDS, vmd_k=VMD_K, label="Simulated bearing",
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
        self.y_lim, self.raw_top, self.env_top = 0.5, 1.0, 0.05

        self.fig = plt.figure(figsize=(15, 9))
        gs = self.fig.add_gridspec(2, 2, height_ratios=[1, 1.25], left=0.06, right=0.985, top=0.80,
                                   bottom=0.075, hspace=0.38, wspace=0.18)
        self.ax_t = self.fig.add_subplot(gs[0, 0])
        self.ax_r = self.fig.add_subplot(gs[0, 1])
        self.ax_e = self.fig.add_subplot(gs[1, :])

        t_axis = np.arange(1, self.dw + 1) / fs
        self.line_t, = self.ax_t.plot(t_axis, np.zeros(self.dw), color="tab:blue", lw=0.6)
        self.ax_t.set_xlim(t_axis[0], t_axis[-1])
        self.ax_t.set_ylim(-self.y_lim, self.y_lim)
        self.ax_t.set_xlabel("Time (s)")
        self.ax_t.set_ylabel("Acceleration (g)")
        self.ax_t.set_title("TIME DOMAIN - BPFO + BPFI + BSF + FTF combined (rolling "
                            f"{display_seconds:g} s)", fontsize=10)
        self.ax_t.grid(True, alpha=0.4)

        self.line_r, = self.ax_r.plot([0, fs / 2], [0, 0], color="red", lw=0.8)
        self.ax_r.axvspan(band[0], min(band[1], fs / 2), color="gold", alpha=0.25,
                          label=f"demodulation band {band[0] / 1000:.1f}-{band[1] / 1000:.1f} kHz")
        self.ax_r.set_xlim(0, fs / 2)
        self.ax_r.set_ylim(0, self.raw_top)
        self.ax_r.set_xlabel("Frequency (Hz)")
        self.ax_r.set_ylabel("Amplitude (g)")
        self.ax_r.set_title("FFT of the raw signal - the energy sits at the resonance; "
                            "the bearing frequencies are not visible here", fontsize=10)
        self.ax_r.legend(loc="upper right", fontsize=8)
        self.ax_r.grid(True, alpha=0.4)

        self.line_e, = self.ax_e.plot([0, xmax_env], [0, 0], color="black", lw=0.9)
        self.ax_e.set_xlim(0, xmax_env)
        self.ax_e.set_ylim(0, self.env_top)
        mark_fault_lines(self.ax_e, F, xmax_env, self.env_top)
        self.ax_e.set_xlabel("Frequency (Hz)")
        self.ax_e.set_ylabel("Envelope amplitude (g)")
        self.ax_e.set_title("FFT OF THE ENVELOPE (demodulated signal) - where BPFO, BPFI, BSF and FTF appear",
                            fontsize=10, fontweight="bold")
        self.ax_e.grid(True, alpha=0.4)

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
        buf = np.zeros(self.dw)
        seg = self.x[max(0, n - self.dw):n]
        if len(seg):
            buf[self.dw - len(seg):] = seg
        self.line_t.set_ydata(buf)
        target = max(1.15 * float(np.max(np.abs(buf))), 0.3)
        self.y_lim = 0.8 * self.y_lim + 0.2 * target
        self.ax_t.set_ylim(-self.y_lim, self.y_lim)

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
            self.line_r.set_data([0, self.fs / 2], [0, 0])
            self.line_e.set_data([0, self.xmax_env], [0, 0])
            self.banner.set_text("CONDITION: ACQUIRING...")
            self.banner.set_color("gray")
            self.sub.set_text("")
            return
        seg = self.x[n - have:n]
        fr_, ar_ = welch_amplitude_spectrum(seg, self.fs, 6)
        self.line_r.set_data(fr_, ar_)
        self.raw_top = 0.8 * self.raw_top + 0.2 * max(1.25 * float(ar_.max()), 1e-3)
        self.ax_r.set_ylim(0, self.raw_top)

        env, fe = demodulate(seg, self.fs, self.band)
        fq, ae = welch_amplitude_spectrum(env, fe, 3)
        m = fq <= self.xmax_env
        self.line_e.set_data(fq[m], ae[m])
        top = max(1.25 * float(ae[m & (fq > 3)].max()), 1e-3)
        self.env_top = 0.8 * self.env_top + 0.2 * top
        self.ax_e.set_ylim(0, self.env_top)
        for t in list(self.ax_e.texts):
            t.remove()
        for ln in list(self.ax_e.lines)[1:]:
            ln.remove()
        mark_fault_lines(self.ax_e, self.F, self.xmax_env, self.env_top)

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
            self.rows.append(result_to_row(res, self.label, verdict))
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
                              left=0.045, right=0.985, top=0.83, bottom=0.07)
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
                       "numbers = actual values)", fontsize=10, pad=40)

        # ---- right column: envelope spectrum with the VMD mode centres, and MSE ----
        gs_r = gs[:, 2].subgridspec(2, 1, hspace=0.55)
        ax_s = fig.add_subplot(gs_r[0])
        f, a = res["spec_f"], res["spec_a"]
        m = f <= 450
        ax_s.plot(f[m], a[m], color="black", lw=0.8)
        ymax = float(a[m & (f > 3)].max()) * 1.3
        ax_s.set_ylim(0, ymax)
        mark_fault_lines(ax_s, self.F, 450, ymax, label_size=6, show_text=False)
        for k in range(min(n_show, K)):
            fault = res["labels"][k][0]
            ax_s.plot(res["centers"][k], ymax * 0.93, "v", color=FAULT_COLORS.get(fault, "gray"), ms=6)
            ax_s.text(res["centers"][k], ymax * 0.99, f"{k + 1}", ha="center", fontsize=7)
        ax_s.set_xlabel("Frequency (Hz)", fontsize=9)
        ax_s.set_ylabel("Envelope amplitude (g)", fontsize=9)
        ax_s.set_title("Envelope spectrum; numbers = IMF centre frequencies", fontsize=9)
        ax_s.grid(True, alpha=0.4)
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
        fig.text(0.5, 0.905, "evidence (dB above noise floor):  " +
                 "   ".join(f"{k} {v:.0f}" for k, v in res["evidence"].items()) + extra +
                 f"   |   window {res['t_start']:.0f}-{res['t_end']:.0f} s   |   analysis #{self.n_results}"
                 f"   |   computed in {res['elapsed']:.1f} s", ha="center", fontsize=10, color="dimgray")
        fig.text(0.5, 0.868, "IMFs are ranked by fault relevance: IMF1 is the strongest fault mode "
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
        env, fe = demodulate(sig[: int(10 * fs)], fs, band)
        f, a = welch_amplitude_spectrum(env, fe, 3)
        m = f <= 450
        ax[i, 1].plot(f[m], a[m], color="black", lw=0.8)
        ymax = float(a[m & (f > 3)].max()) * 1.3
        ax[i, 1].set_ylim(0, ymax)
        mark_fault_lines(ax[i, 1], F, 450, ymax, label_size=7, show_text=(key == "ALL"))
        ax[i, 1].set_title("envelope spectrum", fontsize=9)
    ax[-1, 0].set_xlabel("Time (s)")
    ax[-1, 1].set_xlabel("Frequency (Hz)")
    fig.suptitle("Fault components and their combination", fontweight="bold")
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=130)
    return fig


# ===================== main =====================
def parse_weights(text):
    w = {}
    for part in text.split(","):
        k, v = part.split("=")
        w[k.strip().upper()] = float(v)
    return w


def main():
    ap = argparse.ArgumentParser(description="Bearing fault simulator: BPFO + BPFI + BSF + FTF combined.")
    ap.add_argument("--duration", type=float, default=DURATION_S, help="simulated seconds (default 60)")
    ap.add_argument("--speed", type=float, default=1.0, help="playback speed (1 = real time)")
    ap.add_argument("--dominant", choices=["BPFO", "BPFI", "BSF", "FTF"], default="BPFO",
                    help="fault with weight 1.0; the other three get 0.3 (default BPFO)")
    ap.add_argument("--weights", default=None, help="explicit weights, e.g. BPFO=1,BPFI=0.3,BSF=0.3,FTF=0.3")
    ap.add_argument("--healthy", action="store_true", help="no faults (baseline only)")
    ap.add_argument("--rpm", type=float, default=RPM)
    ap.add_argument("--fs", type=int, default=FS_SIM, help="sampling rate of the simulated signal")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--vmd-k", type=int, default=VMD_K)
    ap.add_argument("--no-analysis", action="store_true", help="disable the VMD + entropy window")
    ap.add_argument("--components", action="store_true", help="also show the four components and their sum")
    ap.add_argument("--log", default=None, help="Excel file for the VMD / entropy results")
    ap.add_argument("--interval", type=int, default=100, help="screen refresh interval (ms)")
    args = ap.parse_args()

    F = bearing_frequencies(BEARING["n"], BEARING["d"], BEARING["D"], BEARING["phi"], args.rpm)
    if args.healthy:
        weights = {k: 0.0 for k in ("BPFO", "BPFI", "BSF", "FTF")}
    elif args.weights:
        weights = parse_weights(args.weights)
    else:
        weights = {k: (1.0 if k == args.dominant else 0.3) for k in ("BPFO", "BPFI", "BSF", "FTF")}

    print("Bearing: SKF 6205-2RS JEM (9 balls, d = 7.94 mm, D = 39.04 mm, phi = 0 deg)")
    print(f"Shaft {args.rpm:.0f} rpm (fr = {F['fr']:.2f} Hz):  FTF {F['FTF']:.2f}  BSF {F['BSF']:.2f} "
          f"(ball defect at 2xBSF = {2 * F['BSF']:.2f})  BPFO {F['BPFO']:.2f}  BPFI {F['BPFI']:.2f} Hz")
    print("Fault weights:", weights)
    print(f"Simulating {args.duration:g} s at {args.fs} Hz ...")
    t, x, comps = simulate_bearing(args.duration, args.fs, F, weights, seed=args.seed)

    if args.components:
        plot_components(comps, args.fs, F)
    sim = BearingSimulator(x, args.fs, F, comps, speed=args.speed, analysis=not args.no_analysis,
                           vmd_k=args.vmd_k, log_file=args.log,
                           label="healthy" if args.healthy else "weights " + ",".join(f"{k}={v:g}" for k, v in weights.items()))
    print(f"Starting playback - about {fmt_clock(args.duration / args.speed)} (min:sec). SPACE = pause, R = restart.")
    sim.run(interval_ms=args.interval)


if __name__ == "__main__":
    main()
