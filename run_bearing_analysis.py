"""
Bearing Vibration & Fault Diagnostics Script (VMD + Hilbert Envelope Demodulation + Entropy Analysis)
------------------------------------------------------------------------------------------------------
This standalone script processes raw vibration data (or simulated signals) to perform:
1. Time Domain & FFT Monitoring (Raw + LP Smoothed + Welch FFT spectrum)
2. Envelope Demodulation (Bandpass resonance filtering + Hilbert transform)
3. Variational Mode Decomposition (VMD) into K Intrinsic Mode Functions (IMFs)
4. Multi-Entropy Analysis (SampEn, ApEn, FuzzyEn, PermEn, DispEn, SpecEn, ShanEn, MSE)
5. Automated Bearing Fault Detection (BPFO, BPFI, BSF, FTF vs noise floor thresholds)

Requirements:
    pip install numpy scipy pandas matplotlib openpyxl
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

# User Adjustable Defaults
DATASET_FILE = "Vibration_Dataset_Healthy_2000Hz.xlsx"
SMOOTH_CUTOFF = 150.0
N_SEG_WELCH = 6
PEAK_RATIO_MIN = 5.0
MIN_SECONDS_FOR_FFT = 0.5
ANALYSIS_SECONDS = 10.0
USE_VMD = True
VMD_K = 7
VMD_ALPHA = 1000.0
ANALYSIS_FS = 500

EN_M = 2
EN_R = 0.2
FUZZY_N = 2
PERM_ORDER = 5
PERM_DELAY = 1
DISP_M = 2
DISP_C = 6
SHANNON_BINS = 32
MSE_MAX_SCALE = 10
FUZZY_MAX_ROWS = 6000

def welch_amplitude_spectrum(x, fs, n_seg=6):
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
    x = np.asarray(x, dtype=float)
    x = x[: len(x) // 2 * 2]
    x = x - x.mean()
    T0 = len(x)
    half = T0 // 2
    f = np.concatenate([x[:half][::-1], x, x[-half:][::-1]])
    L = len(f)
    nb = L // 2
    w = np.arange(nb) / L
    F = np.fft.fft(f)[:nb]

    if init == "peaks":
        mag = np.abs(F)
        sm = np.convolve(mag, np.ones(21) / 21.0, mode="same")
        pk, _ = find_peaks(sm, distance=max(nb // (8 * K), 1))
        pk = pk[np.argsort(sm[pk])[::-1]][:K]
        om = np.sort(w[pk])
        if len(om) < K:
            extra = (0.5 / K) * np.arange(K - len(om)) + 0.5 / (2 * K)
            om = np.sort(np.concatenate([om, extra]))
    else:
        om = (0.5 / K) * np.arange(K)
    om = om.astype(float)

    u = np.zeros((K, nb), dtype=complex)
    lam = np.zeros(nb, dtype=complex)
    total = np.zeros(nb, dtype=complex)
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
                om[k] = np.dot(w, p) / s
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
        full[L - nb + 1:] = np.conj(u[k][1:])[::-1]
        modes[k] = np.real(np.fft.ifft(full))[half:half + T0]
    return modes, om, n_iter

def main():
    parser = argparse.ArgumentParser(description="Bearing Vibration Analysis & Diagnostic Simulator")
    parser.add_argument("--file", default=DATASET_FILE, help="Excel dataset file path")
    parser.add_argument("--bearing", action="store_true", help="Enable bearing diagnostic mode")
    parser.add_argument("--bearing-sim", action="store_true", help="Enable synthetic bearing fault simulation")
    parser.add_argument("--speed", type=float, default=1.0, help="Playback speed multiplier")
    args = parser.parse_args()

    print("=== Bearing Vibration & VMD Diagnostic Tool ===")
    print(f"Dataset file: {args.file}")
    print(f"Bearing diagnostic mode: {args.bearing}")
    print(f"Simulation mode: {args.bearing-sim}")
    print("For full options, run: python run_bearing_analysis.py --help")

if __name__ == "__main__":
    main()
