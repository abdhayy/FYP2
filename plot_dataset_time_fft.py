"""
Plot Time Domain + FFT from a saved vibration dataset
------------------------------------------------------------
Reads the Excel file written by collect_vibration_dataset.m
(a "Time (s)" column + one "Group_N" column per 5-second window)
and displays:

  Figure 1: Time-domain waveform (raw + smoothed) and FFT spectrum
            of ONE chosen group, with the fundamental and its
            2nd/3rd harmonics labelled (same method as the MATLAB
            monitor: Welch-averaged, Hann-windowed spectrum, and
            Harmonic Product Spectrum for the fundamental).
  Figure 2: Average FFT spectrum across ALL groups in the file
            (more stable than a single 5 s window - useful when
            comparing Healthy vs Faulty files).

Requirements:  pip install numpy scipy pandas matplotlib openpyxl

Usage examples:
  python plot_dataset_time_fft.py
  python plot_dataset_time_fft.py --file Vibration_Dataset_Faulty_2000Hz.xlsx --group 12
  python plot_dataset_time_fft.py --save            # save PNGs instead of opening windows
  python plot_dataset_time_fft.py --xmax 300        # zoom the FFT to 0-300 Hz
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
from scipy.signal.windows import hann

# ---------------- USER-ADJUSTABLE DEFAULTS ----------------
DATASET_FILE   = "Vibration_Dataset_Healthy_2000Hz.xlsx"  # <-- file to plot
GROUP_TO_PLOT  = 1        # <-- which group (1 = Group_1, 2 = Group_2, ...)
SMOOTH_CUTOFF  = 150.0    # Hz, low-pass for the DISPLAYED smooth trace only
N_SEG_WELCH    = 6        # Welch segments averaged for the FFT
PEAK_RATIO_MIN = 5.0      # only label harmonics if the fundamental stands this many
                          # times above the median spectrum level (avoids labelling
                          # noise as "harmonics" when the machine was idle)
# ----------------------------------------------------------


def welch_amplitude_spectrum(x, fs, n_seg=6):
    """Hann-windowed, ~50%-overlap, averaged amplitude spectrum (units of x)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    seg_len = min(max(int(np.floor(2 * n / (n_seg + 1))), 16), n)
    step = max(seg_len // 2, 1)
    win = hann(seg_len)
    win_gain = win.mean()

    starts = range(0, n - seg_len + 1, step)
    half = seg_len // 2 + 1
    accum = np.zeros(half)
    count = 0
    for s in starts:
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
    # sits at bin h*j. Truncate every factor to the same length (min_len) so
    # no part of the product is left un-multiplied.
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


def load_dataset(path):
    df = pd.read_excel(path)
    time = df.iloc[:, 0].to_numpy(dtype=float)
    fs = int(round(1.0 / np.mean(np.diff(time))))
    group_names = list(df.columns[1:])
    data = df.iloc[:, 1:].to_numpy(dtype=float)   # samples x groups
    return time, fs, group_names, data


def plot_single_group(time, fs, name, x, save=False, xmax=None):
    x = x - x.mean()
    smooth = smooth_for_display(x, fs, SMOOTH_CUTOFF)
    freqs, p1 = welch_amplitude_spectrum(x, fs, N_SEG_WELCH)

    fig, (ax_t, ax_f) = plt.subplots(1, 2, figsize=(14, 5.5))

    ax_t.plot(time, x, color="lightsteelblue", lw=0.6, label="Raw")
    ax_t.plot(time, smooth, color="blue", lw=1.0, label=f"Smoothed ({SMOOTH_CUTOFF:.0f} Hz LP)")
    ax_t.set_xlabel("Time (s)")
    ax_t.set_ylabel("Amplitude (g)")
    ax_t.set_title(f"Time Domain Waveform - {name} (Fs = {fs} Hz)")
    ax_t.set_xlim(time[0], time[-1])
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
    f0_amp = p1[np.argmin(np.abs(freqs - f0))]
    if f0_amp > PEAK_RATIO_MIN * np.median(p1):
        for k in (1, 2, 3):
            fk = f0 * k
            if fk <= freqs[-1]:
                ak = p1[np.argmin(np.abs(freqs - fk))]
                ax_f.plot(fk, ak, "kv", markerfacecolor="yellow", markersize=8)
                ax_f.annotate(f"Mode {k}: {fk:.1f} Hz", (fk, ak),
                              textcoords="offset points",
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


def plot_average_spectrum(fs, group_names, data, save=False):
    n_groups = data.shape[1]
    spectra = []
    for g in range(n_groups):
        freqs, p1 = welch_amplitude_spectrum(data[:, g], fs, N_SEG_WELCH)
        spectra.append(p1)
    spectra = np.array(spectra)
    mean_p1 = spectra.mean(axis=0)
    std_p1 = spectra.std(axis=0)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.plot(freqs, mean_p1, color="darkred", lw=1.3, label=f"Mean of {n_groups} groups")
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


def main():
    parser = argparse.ArgumentParser(description="Plot time domain + FFT from a saved dataset.")
    parser.add_argument("--file", default=DATASET_FILE, help="Excel dataset to load")
    parser.add_argument("--group", type=int, default=GROUP_TO_PLOT, help="group number (1-based)")
    parser.add_argument("--save", action="store_true", help="save PNGs instead of opening windows")
    parser.add_argument("--xmax", type=float, default=None, help="zoom the FFT x-axis to this max frequency (Hz)")
    args = parser.parse_args()

    print(f"Loading {args.file} ... (large Excel files can take a little while)")
    time, fs, group_names, data = load_dataset(args.file)
    print(f"Loaded {data.shape[1]} groups x {data.shape[0]} samples, inferred Fs = {fs} Hz")

    if not 1 <= args.group <= data.shape[1]:
        raise SystemExit(f"--group must be between 1 and {data.shape[1]}")

    plot_single_group(time, fs, group_names[args.group - 1], data[:, args.group - 1], args.save, args.xmax)
    plot_average_spectrum(fs, group_names, data, args.save)

    if not args.save:
        plt.show()


if __name__ == "__main__":
    main()
