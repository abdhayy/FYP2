"""
Vibration Dataset Playback Simulator  (Time Domain + FFT)
------------------------------------------------------------
Replays a saved dataset (the Excel file written by
collect_vibration_dataset.m: a "Time (s)" column + one "Group_N"
column per 5-second window) as if it were being acquired live.

The groups were cut from ONE continuous, gapless capture, so they
are joined back together in order and "streamed" at the original
sampling rate. A 5-minute dataset therefore plays for 5 minutes
(use --speed to go faster). Just like the MATLAB real-time monitor:

  * the display is a rolling 5-second window - it starts flat at
    zero and fills up as data "arrives"
  * Time Domain: raw trace + smoothed trace (150 Hz low-pass)
  * FFT: Welch-averaged, Hann-windowed spectrum, with the
    fundamental and its 2nd/3rd harmonics labelled (Harmonic
    Product Spectrum) only when a clear fundamental exists
  * status text: VIBRATION DETECTED / IDLE

Keyboard:  SPACE = pause / resume      R = restart

Requirements:  pip install numpy scipy pandas matplotlib openpyxl

Usage examples:
  python plot_dataset_time_fft.py                          # 5-minute playback
  python plot_dataset_time_fft.py --file Vibration_Dataset_Faulty_2000Hz.xlsx
  python plot_dataset_time_fft.py --speed 10               # 10x faster (30 s)
  python plot_dataset_time_fft.py --start-group 20         # begin at Group_20
  python plot_dataset_time_fft.py --xmax 300               # zoom FFT to 0-300 Hz
  python plot_dataset_time_fft.py --static --group 5        # old still plots instead
"""
import argparse
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.signal import butter, filtfilt
from scipy.signal.windows import hann

# ---------------- USER-ADJUSTABLE DEFAULTS ----------------
DATASET_FILE   = "Vibration_Dataset_Healthy_2000Hz.xlsx"  # <-- file to play back
SMOOTH_CUTOFF  = 150.0    # Hz, low-pass for the DISPLAYED smooth trace only
N_SEG_WELCH    = 6        # Welch segments averaged for the FFT
PEAK_RATIO_MIN = 5.0      # a fundamental must stand this many times above the median
                          # spectrum level to be labelled (avoids labelling noise as
                          # "harmonics" when the machine is idle)
MIN_SECONDS_FOR_FFT = 0.5 # need at least this much data before an FFT is drawn
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
    """Replays a continuous stream as a rolling window, like the live monitor."""

    def __init__(self, data, fs, xmax=None, speed=1.0):
        self.fs = fs
        self.group_len = data.shape[0]                 # samples per 5 s group
        self.num_groups = data.shape[1]
        # Groups are consecutive windows of one gapless capture -> join in order.
        stream = data.T.reshape(-1)
        # Remove the sensor's constant zero-g offset once (like the Teensy's
        # calibrated bias), so "no vibration" really sits at zero.
        self.stream = stream - np.median(stream)
        self.total = len(self.stream)
        self.W = self.group_len                        # rolling window = one group
        self.speed = speed
        self.xmax = xmax if xmax else fs / 2

        cutoff = min(SMOOTH_CUTOFF, 0.45 * fs)
        self.b, self.a = butter(6, cutoff / (fs / 2), btype="low")

        self.y_lim = 0.02            # smoothed auto-scale state (time axis)
        self.fft_top = 1e-3          # smoothed auto-scale state (FFT axis)
        self.labels = []
        self.ani = None
        self._tick = None

        t_axis = np.arange(1, self.W + 1) / fs
        self.fig, (self.ax_t, self.ax_f) = plt.subplots(1, 2, figsize=(14, 6.2))

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

        self.fig.tight_layout(rect=(0, 0, 1, 0.94))
        self.render(0)

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

    # ---- real-time driver (wall-clock based, so 5 min really is 5 min) ----
    def run(self, start_group=1, interval_ms=100):
        state = dict(t0=time.perf_counter(), paused=False, pause_t=0.0,
                     paused_total=0.0, done=False,
                     start_n=(max(1, start_group) - 1) * self.group_len)

        def tick(_frame):
            if state["paused"] or state["done"]:
                return
            elapsed = (time.perf_counter() - state["t0"] - state["paused_total"]) * self.speed
            n = min(self.total, state["start_n"] + int(elapsed * self.fs))
            self.render(n)
            if n >= self.total:
                state["done"] = True
                self.status.set_text("SIMULATION COMPLETE")
                self.status.set_color("navy")
                self.ani.event_source.stop()

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
    sim = PlaybackSimulator(data, fs, xmax=args.xmax, speed=args.speed)
    remaining = (data.shape[1] - args.start_group + 1) * data.shape[0] / fs / args.speed
    print(f"Starting playback - about {fmt_clock(remaining)} (min:sec) at x{args.speed:g} speed. "
          f"SPACE = pause/resume, R = restart.")
    sim.run(start_group=args.start_group, interval_ms=args.interval)


if __name__ == "__main__":
    main()
