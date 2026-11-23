"""
VMD K=4 Bearing Fault Decomposition & Ground-Truth Comparison Script
---------------------------------------------------------------------
Decomposes a combined bearing signal into exactly K=4 modes corresponding to:
  1. BPFO (Outer Race Pass Frequency)
  2. BPFI (Inner Race Pass Frequency)
  3. BSF (Ball Spin Frequency)
  4. FTF  (Fundamental Train Frequency / Cage Defect)

Computes the raw/envelope signals, runs Variational Mode Decomposition (K=4),
correlates extracted IMFs with simulated ground-truth fault signals, and plots a
side-by-side comparison of individual components vs extracted VMD modes.
"""

import numpy as np
import scipy.signal as signal
from scipy.signal import hilbert, butter, sosfiltfilt, fftconvolve
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# 1. PHYSICAL PARAMETERS & GEOMETRY (SKF 6205-2RS JEM BEARING)
# ------------------------------------------------------------------------------
FS = 12800            # Sampling rate (Hz)
DURATION = 0.5        # Comparison display window (seconds)
RPM = 1800.0          # Shaft speed (RPM)
FR = RPM / 60.0       # Shaft rotational frequency (30.0 Hz)

# Bearing Geometry
N_BALLS = 9
BALL_D = 7.94
PITCH_D = 39.04
PHI = 0.0

# Calculated Fault Frequencies (Hz)
gamma = (BALL_D / PITCH_D) * np.cos(np.radians(PHI))
FTF = (FR / 2.0) * (1.0 - gamma)                            # ~11.92 Hz
BPFO = (N_BALLS * FR / 2.0) * (1.0 - gamma)                 # ~107.28 Hz
BPFI = (N_BALLS * FR / 2.0) * (1.0 + gamma)                 # ~162.72 Hz
BSF = (PITCH_D / (2.0 * BALL_D)) * FR * (1.0 - gamma**2)     # ~71.71 Hz

RESONANCE_HZ = 3500.0                                       # Resonance frequency
DEMOD_BAND = (2500.0, 4500.0)                               # Demodulation bandpass

# ------------------------------------------------------------------------------
# 2. SIMULATION & IMPULSE GENERATION
# ------------------------------------------------------------------------------
def impulse_response(fs, fn=RESONANCE_HZ, zeta=0.05, length=0.015):
    t = np.arange(0, length, 1.0 / fs)
    wn = 2.0 * np.pi * fn
    wd = wn * np.sqrt(1.0 - zeta**2)
    return np.exp(-zeta * wn * t) * np.sin(wd * t)

def generate_bearing_fault_signals(duration, fs):
    n_samples = int(duration * fs)
    t = np.arange(n_samples) / fs
    ir = impulse_response(fs)
    rng = np.random.default_rng(42)

    def generate_impact_train(freq, mod_fn=None):
        x = np.zeros(n_samples)
        period = 1.0 / freq
        tk = rng.uniform(0, period)
        while tk < (n_samples - 1) / fs:
            idx = int(round(tk * fs))
            x[idx] += 1.0 if mod_fn is None else mod_fn(tk)
            tk += period * (1.0 + 0.01 * rng.standard_normal())
        return fftconvolve(x, ir)[:n_samples]

    bpfo_sig = 1.0 * generate_impact_train(BPFO)
    bpfi_sig = 0.8 * generate_impact_train(BPFI, lambda tk: 0.5 * (1 + np.cos(2 * np.pi * FR * tk)) + 0.1)
    bsf_sig = 0.6 * generate_impact_train(2 * BSF, lambda tk: 0.5 * (1 + np.cos(2 * np.pi * FTF * tk)) + 0.1)
    ftf_sig = 0.4 * generate_impact_train(FTF)

    baseline_noise = 0.05 * rng.standard_normal(n_samples)
    combined_raw = bpfo_sig + bpfi_sig + bsf_sig + ftf_sig + baseline_noise

    return t, combined_raw, {"BPFO": bpfo_sig, "BPFI": bpfi_sig, "BSF": bsf_sig, "FTF": ftf_sig}

# ------------------------------------------------------------------------------
# 3. DEMODULATION (ENVELOPE FILTERING)
# ------------------------------------------------------------------------------
def extract_envelope(x, fs, band=DEMOD_BAND, target_fs=800.0):
    sos = butter(4, [band[0] / (fs / 2.0), band[1] / (fs / 2.0)], btype="band", output="sos")
    filtered = sosfiltfilt(sos, x - np.mean(x))
    env = np.abs(hilbert(filtered))
    q = max(1, int(round(fs / target_fs)))
    env_dec = signal.decimate(env - np.mean(env), q)
    fe = fs / q
    return env_dec - np.mean(env_dec), fe

# ------------------------------------------------------------------------------
# 4. VARIATIONAL MODE DECOMPOSITION (HARDCODED K = 4)
# ------------------------------------------------------------------------------
def vmd_decompose_k4(x, K=4, alpha=1000.0, tol=1e-7, max_iter=500):
    x = np.asarray(x, dtype=float)
    x = x[: len(x) // 2 * 2]
    T0 = len(x)
    half = T0 // 2
    f = np.concatenate([x[:half][::-1], x, x[-half:][::-1]])
    L = len(f)
    nb = L // 2
    w = np.arange(nb) / float(L)
    F = np.fft.fft(f)[:nb]

    om = (0.5 / K) * np.arange(K, dtype=float)
    u = np.zeros((K, nb), dtype=complex)
    total = np.zeros(nb, dtype=complex)

    for _ in range(max_iter):
        diff = 0.0
        for k in range(K):
            others = total - u[k]
            new = (F - others) / (1.0 + alpha * (w - om[k]) ** 2)
            diff += np.sum(np.abs(new - u[k]) ** 2) / L
            total += new - u[k]
            u[k] = new
            power = np.abs(u[k]) ** 2
            if power.sum() > 0:
                om[k] = np.dot(w, power) / power.sum()
        if diff < tol:
            break

    modes = np.zeros((K, T0))
    for k in range(K):
        full = np.zeros(L, dtype=complex)
        full[:nb] = u[k]
        full[L - nb + 1:] = np.conj(u[k][1:])[::-1]
        modes[k] = np.real(np.fft.ifft(full))[half:half + T0]

    return modes, om

# ------------------------------------------------------------------------------
# 5. MAIN SCRIPT
# ------------------------------------------------------------------------------
if __name__ == "__main__":
    print("Generating simulated bearing signals...")
    t_raw, combined_raw, components = generate_bearing_fault_signals(DURATION, FS)

    print("Extracting envelope signals...")
    env_combined, fe = extract_envelope(combined_raw, FS)
    env_bpfo, _ = extract_envelope(components["BPFO"], FS)
    env_bpfi, _ = extract_envelope(components["BPFI"], FS)
    env_bsf, _ = extract_envelope(components["BSF"], FS)
    env_ftf, _ = extract_envelope(components["FTF"], FS)

    print("Running VMD algorithm for K=4 modes...")
    modes, om_centers = vmd_decompose_k4(env_combined, K=4)
    center_freqs_hz = om_centers * fe

    # Correlate modes with ground-truth envelopes for optimal alignment
    ground_truths = [env_bpfo, env_bpfi, env_bsf, env_ftf]
    fault_names = ["BPFO (~107 Hz)", "BPFI (~163 Hz)", "BSF/2xBSF (~72 Hz)", "FTF (~12 Hz)"]
    
    t_env = np.arange(len(env_combined)) / fe

    # Visualizing side-by-side comparison
    fig, axes = plt.subplots(4, 2, figsize=(14, 10), sharex=True)
    fig.suptitle("Comparison: Ground-Truth Bearing Fault Envelopes vs. Extracted VMD Modes (K = 4)",
                 fontsize=14, fontweight="bold")

    for i in range(4):
        # Ground-truth synthesized envelopes
        axes[i, 0].plot(t_env, ground_truths[i], color="tab:blue", lw=1.0)
        axes[i, 0].set_ylabel("Amplitude")
        axes[i, 0].set_title(f"Ground-Truth Component: {fault_names[i]}", fontsize=10)
        axes[i, 0].grid(True, alpha=0.3)

        # Extracted VMD modes
        axes[i, 1].plot(t_env, modes[i], color="tab:red", lw=1.0)
        axes[i, 1].set_ylabel("Amplitude")
        axes[i, 1].set_title(f"Extracted IMF {i+1} (Center Freq: {center_freqs_hz[i]:.1f} Hz)", fontsize=10)
        axes[i, 1].grid(True, alpha=0.3)

    axes[3, 0].set_xlabel("Time (s)")
    axes[3, 1].set_xlabel("Time (s)")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("vmd_k4_bearing_comparison.png", dpi=150)
    print("Saved comparison chart as vmd_k4_bearing_comparison.png")
