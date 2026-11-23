"""
Synthetic rolling-element bearing vibration signal generator
============================================================
Generates acceleration signals (in g) containing the four characteristic
bearing defect frequencies: BPFO, BPFI, BSF and FTF. The signals can be used
as sample "measured" vibration data to test envelope analysis, fault
detection or AI models.

Bearing: SKF 6205-2RS JEM deep-groove ball bearing (same type as the CWRU
bearing data set, drive end), so you can cross-check the results.

Signal model (impulse-response / amplitude-modulation model):
  - Each defect produces a train of impacts at its defect frequency.
  - Each impact excites a structural resonance (decaying sinusoid).
  - Outer race (fixed, load zone) : no modulation.
  - Inner race (rotates with shaft): amplitude modulated at shaft speed fr.
  - Rolling element (spins, orbits): impacts at 2*BSF, modulated at FTF.
  - Cage (FTF)                     : low-frequency component at FTF and
                                     its harmonics + modulation at FTF.
  - Random slip (1-2 %) jitter on impact timing, as in real bearings.
  - Shaft components (1x, 2x, 3x) and Gaussian noise added.
"""

import numpy as np

# -------------------------------------------------------------------------
# 1. BEARING INFORMATION (for analytical calculation)
# -------------------------------------------------------------------------
BEARING = {
    "designation": "SKF 6205-2RS JEM (deep-groove ball bearing)",
    "inner_diameter_mm": 25.0,
    "outer_diameter_mm": 52.0,
    "width_mm": 15.0,
    "n_balls": 9,                  # Nb
    "ball_diameter_mm": 7.94,      # d  (0.3126 in)
    "pitch_diameter_mm": 39.04,    # D  (1.537 in)
    "contact_angle_deg": 0.0,      # phi
}

SHAFT_RPM = 1800.0                 # outer race fixed, inner race rotating
FS = 25600                         # sampling frequency, Hz
DURATION = 10.0                    # s
RESONANCE_HZ = 3500.0              # excited structural resonance
DAMPING_RATIO = 0.05
SLIP = 0.01                        # 1 % random slip
NOISE_STD = 0.05                   # g
SEED = 42


def defect_frequencies(b=BEARING, rpm=SHAFT_RPM):
    """Analytical bearing defect frequencies (Hz)."""
    fr = rpm / 60.0
    Nb = b["n_balls"]
    d = b["ball_diameter_mm"]
    D = b["pitch_diameter_mm"]
    c = np.cos(np.radians(b["contact_angle_deg"]))
    r = d / D * c
    return {
        "fr":   fr,
        "FTF":  fr / 2 * (1 - r),
        "BPFO": Nb * fr / 2 * (1 - r),
        "BPFI": Nb * fr / 2 * (1 + r),
        "BSF":  D / (2 * d) * fr * (1 - r ** 2),
    }


# -------------------------------------------------------------------------
# 2. SIGNAL BUILDING BLOCKS
# -------------------------------------------------------------------------
def impulse_response(fs, fn=RESONANCE_HZ, zeta=DAMPING_RATIO, length=0.02):
    t = np.arange(0, length, 1 / fs)
    wn = 2 * np.pi * fn
    wd = wn * np.sqrt(1 - zeta ** 2)
    return np.exp(-zeta * wn * t) * np.sin(wd * t)


def impact_train(f_defect, t, fs, rng, slip=SLIP, modulation=None):
    """Impact train at f_defect with random slip; optional amplitude modulation."""
    x = np.zeros_like(t)
    T = 1.0 / f_defect
    tk = rng.uniform(0, T)
    while tk < t[-1]:
        idx = int(round(tk * fs))
        amp = 1.0 if modulation is None else modulation(tk)
        x[idx] += amp
        tk += T * (1 + slip * rng.standard_normal())
    return np.convolve(x, impulse_response(fs))[: len(t)]


def generate(rpm=SHAFT_RPM, fs=FS, duration=DURATION, seed=SEED):
    rng = np.random.default_rng(seed)
    f = defect_frequencies(rpm=rpm)
    t = np.arange(0, duration, 1 / fs)
    fr, ftf = f["fr"], f["FTF"]

    # Baseline: shaft harmonics (imbalance/misalignment) + noise
    shaft = (0.20 * np.sin(2 * np.pi * fr * t)
             + 0.08 * np.sin(2 * np.pi * 2 * fr * t + 0.5)
             + 0.03 * np.sin(2 * np.pi * 3 * fr * t + 1.1))
    healthy = shaft + NOISE_STD * rng.standard_normal(len(t))

    # Outer race: fixed in load zone -> constant amplitude
    outer = 1.0 * impact_train(f["BPFO"], t, fs, rng)

    # Inner race: passes through load zone once per rev -> modulated at fr
    inner = 0.8 * impact_train(
        f["BPFI"], t, fs, rng,
        modulation=lambda tk: 0.5 * (1 + np.cos(2 * np.pi * fr * tk)) + 0.1)

    # Ball: defect hits inner & outer race each spin (2*BSF), carried by cage (FTF)
    ball = 0.6 * impact_train(
        2 * f["BSF"], t, fs, rng,
        modulation=lambda tk: 0.5 * (1 + np.cos(2 * np.pi * ftf * tk)) + 0.1)

    # Cage: low-frequency FTF content + FTF-modulated impacts
    # (a broken/worn cage pocket knocks once per cage revolution)
    cage = (0.10 * np.sin(2 * np.pi * ftf * t)
            + 0.05 * np.sin(2 * np.pi * 2 * ftf * t + 0.3)
            + 0.4 * impact_train(ftf, t, fs, rng))

    signals = {
        "healthy":  healthy,
        "outer":    healthy + outer,
        "inner":    healthy + inner,
        "ball":     healthy + ball,
        "cage":     healthy + cage,
        "combined": healthy + outer + inner + ball + cage,
    }
    return t, signals, f


# -------------------------------------------------------------------------
# 3. ENVELOPE SPECTRUM (to verify)
# -------------------------------------------------------------------------
def envelope_spectrum(x, fs, band=(2500, 4500)):
    from scipy.signal import butter, filtfilt, hilbert
    b, a = butter(4, [band[0] / (fs / 2), band[1] / (fs / 2)], btype="band")
    env = np.abs(hilbert(filtfilt(b, a, x)))
    env -= env.mean()
    spec = np.abs(np.fft.rfft(env)) / len(env) * 2
    freq = np.fft.rfftfreq(len(env), 1 / fs)
    return freq, spec


if __name__ == "__main__":
    import pandas as pd
    import matplotlib.pyplot as plt

    t, sig, f = generate()

    print("Bearing:", BEARING["designation"])
    print(f"Shaft speed: {SHAFT_RPM:.0f} rpm  (fr = {f['fr']:.3f} Hz)")
    for k in ("FTF", "BPFO", "BPFI", "BSF"):
        print(f"  {k:5s} = {f[k]:8.3f} Hz  ({f[k] / f['fr']:.4f} x fr)")

    df = pd.DataFrame({"time_s": t, **{f"{k}_g": v for k, v in sig.items()}})
    df.to_csv("bearing_signals.csv", index=False, float_format="%.6f")
    print("Saved bearing_signals.csv", df.shape)

    # Verification plot
    marks = {"outer": ("BPFO", f["BPFO"]), "inner": ("BPFI", f["BPFI"]),
             "ball": ("2xBSF", 2 * f["BSF"]), "cage": ("FTF", f["FTF"])}
    fig, ax = plt.subplots(4, 2, figsize=(13, 11))
    for i, (k, (lab, fd)) in enumerate(marks.items()):
        n = int(0.2 * FS)
        ax[i, 0].plot(t[:n], sig[k][:n], lw=0.6)
        ax[i, 0].set_title(f"{k} race / element fault - time signal")
        ax[i, 0].set_ylabel("Acc (g)")
        if k == "cage":
            fr_, X = np.fft.rfftfreq(len(t), 1 / FS), np.abs(np.fft.rfft(sig[k])) / len(t) * 2
            ax[i, 1].plot(fr_, X, lw=0.8)
            ax[i, 1].set_xlim(0, 100)
            ax[i, 1].set_title("cage fault - FFT spectrum")
        else:
            fe, E = envelope_spectrum(sig[k], FS)
            ax[i, 1].plot(fe, E, lw=0.8)
            ax[i, 1].set_xlim(0, 500)
            ax[i, 1].set_title(f"{k} fault - envelope spectrum")
        for h in (1, 2, 3):
            ax[i, 1].axvline(h * fd, color="r", ls="--", lw=0.8)
        ax[i, 1].text(fd, ax[i, 1].get_ylim()[1] * 0.9, f" {lab}={fd:.1f} Hz", color="r")
    ax[-1, 0].set_xlabel("Time (s)")
    ax[-1, 1].set_xlabel("Frequency (Hz)")
    plt.tight_layout()
    plt.savefig("bearing_signals_check.png", dpi=130)
    print("Saved bearing_signals_check.png")
