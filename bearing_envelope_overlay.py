"""
Bearing Envelope-Spectrum Overlay  (Healthy vs Faulty)
------------------------------------------------------------
Follows the MathWorks "Rolling Element Bearing Fault Diagnosis" workflow
on your own recordings:

  1. band-pass the vibration around a resonance band (default 600-900 Hz)
  2. take the envelope (Hilbert) and its spectrum  -> the "envelope spectrum"
  3. a bearing defect shows up as peaks at the defect repetition rate and
     its harmonics, which are ABSENT in the healthy bearing
  4. scan for the strongest harmonic comb that is present in the faulty
     recordings but not in the healthy ones (no geometry needed), then
  5. overlay the bearing frequencies on the spectra.

Bearing frequencies for shaft speed fr, n balls, ball diameter d, pitch
diameter D, contact angle phi   (c = (d/D) cos(phi)):
    FTF  (cage)        = fr/2 * (1 - c)
    BPFO (outer race)  = n * FTF
    BPFI (inner race)  = n * (fr - FTF)
    BSF  (ball spin)   = fr * D/(2d) * (1 - c^2)         (ball defect appears at 2 x BSF)
Because FTF depends only on c, the scanned FTF gives you d/D directly
(for phi = 0:  d/D = 1 - 2*FTF/fr).  Once you know n, BPFO and BPFI follow.

Requirements:  pip install numpy scipy pandas matplotlib openpyxl

Usage:
  python bearing_envelope_overlay.py
  python bearing_envelope_overlay.py --band 500 700
  python bearing_envelope_overlay.py --n 8                 # also mark BPFO / BPFI for n = 8 balls
  python bearing_envelope_overlay.py --healthy H1.xlsx H2.xlsx --faulty F1.xlsx F2.xlsx --fr 24.8
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt, hilbert, welch, find_peaks


def load_stream(path):
    df = pd.read_excel(path)
    tc = df.iloc[:, 0].to_numpy(float)
    fs = int(round(1.0 / np.mean(np.diff(tc))))
    return df.iloc[:, 1:].to_numpy(float).T.reshape(-1), fs     # groups joined in order


def envelope_spectrum(x, fs, band, nperseg_s=20.0):
    sos = butter(4, [band[0] / (fs / 2), band[1] / (fs / 2)], btype="band", output="sos")
    y = sosfiltfilt(sos, x - x.mean())
    env = np.abs(hilbert(y))
    env -= env.mean()
    nper = int(nperseg_s * fs)
    f, p = welch(env, fs=fs, nperseg=nper, window="hann", noverlap=nper // 2, scaling="spectrum")
    return f, np.sqrt(2 * p)


def estimate_shaft_speed(x, fs):
    nper = 20 * fs
    f, p = welch(x - x.mean(), fs=fs, nperseg=nper, window="hann", noverlap=nper // 2, scaling="spectrum")
    amp = np.sqrt(2 * p)
    m = (f > 15) & (f < 40)
    cand = f[m][np.argsort(amp[m])[::-1][:3]]
    # the true shaft speed is the candidate whose 2x and 3x are also strong
    best, score = cand[0], -1
    for c in cand:
        s = sum(amp[np.argmin(abs(f - k * c))] for k in (1, 2, 3))
        if s > score:
            best, score = c, s
    return float(best)


def comb_score(f, num, den, fc, K=6, tol=0.12, fmax=200):
    vals = []
    for k in range(1, K + 1):
        fq = k * fc
        if fq > fmax:
            break
        w = (f > fq - tol) & (f < fq + tol)
        vals.append(np.log(num[w].max() / max(den[w].max(), 1e-12)))
    return float(np.mean(vals))


def main():
    ap = argparse.ArgumentParser(description="Envelope-spectrum overlay of bearing fault frequencies.")
    ap.add_argument("--healthy", nargs="+", default=["Vibration_Dataset_Healthy_2000Hz.xlsx",
                                                      "Vibration_Dataset_Healthy2_2000Hz.xlsx"])
    ap.add_argument("--faulty", nargs="+", default=["Vibration_Dataset_Faulty1_2000Hz.xlsx",
                                                     "Vibration_Dataset_Faulty2_2000Hz.xlsx"])
    ap.add_argument("--band", nargs=2, type=float, default=[600, 900], help="resonance band (Hz)")
    ap.add_argument("--fr", type=float, default=None, help="shaft speed in Hz (default: estimated)")
    ap.add_argument("--phi", type=float, default=0.0, help="contact angle in degrees (default 0)")
    ap.add_argument("--n", type=int, default=None, help="number of balls (adds BPFO / BPFI lines)")
    ap.add_argument("--xmax", type=float, default=130.0, help="max frequency of the plot (Hz)")
    ap.add_argument("--save", default="bearing_envelope_overlay.png")
    args = ap.parse_args()

    print("Loading recordings ...")
    H = [load_stream(p) for p in args.healthy]
    F = [load_stream(p) for p in args.faulty]
    fs = H[0][1]
    band = tuple(args.band)

    fr = args.fr or estimate_shaft_speed(H[0][0], fs)
    print(f"Shaft speed fr = {fr:.2f} Hz  ({fr * 60:.0f} rpm){'' if args.fr else '  [estimated]'}")

    Hs = [envelope_spectrum(x, fs, band)[1] for x, _ in H]
    Fs_ = [envelope_spectrum(x, fs, band)[1] for x, _ in F]
    f = envelope_spectrum(H[0][0], fs, band)[0]
    H_ref = np.maximum.reduce(Hs)                      # conservative healthy reference
    H_mean = np.mean(Hs, axis=0)

    # ---- scan for the strongest periodic comb that the faulty files have and the healthy ones don't ----
    grid = np.arange(3.0, 20.0, 0.01)
    null = max(max(comb_score(f, Hs[0], Hs[1], g) for g in grid),
               max(comb_score(f, Hs[1], Hs[0], g) for g in grid)) if len(Hs) > 1 else float("nan")
    scores = [np.array([comb_score(f, fs_, H_mean, g) for g in grid]) for fs_ in Fs_]
    best = [grid[np.argmax(s)] for s in scores]
    ftf_scan = float(np.mean(best))
    # refine: locate the real peaks near k x FTF in the mean faulty spectrum and fit f_k = k * FTF
    F_mean = np.mean(Fs_, axis=0)
    ks_fit, f_fit = [], []
    for k in range(1, 8):
        w = (f > k * ftf_scan - 0.25) & (f < k * ftf_scan + 0.25)
        if w.any():
            ks_fit.append(k)
            f_fit.append(f[w][np.argmax(F_mean[w])])
    ks_fit, f_fit = np.array(ks_fit), np.array(f_fit)
    ftf = float(np.sum(ks_fit * f_fit) / np.sum(ks_fit ** 2))
    print(f"\nHarmonic-comb scan (band {band[0]:.0f}-{band[1]:.0f} Hz):")
    for p, b, s in zip(args.faulty, best, scores):
        print(f"   {p}: best periodic frequency = {b:.2f} Hz  (score {s.max():.2f})")
    print(f"   Healthy-vs-Healthy best score (null level): {null:.2f}  -> a faulty score far above this is significant")
    print(f"   Refined by fitting the measured harmonic positions (k = 1..7): FTF = {ftf:.3f} Hz")

    c = 1.0 - 2.0 * ftf / fr                            # c = (d/D) cos(phi)
    d_over_D = c / np.cos(np.radians(args.phi))
    print(f"\nTreating {ftf:.2f} Hz as the cage frequency (FTF):  FTF/fr = {ftf / fr:.4f}")
    print(f"   => (d/D)cos(phi) = {c:.4f}   ->  d/D = {d_over_D:.4f}  (phi = {args.phi:g} deg)")
    bsf = fr * (1 / (2 * d_over_D)) * (1 - c ** 2)
    lines = {"FTF": ftf, "BSF": bsf, "2xBSF": 2 * bsf}
    print(f"   BSF = {bsf:.2f} Hz,  2xBSF = {2 * bsf:.2f} Hz")
    if args.n:
        lines["BPFO"] = args.n * ftf
        lines["BPFI"] = args.n * (fr - ftf)
        print(f"   n = {args.n}:  BPFO = {lines['BPFO']:.2f} Hz,  BPFI = {lines['BPFI']:.2f} Hz")
    else:
        print("   BPFO = n x FTF and BPFI = n x (fr - FTF); give --n once you have counted the balls.")
        for n in range(6, 11):
            print(f"      n = {n:2d}:  BPFO {n * ftf:6.1f} Hz   BPFI {n * (fr - ftf):6.1f} Hz")

    # ---- figure ----
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(13, 9), gridspec_kw={"height_ratios": [2.2, 1]})
    m = f <= args.xmax
    for i, h in enumerate(Hs):
        ax.plot(f[m], h[m], color="tab:blue", lw=1.0, alpha=0.85, label="Healthy" if i == 0 else None)
    for i, fa in enumerate(Fs_):
        ax.plot(f[m], fa[m], color="tab:red", lw=1.0, alpha=0.85, label="Faulty" if i == 0 else None)
    ymax = max(fa[m].max() for fa in Fs_) * 1.25
    ax.set_ylim(0, ymax)
    for k in range(1, int(args.xmax // ftf) + 1):
        ax.axvline(k * ftf, color="green", ls="--", lw=0.8, alpha=0.7)
        ax.text(k * ftf, ymax * 0.98, f"{k}", color="green", ha="center", va="top", fontsize=8)
    for k in range(1, int(args.xmax // fr) + 1):
        ax.axvline(k * fr, color="gray", ls=":", lw=0.9)
        ax.text(k * fr, ymax * 0.90, f"{k}X", color="dimgray", ha="center", va="top", fontsize=8)
    for name, val in lines.items():
        if name != "FTF" and val < args.xmax:
            ax.axvline(val, color="purple", ls="-.", lw=1.0)
            ax.text(val, ymax * 0.80, name, color="purple", rotation=90, ha="right", va="top", fontsize=8)
    ax.plot([], [], color="green", ls="--", label=f"harmonics of FTF = {ftf:.2f} Hz (numbered)")
    ax.plot([], [], color="gray", ls=":", label=f"shaft orders (fr = {fr:.1f} Hz)")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Envelope amplitude (g)")
    ax.set_title(f"Envelope spectrum, band {band[0]:.0f}-{band[1]:.0f} Hz - Healthy vs Faulty bearing")
    ax.set_xlim(0, args.xmax)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)

    ks = np.arange(1, 12)
    ratios = []
    for k in ks:
        tol = 0.10 + 0.015 * k                       # position error grows with k
        w = (f > k * ftf - tol) & (f < k * ftf + tol)
        ratios.append(min(fa[w].max() for fa in Fs_) / max(H_ref[w].max(), 1e-12))
    ax2.bar(ks, ratios, color="tab:green", edgecolor="black", linewidth=0.5)
    ax2.axhline(1, color="black", lw=1)
    ax2.set_xticks(ks)
    ax2.set_xlabel("Harmonic number k  (frequency = k x FTF)")
    ax2.set_ylabel("Faulty / Healthy\n(weaker faulty file vs stronger healthy file)")
    ax2.set_title("Faulty / Healthy amplitude at each harmonic of the fitted frequency (1 = no difference)")
    ax2.grid(True, alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(args.save, dpi=140)
    print(f"\nSaved {args.save}")
    plt.show()


if __name__ == "__main__":
    main()
