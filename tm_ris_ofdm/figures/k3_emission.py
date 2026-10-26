"""K3: adjacent-channel leakage of the signal reflected by a switched (time-modulated) RIS element.

Model (results/k3_oob.json): ideal +-1 element, duty D = 0.75 in one block per OFDM symbol, K = 8 slots, uniform loading
of 3276 subcarriers (30 kHz), neighbouring channels at 100 and 200 MHz spacing, linear switching ramps of rise time
rho * (slot length). The reflected spectrum is sum_h |c_h|^2 S(f - h df) with the Fourier coefficients of the switching
waveform: ideal edges give c_0 = 2D - 1, c_h = 2 (1 - exp(-j 2 pi h D)) / (j 2 pi h); a linear ramp of duration
rho T/K multiplies c_h by sinc(h rho / K). ACLR = in-channel power / power in the adjacent channel (worse side).
Run:  .venv/bin/python -I tm_ris_ofdm/figures/k3_emission.py   (prints and compares with the saved JSON)
"""
import json
import os

import numpy as np

D, K, M, SPACING_HZ, DF = 0.75, 8, 3276, 100e6, 30e3
HMAX = 2_000_000                                  # harmonics kept; the 1/h^2 tail beyond is negligible


def coeffs(rho):
    h = np.arange(-HMAX, HMAX + 1)
    c = np.empty(h.size, complex)
    nz = h != 0
    c[nz] = 2 * (1 - np.exp(-2j * np.pi * h[nz] * D)) / (2j * np.pi * h[nz])
    c[~nz] = 2 * D - 1
    return h, np.abs(c * np.sinc(h * rho / K)) ** 2


def band_power(h, p, lo, hi):
    """Power landing in subcarrier interval [lo, hi) from a flat band m in [-M/2, M/2) shifted by each harmonic h."""
    a, b = -M // 2 + h, M // 2 + h                # shifted band [a, b)
    overlap = np.clip(np.minimum(b, hi) - np.maximum(a, lo), 0, None)
    return float(np.sum(p * overlap))


def main():
    off1, off2 = SPACING_HZ / DF, 2 * SPACING_HZ / DF
    out = {}
    for rho in (0.0, 0.001, 0.01, 0.05, 0.1, 0.3):
        h, p = coeffs(rho)
        main_ = band_power(h, p, -M // 2, M // 2)
        adj = [band_power(h, p, int(round(s * off1)) - M // 2, int(round(s * off1)) + M // 2) for s in (1, -1)]
        adj2 = [band_power(h, p, int(round(s * off2)) - M // 2, int(round(s * off2)) + M // 2) for s in (1, -1)]
        total = float(np.sum(p)) * M
        out[f"rho_{rho}"] = {"ACLR_dB": 10 * np.log10(main_ / max(adj)), "ACLR2_dB": 10 * np.log10(main_ / max(adj2)),
                             "outside_band_frac_dB": 10 * np.log10((total - main_) / total)}
    saved = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "k3_oob.json")
    ref = json.load(open(saved))["results"] if os.path.exists(saved) else {}
    print(f"{'rise time':>10} {'ACLR':>8} {'(saved)':>8} {'ACLR2':>8} {'(saved)':>8} {'outside':>8} {'(saved)':>8}")
    for k, v in out.items():
        r = ref.get(k, {})
        print(f"{k[4:]:>10} {v['ACLR_dB']:8.2f} {r.get('ACLR_dB', np.nan):8.2f} {v['ACLR2_dB']:8.2f} "
              f"{r.get('ACLR2_dB', np.nan):8.2f} {v['outside_band_frac_dB']:8.2f} {r.get('outside_band_frac_dB', np.nan):8.2f}")


if __name__ == "__main__":
    main()
