"""Operation counts for the timing-only estimator (est_general), matrix-free form.

Per alternating-least-squares iteration, with L delay basis functions per hop,
rank r, G distinct switching patterns, Np pilot symbols, FFT length Nf = 2^ceil(log2 4M):
  regressors  : 2 states x G x Np x L Toeplitz products by FFT -> 2 G Np L (Nf log2 Nf + 2 Nf)
                (computed once, reused by every iteration)
  U step      : normal equations (M Np) x (L r) -> M Np (L r)^2 + (L r)^3 / 3
  V step      : (M Np) x (2 L r)                -> M Np (2 L r)^2 + (2 L r)^3 / 3
  Lam step    : (M Np) x (G r^2)                -> M Np (G r^2)^2 + (G r^2)^3 / 3
L grows with bandwidth: L = ceil(M df tau) + 10 (tau = 0.2 us at full band).
Compared against the dense alternative of estimating every band entry of H.
Writes results/week6_complexity.json.
"""

import json
from pathlib import Path

import numpy as np

df, tau, r, Np, iters = 30e3, 0.2e-6, 8, 2, 20
Ms = [256, 512, 1024, 2048, 3276]
out = {"M": Ms, "assumptions": {"df": df, "tau_s": tau, "rank": r, "Np": Np, "iterations": iters}}
for G in (1, 8, 16):
    ops = []
    for M in Ms:
        L = int(np.ceil(M * df * tau)) + 10
        nf = 2 ** int(np.ceil(np.log2(4 * M)))
        reg = 2 * G * Np * L * (nf * np.log2(nf) + 2 * nf)
        per_it = sum(M * Np * n ** 2 + n ** 3 / 3 for n in (L * r, 2 * L * r, G * r * r))
        ops.append(float(reg + iters * per_it))
    out[f"ops_G{G}"] = ops
# unstructured alternative: LS of every band entry with 2*hmax+1 unknowns per row
hmax = 16
out["ops_banded_unstructured"] = [float(M * ((2 * hmax + 1) ** 3 / 3 + (2 * hmax + 1) ** 2 * (2 * hmax + 1)))
                                  for M in Ms]
out["pilots_needed_banded_unstructured"] = 2 * hmax + 1
out["measured_matlab_s"] = {"M256_G1": 0.4, "M256_G16": 2.1, "M3276_G1_lowrank_with_DD": 27.0,
                            "note": "est_general / week4_fullband_trials on this laptop, MATLAB R2026a"}
Path(__file__).resolve().parent.parent.joinpath("results", "week6_complexity.json").write_text(json.dumps(out, indent=1))
for k in ("ops_G1", "ops_G8", "ops_G16"):
    print(k, ["%.2e" % v for v in out[k]])
