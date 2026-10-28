"""Step 13a (week13_review_baselines): pre-registered checks C13-1, C13-2, C13-S and the reported quantities.
Run:  .venv/bin/python -I tm_ris_ofdm/figures/analyse_week13.py
"""
import json
import pathlib

import numpy as np

RES = pathlib.Path(__file__).resolve().parent.parent / "results"


def ci(x):
    x = np.asarray(x, float); n = len(x)
    from scipy import stats
    h = stats.t.ppf(0.975, n - 1) * x.std(ddof=1) / np.sqrt(n)
    return x.mean(), h


def load(c):
    return json.loads((RES / f"week13_review_baselines_{c}.json").read_text())[c]


PRACT = ["linear", "flatcal", "flatest", "static", "omp_T", "ideal", "bemall_Np32"]
VARS = ["ours_tw05", "ours_tw2", "ours_rlo", "ours_rhi"]
for c in ("phys30", "TDLC100"):
    s = load(c); e = s["estimators"]; snr = s["snr_db"]
    rt, rp, nm = np.array(s["rate"], float), np.array(s["rate_perfect"], float), np.array(s["nmse"], float)
    loss = rp[:, None, :] - rt                                       # channel x est x snr
    n = loss.shape[0]
    print(f"\n=== {c}: {n} channels, Np={s['Np']}, L={s['L']}, r={s['r']}, unknowns {s['n_unknowns_lowrank']}, obs {s['n_obs']}")
    print(f"{'estimator':12s} " + " ".join(f"loss@{x:2d}dB      " for x in snr) + " " + " ".join(f"NMSE@{x}" for x in snr))
    for i, k in enumerate(e):
        cells = " ".join(f"{m:6.3f}+-{h:5.3f}" for m, h in (ci(loss[:, i, j]) for j in range(len(snr))))
        nmd = " ".join(f"{10*np.log10(nm[:, i, j].mean()):6.1f}" for j in range(len(snr)))
        print(f"{k:12s} {cells}  {nmd}")
    dof = np.array(s["dof_linear"], float); fb = np.array(s["fallback_used"], float)
    print("effective dof of linear stage (mean, min-max):", [f"{dof[:, j].mean():.1f} ({dof[:, j].min():.1f}-{dof[:, j].max():.1f})" for j in range(len(snr))])
    print("fallback used (fraction):", np.round(fb.mean(0), 2))
    io = e.index("ours_nofb")
    for j, x in enumerate(snr):
        if x not in (20, 30):
            continue
        best = min(PRACT, key=lambda k: loss[:, e.index(k), j].mean())
        m_o = loss[:, io, j].mean(); m_b = loss[:, e.index(best), j].mean()
        print(f"C13-2 {x} dB: ours_nofb {m_o:.3f} vs best practical '{best}' {m_b:.3f} -> {'PASS' if m_o <= m_b + 0.02 else 'FAIL'}")
        for v in VARS:
            d = loss[:, e.index(v), j] - loss[:, io, j]; m, h = ci(d)
            tol = 0.10 if x == 20 else 0.20
            print(f"  C13-S {v:10s} minus nominal {m:+.3f}+-{h:.3f} -> {'ok' if abs(m) <= tol else 'OUT'}")
        d = loss[:, io, j] - loss[:, e.index("genie_lr"), j]; m, h = ci(d)
        print(f"  ours_nofb - genie_lr {m:+.3f}+-{h:.3f};", end=" ")
        d = loss[:, e.index("ours"), j] - loss[:, io, j]; m, h = ci(d)
        print(f"ours(fallback) - ours_nofb {m:+.3f}+-{h:.3f}")
        gi, gg, gf = (loss[:, e.index(k), j].mean() for k in ("ideal", "genie", "flatcal"))
        print(f"  flatcal closes {(gi - gf) / (gi - gg) * 100:.0f}% of the ideal-to-genie gap ({gi:.3f} -> {gf:.3f}; genie {gg:.3f})")
        if c == "phys30" and x == 30:
            d = loss[:, e.index("flatcal"), j] - loss[:, io, j]; m, h = ci(d)
            print(f"C13-1 flatcal - ours_nofb at 30 dB: {m:+.3f}+-{h:.3f} -> {'PASS' if m - h > 0 else 'FAIL'}")
