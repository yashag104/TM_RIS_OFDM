"""Print every number the paper quotes from the week-12 runs, with the pre-registered checks.

Reads results/week12_standard_curves.json (Step 8), week12_c6_lowsnr.json (Step 9) and week12_c6b.json (Step 10).
Missing files are skipped. Run:  .venv/bin/python -I tm_ris_ofdm/figures/analyse_week12.py
"""
import json
import os

import numpy as np
from scipy import stats

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")


def load(name):
    p = os.path.join(RES, name + ".json")
    return json.load(open(p)) if os.path.exists(p) else None


def A(x):
    return np.array(x, float)


def ci(x):
    """Mean and 95% t half-width of the finite entries."""
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return x.mean(), stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))


def db(x):
    return 10 * np.log10(np.nanmean(x))


def step8(d):
    print("=" * 30, "Step 8: week12_standard_curves (Figs. 3-5)")
    for ch in ("phys30", "TDLC100"):
        s = d["S1"][ch]; e = s["estimators"]; snr = s["snr_db"]
        nm, rt, be = A(s["nmse"]), A(s["rate"]), A(s["ber"])
        rp, bp = A(s["rate_perfect"]), A(s["ber_perfect"])
        print(f"\n-- S1 {ch}, Np = {s['Np']}, {nm.shape[0]} channels; columns: SNR {snr} dB")
        print("NMSE (dB)")
        for i, n in enumerate(e):
            print(f"  {n:12s}" + "".join(f"{db(nm[:, i, k]):7.1f}" for k in range(len(snr))))
        print("loss to perfect CSI (bit/s/Hz)")
        for i, n in enumerate(e):
            print(f"  {n:12s}" + "".join(f"{np.nanmean(rp[:, k] - rt[:, i, k]):7.3f}" for k in range(len(snr))))
        print("uncoded 16-QAM BER")
        print(f"  {'perfect':12s}" + "".join(f"{np.nanmean(bp[:, k]):9.1e}" for k in range(len(snr))))
        for i, n in enumerate(e):
            print(f"  {n:12s}" + "".join(f"{np.nanmean(be[:, i, k]):9.1e}" for k in range(len(snr))))
    s = d["S2"]; e = s["estimators"]; al = s["alpha"]; rt = A(s["rate"]); rp = A(s["rate_perfect"])
    print("\n-- S2 element non-ideality continuum, loss (mean +- 95% CI)")
    for k, sn in enumerate(s["snr_db"]):
        print(f"  {sn} dB, alpha = {al}")
        for i, n in enumerate(e):
            print(f"    {n:8s}" + "".join("%8.3f+-%.3f" % ci(rp[:, a, k] - rt[:, i, a, k]) for a in range(len(al))))
    g = {n: rp[:, :, 0] - rt[:, e.index(n), :, 0] for n in e}
    cross = [a for a, v in zip(al, g["ideal"].mean(0)) if v > 0.3]
    print(f"  check (b): ours - exact-element at alpha = 1, 20 dB: {np.mean(g['ours'][:, -1] - g['genie'][:, -1]):+.3f} (pass if < 0.1)")
    print(f"  check (c): smallest alpha with ideal-element loss > 0.3 at 20 dB: {cross[0] if cross else 'none'}")
    s = d["S3"]; e = s["estimators"]; nm = A(s["nmse"])
    print(f"\n-- S3 NMSE (dB) vs Np = {s['Np']} at {s['snr_db']} dB")
    for i, n in enumerate(e):
        print(f"  {n:8s}" + "".join(f"{db(nm[:, i, k]):7.1f}" for k in range(len(s["Np"]))))


def step9(d):
    print("\n" + "=" * 30, "Step 9: week12_c6_lowsnr")
    snr = d["P1"]["snr_db"]
    for ch, frozen in (("phys30", 8), ("TDLC100", 16)):
        c = d["P1"][ch]; ranks = c["ranks"]; ifr = ranks.index(frozen)
        gs = A(c["selection"]["rate_perfect"])[:, None, :] - A(c["selection"]["rate"])
        gh = A(c["heldout"]["rate_perfect"])[:, None, :] - A(c["heldout"]["rate"])
        pen = np.array([3 * c["L"] * r * np.log(d["P1"]["M"] * c["Np"]) for r in ranks])
        rB = np.argmin(A(c["heldout"]["resid_over_s2"]) + pen[None, :, None], axis=1)
        print(f"\n-- low-SNR rank, {ch}: ranks {ranks}, frozen {frozen}; improvement over frozen (held-out, paired)")
        for k, s in enumerate(snr):
            ia = int(np.argmin(gs[:, :, k].mean(0)))
            dA = ci(gh[:, ifr, k] - gh[:, ia, k])
            dB = ci(gh[:, ifr, k] - gh[np.arange(gh.shape[0]), rB[:, k], k])
            print(f"  {s:2d} dB: rule A r = {ranks[ia]:2d} {dA[0]:+.3f}+-{dA[1]:.3f} | rule B (BIC) {dB[0]:+.3f}+-{dB[1]:.3f}")
    for ch in ("phys30", "TDLC100"):
        c = d["P2"][ch]; n = c["names"]; rt = A(c["rate"]); rp = A(c["rate_perfect"])
        io, il = n.index("ours"), n.index(n[-2])
        print(f"\n-- C6 learned prior, {ch}, Np = {c['Np']}, components {c['ncomp']}; loss by SNR {c['snr_db']}")
        for i, nn in enumerate(n):
            print(f"  {nn:12s}" + "".join(f"{np.mean(rp[:, k] - rt[:, i, k]):7.3f}" for k in range(len(c["snr_db"]))))
        print("  ours - " + n[il] + " (paired): " +
              "  ".join("%+.3f+-%.3f" % ci(rt[:, il, k] - rt[:, io, k]) for k in range(len(c["snr_db"]))))


def step10(d):
    print("\n" + "=" * 30, "Step 10: week12_c6b (Table: learned baseline)")
    for tn in ("phys30", "phys100K3", "TDLC", "TDLA", "TDLD"):
        if tn not in d:
            continue
        c = d[tn]; n = c["names"]; rt = A(c["rate"]); rp = A(c["rate_perfect"]); io = n.index("ours")
        print(f"\n-- test {tn}, components {c['ncomp']}")
        for k, s in enumerate(c["snr_db"]):
            g = rp[:, None, k] - rt[:, :, k]
            print(f"  {s:2d} dB  " + "  ".join(f"{n[i]} {g[:, i].mean():.3f}" for i in range(len(n))) + "  | " +
                  "; ".join("%s - ours %+.3f+-%.3f" % (n[i], *ci(g[:, i] - g[:, io])) for i in range(len(n)) if i != io))


if __name__ == "__main__":
    for name, fn in (("week12_standard_curves", step8), ("week12_c6_lowsnr", step9), ("week12_c6b", step10)):
        d = load(name)
        if d is None:
            print(f"(missing results/{name}.json, skipped)")
        else:
            fn(d)
