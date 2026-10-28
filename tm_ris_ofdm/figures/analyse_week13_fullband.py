"""Step 13b (week13_fullband): pool the channel files, then C13-3 (the 'dispersive' claim) and C13-4 (full-band pilots).
Run:  .venv/bin/python -I tm_ris_ofdm/figures/analyse_week13_fullband.py
"""
import json
import pathlib

import numpy as np

RES = pathlib.Path(__file__).resolve().parent.parent / "results"
files = ["week13_fullband"] + [f"week13_fullband_part_{p}" for p in "abc"]
ser, nm, tr = [], [], []
for f in files:
    d = json.loads((RES / f"{f}.json").read_text())
    s, n = np.array(d["ser"], float), np.array(d["nmse"], float)
    t = np.atleast_1d(d.get("trials", np.arange(1, s.shape[0] + 1)))[: s.shape[0]]
    keep = t <= 8 if f == "week13_fullband" else np.ones(len(t), bool)
    ser.append(s[keep]); nm.append(n[keep]); tr += list(t[keep])
names, Nps = d["names"], d["Np"]
cfg = [(20, 16), (28, 64)]                                        # as in week13_fullband.m (cell array in the JSON)
ser, nm, tr = np.concatenate(ser), np.concatenate(nm), np.array(tr)
assert sorted(tr) == list(range(1, 31)), sorted(tr)
print(f"{len(tr)} channels; receivers {names}")
rng = np.random.default_rng(0)


def ratio_ci(a, b, B=20000):
    """mean(a)/mean(b) with a paired (channel) bootstrap 95% interval"""
    idx = rng.integers(0, len(a), (B, len(a)))
    r = a[idx].mean(1) / np.maximum(b[idx].mean(1), 1e-12)
    return a.mean() / b.mean(), np.percentile(r, 2.5), np.percentile(r, 97.5)


for ic, (snr, mq) in enumerate(cfg):
    for ip, npil in enumerate(Nps):
        print(f"\n{snr} dB {mq}-QAM Np={npil}:  SER  | NMSE dB")
        for e, k in enumerate(names):
            nmd = "" if e == 0 else f"{10*np.log10(np.nanmean(nm[:, e, ip, ic])):6.1f}"
            print(f"  {k:8s} {ser[:, e, ip, ic].mean():.5f}  {nmd}")
        io = names.index("ours")
        for k in ("flatcal", "flatest", "genie", "ideal", "bem"):
            m, lo, hi = ratio_ci(ser[:, names.index(k), ip, ic], ser[:, io, ip, ic])
            tag = ""
            if snr == 28 and npil == 2 and k in ("flatcal", "flatest"):
                tag = "  C13-3 " + ("PASS" if lo > 1 else "FAIL")
            print(f"  SER {k}/ours = {m:.2f} [{lo:.2f}, {hi:.2f}]{tag}")
        m, lo, hi = ratio_ci(ser[:, io, ip, ic], ser[:, names.index("true"), ip, ic])
        print(f"  SER ours/true = {m:.2f} [{lo:.2f}, {hi:.2f}]" + ("  (C13-4 candidate)" if snr == 20 and m <= 1.2 else ""))
