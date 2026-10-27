"""Steps 11-12 (week12_snr_tang): print NMSE, rate loss and BER (QPSK, 16-QAM) vs SNR 10-40 dB for every estimator.
Run:  .venv/bin/python -I tm_ris_ofdm/figures/analyse_snr_tang.py   (numbers quoted with Fig. snr in main.tex)
"""
import json
import pathlib

import numpy as np

RES = pathlib.Path(__file__).resolve().parent.parent / "results"
R = json.loads((RES / "week12_snr_tang.json").read_text())                    # Step 11: 10:5:40 dB
for c in ("phys30", "TDLC100"):                                               # Step 12: 10:2:30 dB (plotted)
    f = RES / f"week12_snr_tang_fine_{c}.json"
    if f.exists():
        R[c + " (fine grid)"] = json.loads(f.read_text())[c]
for c, s in R.items():
    e = s["estimators"]
    nm, rt, rp = (np.array(s[k], float) for k in ("nmse", "rate", "rate_perfect"))
    be, bp = np.array(s["ber"], float), np.array(s["ber_perfect"], float)
    print(f"\n{c}  Np={s['Np']}  SNR={s['snr_db']}  bits/channel={s['bits_per_channel']}")
    for i, k in enumerate(e):
        print(f"{k:12s} NMSE dB", np.round(10 * np.log10(np.nanmean(nm[:, i], 0)), 1))
        print(f"{'':12s} loss   ", np.round(np.nanmean(rp - rt[:, i], 0), 3))
        for q, Mq in enumerate(s["qam"]):
            print(f"{'':12s} BER {Mq:2d} ", " ".join(f"{x:.1e}" for x in np.nanmean(be[:, i, :, q], 0)))
    for q, Mq in enumerate(s["qam"]):
        print(f"{'perfect':12s} BER {Mq:2d} ", " ".join(f"{x:.1e}" for x in np.nanmean(bp[:, :, q], 0)))
