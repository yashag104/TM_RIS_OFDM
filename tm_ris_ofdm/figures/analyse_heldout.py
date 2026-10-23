"""Review E4: net throughput after pilot overhead with (N, Np[, hmax]) selected on the selection channels and
reported on disjoint held-out channels (week11_overhead_heldout*.json). Main metric: common-phase tracking.

Selection rule (pre-registered): for each receiver, channel model, speed and SNR, pick the configuration that
maximises the MEAN net rate over the 20 selection channels; report the held-out mean with a 95% t-interval and
paired held-out differences (ours - best CE-BEM, ours - one-tap). Also reported: the in-sample optimum on the
held-out channels themselves, to show the selection bias of the earlier table.
Writes results/week11_overhead_heldout_summary.json and paper/journal/tab_overhead_heldout.tex.
"""
import json
import os

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
NS_REFRESH = [4, 6, 8, 10, 15, 20, 30, 40, 60, 80, 120, 160, 240, 320]
CHANS = [("phys30", "Shared-path, 30\\,ns"), ("TDLA100", "TDL-A, 100\\,ns"),
         ("TDLC100", "TDL-C, 100\\,ns"), ("TDLD100", "TDL-D, 100\\,ns")]
METRIC = 1  # 0 sub (genie), 1 cpe (main), 2 none


def load():
    d = json.load(open(os.path.join(RES, "week11_overhead_heldout.json")))
    for ch, _ in CHANS:
        p = os.path.join(RES, f"week11_overhead_heldout_{ch}.json")
        if ch not in d and os.path.exists(p):
            d[ch] = json.load(open(p))[ch]
    return d


def net(curve, grid, Np, N):
    """Net rate (N - Np)/N * mean rate over data symbols Np..N-1; curve on the grid (per channel)."""
    g = np.asarray(grid, float); ok = (g >= Np) & np.isfinite(curve)
    if N <= Np or ok.sum() < 1 or (N - 1) > g[ok][-1]:
        return np.nan
    n = np.arange(Np, N)
    return (N - Np) / N * np.interp(n, g[ok], curve[ok]).mean()


def configs(arr, nps, iv, i_s, bem):
    """Yield (label, per-channel curves [nch x grid]) for every (Np[, hmax]) of one receiver."""
    if bem:
        for ip, np_ in enumerate(nps):
            for ih in range(arr.shape[2]):
                yield (np_, ih), arr[:, ip, ih, iv, i_s, :, METRIC]
    else:
        for ip, np_ in enumerate(nps):
            yield (np_, None), arr[:, ip, iv, i_s, :, METRIC]


def per_channel_net(curves, grid, Np, N):
    return np.array([net(c, grid, Np, N) for c in curves])


def select(arr_sel, nps, iv, i_s, grid, bem):
    best = (-np.inf, None)
    for (np_, ih), cur in configs(arr_sel, nps, iv, i_s, bem):
        for N in NS_REFRESH:
            v = per_channel_net(cur, grid, np_, N)
            if np.all(np.isfinite(v)) and v.mean() > best[0]:
                best = (v.mean(), (np_, ih, N))
    return best[1]


def heldout_net(arr, nps, iv, i_s, grid, bem, cfg):
    np_, ih, N = cfg
    ip = list(nps).index(np_)
    cur = arr[:, ip, ih, iv, i_s, :, METRIC] if bem else arr[:, ip, iv, i_s, :, METRIC]
    return per_channel_net(cur, grid, np_, N)


def tci(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return float(x.mean()), float(stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x)))


def main():
    d = load(); grid = d["grid"]; out = {}
    rows = []
    for ch, name in CHANS:
        if ch not in d:
            continue
        S = {k: np.array(v, float) for k, v in d[ch]["selection"].items()}
        H = {k: np.array(v, float) for k, v in d[ch]["heldout"].items()}
        out[ch] = {}
        for iv, v in enumerate(d["speeds_kmh"]):
            for i_s, snr in enumerate(d["snr_db"]):
                r = {}
                for rec, nps, bem in (("ours", d["Np_ours"], False), ("onetap", d["Np_ours"], False),
                                      ("bem", d["Np_bem"], True)):
                    cfg = select(S[rec], nps, iv, i_s, grid, bem)
                    vh = heldout_net(H[rec], nps, iv, i_s, grid, bem, cfg)
                    insample = select(H[rec], nps, iv, i_s, grid, bem)
                    vi = heldout_net(H[rec], nps, iv, i_s, grid, bem, insample)
                    r[rec] = {"cfg": [cfg[0], None if cfg[1] is None else d["hmax_bem"][cfg[1]], cfg[2]],
                              "heldout": tci(vh), "insample_on_heldout": float(np.mean(vi)), "_v": vh}
                pf = H["perfect"][:, iv, i_s, :]
                r["perfect"] = tci(np.nanmean(pf, axis=1))
                r["ours_minus_bem"] = tci(r["ours"]["_v"] - r["bem"]["_v"])
                r["ours_minus_onetap"] = tci(r["ours"]["_v"] - r["onetap"]["_v"])
                for k in ("ours", "onetap", "bem"):
                    del r[k]["_v"]
                out[ch][f"{v}kmh_{snr}dB"] = r
                o, b = r["ours"], r["bem"]
                rows.append(f"{name} & {v} & {snr} & ${o['heldout'][0]:.2f}\\pm{o['heldout'][1]:.2f}$ ({o['cfg'][0]}/{o['cfg'][2]}) & "
                            f"${b['heldout'][0]:.2f}\\pm{b['heldout'][1]:.2f}$ ({b['cfg'][0]}/{b['cfg'][2]}) & "
                            f"${r['ours_minus_bem'][0]:.2f}\\pm{r['ours_minus_bem'][1]:.2f}$ & "
                            f"{r['onetap']['heldout'][0]:.2f} & {r['perfect'][0]:.2f}\\\\")
        rows.append("\\midrule")
    rows = rows[:-1]
    json.dump(out, open(os.path.join(RES, "week11_overhead_heldout_summary.json"), "w"), indent=1)
    tex = ("\\begin{table}[!htbp]\\centering\n\\caption{Net throughput (bit/s/Hz) after pilot overhead, pre-registered "
           "selection/held-out protocol: pilot count, refresh interval (and CE-BEM harmonic range) chosen on 20 selection "
           "channels, reported on 20 disjoint held-out channels (mean $\\pm$ 95\\% $t$-interval); in brackets: pilot symbols / "
           "block length. Common-phase tracking; robust MMSE for every estimated receiver; CE-BEM up to all harmonics.}"
           "\\label{tab:overhead}\n\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}{llccccccc}\\toprule\n"
           "Channel & km/h & SNR & ours & best CE-BEM & ours $-$ CE-BEM (paired) & one-tap & perfect\\\\\\midrule\n"
           + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}}\n\\end{table}\n")
    open(os.path.join(HERE, "..", "paper", "journal", "tab_overhead_heldout.tex"), "w").write(tex)
    for ch in out:
        for k, r in out[ch].items():
            print(f"{ch:8s} {k:12s} ours {r['ours']['heldout'][0]:.2f}+-{r['ours']['heldout'][1]:.2f} cfg {r['ours']['cfg']} "
                  f"(in-sample {r['ours']['insample_on_heldout']:.2f}) | bem {r['bem']['heldout'][0]:.2f} cfg {r['bem']['cfg']} "
                  f"| diff {r['ours_minus_bem'][0]:.2f}+-{r['ours_minus_bem'][1]:.2f} | 1tap {r['onetap']['heldout'][0]:.2f} "
                  f"| perf {r['perfect'][0]:.2f}")


if __name__ == "__main__":
    main()
