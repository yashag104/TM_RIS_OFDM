"""Paper figures for the TM-RIS receiver, drawn only from saved results/*.json.

No simulation runs here. Each figure names the experiment script that produced
its data. Output: PNG (screen) and PDF (paper) next to this file.

Run:  .venv/bin/python tm_ris_ofdm/figures/make_figures.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE.parent / "results"

# Validated categorical slots (reference palette, first three: all-pairs pass, light mode).
OURS, BEM, TOEP = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
ONETAP = "#8a8985"          # neutral: the conventional receiver is a reference, not a series hue

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 9.5,
    "axes.labelsize": 9, "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 200,
    "lines.linewidth": 2.0, "lines.markersize": 6, "figure.facecolor": "white",
})


def load(name):
    return json.loads((RES / f"{name}.json").read_text())


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------------------
# Fig. 1 - estimators vs pilot budget (week4_bem_compare.m)
# ---------------------------------------------------------------------------
def fig_estimators():
    d = load("week4_bem_compare")
    names, nps, snrs = d["names"], d["Np"], d["snr_db"]
    mean, ci = np.array(d["rate_mean"]), np.array(d["rate_ci95"])     # names x Np x snr
    style = {
        "ours": ("Ours (timing only)", OURS, "o", "-"),
        "bem12": ("CE-BEM, 25 harmonics", BEM, "s", "-"),
        "toeplitz": ("Toeplitz mixing model", TOEP, "^", "-"),
        "onetap": ("One-tap (conventional)", ONETAP, "D", "-"),
        "perfect": ("Perfect channel", INK, None, "--"),
    }
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8), sharey=False)
    for j, (ax, snr) in enumerate(zip(axes, snrs)):
        for key, (lab, col, mk, ls) in style.items():
            i = names.index(key)
            ax.errorbar(nps, mean[i, :, j], yerr=ci[i, :, j], color=col, marker=mk, ls=ls,
                        capsize=2.5, label=lab, lw=1.6 if key == "perfect" else 2.0,
                        mec="white", mew=0.8)
        ax.set_xscale("log", base=2)
        ax.set_xticks(nps, [str(n) for n in nps])
        ax.set_xlabel("Pilot OFDM symbols $N_p$")
        ax.set_title(f"SNR {snr} dB", loc="left", color=INK)
        ax.set_ylim(0, 10.5 if snr == 30 else 7.2)
    axes[0].set_ylabel("Rate after MMSE (bit/s/Hz)")
    axes[1].legend(loc="center right", frameon=False)
    fig.suptitle("Ours approaches the perfect-channel rate with 2 pilots; the literature's estimators stall",
                 x=0.01, ha="left", fontsize=10, color=INK)
    fig.tight_layout()
    save(fig, "fig1_estimators_vs_pilots")


# ---------------------------------------------------------------------------
# Fig. 2 - what the receiver must know about the element (week3_modelfree_confirm.m)
# ---------------------------------------------------------------------------
def fig_model_knowledge():
    d = load("week3_modelfree_confirm")
    snrs = d["snr_db"]
    loss = np.array(d["S1"]["loss_db"])                 # rows: exact, df0 10 MHz, ideal, timing-only
    slope = np.array([0.30, 0.31, 0.33])
    ci = np.array(d["S1"]["loss_bits_ci95"]) / slope
    rows = [("Exact element model", 0, "#9a9993"), ("Model off by 10 MHz", 1, BEM),
            ("Ideal ±1 switch (prior papers)", 2, TOEP), ("Ours: timing only", 3, OURS)]
    fig, ax = plt.subplots(figsize=(4.6, 2.9))
    w = 0.19
    x = np.arange(len(snrs))
    for k, (lab, r, col) in enumerate(rows):
        pos = x + (k - 1.5) * (w + 0.02)
        ax.bar(pos, loss[r], w, color=col, label=lab, yerr=ci[r], capsize=2,
               error_kw={"ecolor": MUTED, "lw": 0.8})
        for p, v in zip(pos, loss[r]):
            if r in (2, 3):
                ax.text(p, v + 0.15, f"{v:.2f}", ha="center", va="bottom", fontsize=7, color=INK)
    ax.set_xticks(x, [f"{s} dB" for s in snrs])
    ax.set_ylabel("Loss vs perfect channel (dB)")
    ax.set_xlabel("SNR")
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("No element model needed: timing-only matches the exact model",
                 loc="left", color=INK)
    fig.tight_layout()
    save(fig, "fig2_element_model_knowledge")


# ---------------------------------------------------------------------------
# Fig. 3 - full band, M = 3276 (week4_fullband_trials.m, week4_fullband_s2_shift.m)
# ---------------------------------------------------------------------------
def fig_fullband():
    d = load("week4_fullband_trials")
    s2 = load("week4_fullband_s2_shift")
    names = d["names"]
    cases = [("snr20_qam16", "20 dB, 16-QAM"), ("snr28_qam64", "28 dB, 64-QAM")]
    order = [("onetap", "One-tap", ONETAP), ("bem_dd", "CE-BEM", BEM),
             ("timing_dd", "Ours (timing only)", OURS), ("true", "Perfect channel", INK)]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8), gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    w = 0.19
    x = np.arange(len(cases))
    for k, (key, lab, col) in enumerate(order):
        i = names.index(key)
        m = [d[c]["ser_mean"][i] for c, _ in cases]
        e = [d[c]["ser_ci95"][i] for c, _ in cases]
        pos = x + (k - 1.5) * (w + 0.02)
        lo = [min(ei, mi * 0.9) for mi, ei in zip(m, e)]            # keep log-axis error bars positive
        ax.bar(pos, m, w, color=col, label=lab, yerr=[lo, e], capsize=2,
               error_kw={"ecolor": MUTED, "lw": 0.8})
    ax.set_yscale("log"); ax.set_ylim(5e-3, 1.5)
    ax.set_xticks(x, [lab for _, lab in cases])
    ax.set_ylabel("Symbol error rate")
    ax.set_title("Identical patterns, 10 channels", loc="left", color=INK)
    ax.legend(frameon=False, loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.14))
    # right: per-element patterns, shift-structure estimator
    ax = axes[1]
    order2 = [(0, "One-tap", ONETAP), (2, "Ours (timing only)", OURS), (1, "Perfect channel", INK)]
    for k, (i, lab, col) in enumerate(order2):
        m = [s2[c]["ser_mean"][i] for c, _ in cases]
        e = [s2[c]["ser_ci95"][i] for c, _ in cases]
        pos = x + (k - 1) * (0.25 + 0.02)
        lo = [min(ei, mi * 0.9) for mi, ei in zip(m, e)]
        ax.bar(pos, m, 0.25, color=col, yerr=[lo, e], capsize=2, error_kw={"ecolor": MUTED, "lw": 0.8})
    ax.set_yscale("log"); ax.set_ylim(5e-3, 1.5)
    ax.set_xticks(x, [lab for _, lab in cases])
    ax.set_title("Per-element patterns, 5 channels", loc="left", color=INK)
    for xi, (c, _) in zip(x, cases):
        ax.text(xi, 1.15, f"{s2[c]['op_nmse_db_mean']:.1f} dB NMSE", ha="center", fontsize=7, color=MUTED)
    fig.suptitle("Full n78 band (3276 subcarriers), 2 pilots + 1 refinement pass",
                 x=0.01, ha="left", fontsize=10, color=INK)
    fig.tight_layout()
    save(fig, "fig3_full_band")


# ---------------------------------------------------------------------------
# Fig. 4 - fair comparison: pattern designed to favour the one-tap receiver (week3_designed_user.m)
# ---------------------------------------------------------------------------
def fig_fair():
    d = load("week3_designed_user")
    snrs = d["snr_db"]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7), sharey=True)
    for ax, key, title in [(axes[0], "D875_designed", "Pattern designed for the user (duty 0.875)"),
                           (axes[1], "D750_generic", "Generic user position (duty 0.75)")]:
        r = d[key]
        ax.plot(snrs, r["rate_mmse_perfect"], color=INK, ls="--", lw=1.6, label="Perfect channel")
        ax.plot(snrs, r["rate_mmse_est_Np2"], color=OURS, marker="o", label="Ours, 2 pilots", mec="white")
        ax.plot(snrs, r["rate_tin"], color=ONETAP, marker="D", label="One-tap (conventional)", mec="white")
        ax.set_xticks(snrs); ax.set_xlabel("SNR (dB)")
        ax.set_title(f"{title}\none-tap SIR {r['sir_db']:+.1f} dB", loc="left", color=INK)
    axes[0].set_ylabel("Rate (bit/s/Hz)")
    axes[0].legend(frameon=False, loc="upper left")
    fig.tight_layout()
    save(fig, "fig4_fair_comparison")


# ---------------------------------------------------------------------------
# Fig. 5 - robustness to longer delay spread (week3_delayspread.m)
# ---------------------------------------------------------------------------
def fig_delay():
    d = load("week3_delayspread")
    keys = ["rms30ns", "rms100ns", "rms200ns"]
    labels = ["30 ns", "100 ns", "200 ns"]
    fig, ax = plt.subplots(figsize=(4.2, 2.6))
    x = np.arange(len(keys)); w = 0.32
    for j, (snr, col) in enumerate(zip(d["snr_db"], [OURS, "#1a4f94"])):
        m = [d[k]["loss_db_mean"][j] for k in keys]
        e = [d[k]["loss_db_ci95"][j] for k in keys]
        ax.bar(x + (j - 0.5) * (w + 0.03), m, w, color=col, yerr=e, capsize=2,
               error_kw={"ecolor": MUTED, "lw": 0.8}, label=f"{snr} dB")
    ax.axhline(1.0, color=MUTED, lw=0.9, ls=":")
    ax.text(-0.45, 1.03, "1 dB", fontsize=7, color=MUTED, ha="left", va="bottom")
    ax.set_xticks(x, labels)
    ax.set_xlabel("RMS delay spread per hop (estimator assumes ≤ 0.3 µs support)")
    ax.set_ylabel("Loss vs perfect channel (dB)")
    ax.legend(frameon=False, title="SNR", loc="upper left")
    ax.set_title("Longer delay spreads, estimator settings unchanged", loc="left", color=INK)
    fig.tight_layout()
    save(fig, "fig5_delay_spread")


# ---------------------------------------------------------------------------
# Fig. 0 - system model block diagram
# ---------------------------------------------------------------------------
def fig_system():
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    fig, ax = plt.subplots(figsize=(7.0, 2.7))
    ax.set_xlim(0, 100); ax.set_ylim(0, 40); ax.axis("off")
    soft = "#e8f0fb"

    def box(x, y, w, h, text, fc="white", ec=MUTED):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2",
                                    fc=fc, ec=ec, lw=1.0))
        ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=7.3, color=INK,
                linespacing=1.25)

    def arrow(p0, p1, col=INK, ls="-", rad=0.0):
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=9, lw=1.0, color=col,
                                     ls=ls, shrinkA=0, shrinkB=0, connectionstyle=f"arc3,rad={rad}"))

    box(1, 21, 12, 9, "OFDM Tx\nIFFT + CP")
    box(25, 19, 17, 13, "TM-RIS\n$Q$ elements\n$K$ slots / symbol", fc=soft, ec=OURS)
    box(25, 2, 17, 8, "Controller\n(switching timing)", fc=soft, ec=OURS)
    box(53, 21, 11, 9, "Remove CP\nFFT")
    box(68, 25, 17, 9, "Matrix-free\nCG-MMSE", fc=soft, ec=OURS)
    box(68, 11, 17, 9, "Timing-only\nestimator", fc=soft, ec=OURS)
    box(89, 25, 10, 9, "Detector")
    arrow((13.3, 25.5), (24.7, 25.5)); ax.text(19, 27.2, r"$\mathbf{A}^{\rm s}_q(f)$", ha="center", fontsize=8)
    arrow((42.3, 25.5), (52.7, 25.5)); ax.text(47.5, 27.2, r"$\mathbf{A}^{\rm d}_q(f)$", ha="center", fontsize=8)
    arrow((33.5, 10.3), (33.5, 18.7), col=OURS)
    arrow((64.3, 27.5), (67.7, 29.5)); ax.text(65.8, 31.2, "data", fontsize=6.5, color=MUTED, ha="center")
    arrow((64.3, 23.5), (67.7, 16.5)); ax.text(64.2, 18.0, "pilots", fontsize=6.5, color=MUTED, ha="right")
    arrow((76.5, 20.3), (76.5, 24.7)); ax.text(77.5, 22.2, r"$\hat{\mathbf{U}},\hat{\mathbf{V}}_s$", fontsize=7, color=INK)
    arrow((85.3, 29.5), (88.7, 29.5))
    ax.plot([42.3, 76.5], [4.5, 4.5], color=MUTED, lw=0.9, ls=(0, (3, 2)))
    arrow((76.5, 4.5), (76.5, 10.7), col=MUTED, ls=(0, (3, 2)))
    ax.text(61, 5.6, "timing $\\mathbf{T}_s$ known to receiver", fontsize=7, color=MUTED, ha="center")
    ax.plot([14, 52], [16.5, 16.5], color=MUTED, lw=0.8, ls=(0, (4, 3)))
    ax.text(18, 14.6, "direct path blocked", fontsize=6.8, color=MUTED, ha="left")
    ax.text(50, 37.5, r"$\mathbf{y} = \mathbf{H}\mathbf{x} + \mathbf{w}$,     "
            r"$\mathbf{H} = \sum_s \mathbf{T}_s \circ (\mathbf{F}\mathbf{U})(\mathbf{F}\mathbf{V}_s)^{\top}$",
            fontsize=8, color=INK, ha="center")
    save(fig, "fig0_system_model")


# ---------------------------------------------------------------------------
# Fig. 6 - NMSE vs SNR, and Fig. 7 - BER vs SNR (week4_curves.m)
# ---------------------------------------------------------------------------
def fig_curves():
    path = RES / "week4_curves.json"
    if not path.exists():
        print("skip curves: run matlab/week4_curves.m first"); return
    d = load("week4_curves")
    snr = np.array(d["snr_db"])
    nm = np.array(d["nmse_db_mean"]); ci = np.array(d["nmse_db_ci95"])
    sty = {"ours": ("Ours (timing only)", OURS, "o"), "bem": ("CE-BEM, 25 harmonics", BEM, "s"),
           "toeplitz": ("Toeplitz mixing model", TOEP, "^")}
    fig, ax = plt.subplots(figsize=(3.6, 2.9))
    for key in ["ours", "bem", "toeplitz"]:
        i = d["est_names"].index(key); lab, col, mk = sty[key]
        ax.errorbar(snr, nm[i], yerr=ci[i], color=col, marker=mk, markevery=2, capsize=2,
                    label=lab, mec="white", mew=0.8)
    ax.set_xlabel("SNR (dB)"); ax.set_ylabel("NMSE of channel estimate (dB)")
    ax.set_title(f"Channel estimation, $N_p$ = {d['Np']} pilots", loc="left", color=INK)
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout(); save(fig, "fig6_nmse_vs_snr")

    ber = np.array(d["ber"]); floor = 0.5 / d["bits_per_point"]
    rsty = {"onetap": ("One-tap (conventional)", ONETAP, "D", "-"),
            "toeplitz": ("MMSE, Toeplitz estimate", TOEP, "^", "-"),
            "bem": ("MMSE, CE-BEM estimate", BEM, "s", "-"),
            "ours": ("MMSE, ours (2 pilots)", OURS, "o", "-"),
            "ours_dd": ("MMSE, ours + 1 DD pass", "#1a4f94", "v", "-"),
            "perfect": ("MMSE, perfect channel", INK, None, "--")}
    fig, ax = plt.subplots(figsize=(3.6, 3.1))
    for key in ["onetap", "toeplitz", "bem", "ours", "ours_dd", "perfect"]:
        i = d["rx_names"].index(key); lab, col, mk, ls = rsty[key]
        y = np.where(ber[i] > 0, ber[i], np.nan)
        ax.semilogy(snr, y, color=col, marker=mk, markevery=2, ls=ls, label=lab, mec="white", mew=0.8,
                    lw=1.6 if key == "perfect" else 2.0)
    ax.set_ylim(max(floor, 5e-5), 1)
    ax.set_xlabel("SNR (dB)"); ax.set_ylabel("Bit error rate (16-QAM, uncoded)")
    ax.set_title("Detection", loc="left", color=INK)
    ax.legend(frameon=False, loc="lower left", fontsize=7)
    fig.tight_layout(); save(fig, "fig7_ber_vs_snr")

    c = d["constellation_20dB"]
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.4), sharex=True, sharey=True)
    for ax, key, title, col in [(axes[0], "onetap", "One-tap", ONETAP), (axes[1], "ours", "Ours (2 pilots)", OURS),
                                (axes[2], "perfect", "Perfect channel", INK)]:
        ax.scatter(c[key]["re"], c[key]["im"], s=3, color=col, alpha=0.55, lw=0)
        lv = np.array([-3, -1, 1, 3]) / np.sqrt(10)
        gx, gy = np.meshgrid(lv, lv)
        ax.scatter(gx, gy, s=14, facecolors="none", edgecolors=INK, lw=0.6)
        ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.6, 1.6); ax.set_aspect("equal")
        ax.set_title(title, loc="left", color=INK); ax.set_xlabel("In-phase")
    axes[0].set_ylabel("Quadrature")
    fig.suptitle("Equalised 16-QAM symbols at 20 dB (circles: ideal points)", x=0.01, ha="left",
                 fontsize=10, color=INK)
    fig.tight_layout(); save(fig, "fig8_constellation")


# ---------------------------------------------------------------------------
# Fig. 9 - equaliser complexity and memory vs bandwidth (operation counts)
# ---------------------------------------------------------------------------
def fig_complexity():
    Ms = np.array([256, 512, 1024, 2048, 3276, 4096, 8192])
    Q, iters = 16, 25
    dense = (4 / 3) * Ms.astype(float) ** 3                        # form H^H H + Cholesky, complex MACs
    nf = 2 ** np.ceil(np.log2(4 * Ms))
    cg = iters * 2 * (Q * 2) * (nf * np.log2(nf) + 2 * nf)       # 2 operator applies per CG iteration
    mem_dense = 16 * Ms.astype(float) ** 2 / 1e6                    # MB, complex double
    mem_free = 16 * (Q * 2 * nf + 4 * Ms) / 1e6
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))
    ax = axes[0]
    ax.loglog(Ms, dense, color=ONETAP, marker="D", label="Dense MMSE solve", mec="white")
    ax.loglog(Ms, cg, color=OURS, marker="o", label="Matrix-free CG (25 it.)", mec="white")
    ax.axvline(3276, color=MUTED, lw=0.8, ls=":"); ax.text(3276, 3e7, " n78\n 100 MHz", fontsize=7, color=MUTED)
    ax.set_xlabel("Active subcarriers $M$"); ax.set_ylabel("Complex multiply-adds per symbol")
    ax.set_title("Computation", loc="left", color=INK); ax.legend(frameon=False, loc="upper left")
    ax = axes[1]
    ax.loglog(Ms, mem_dense, color=ONETAP, marker="D", label="Dense $\\mathbf{H}$", mec="white")
    ax.loglog(Ms, mem_free, color=OURS, marker="o", label="Matrix-free operator", mec="white")
    ax.axvline(3276, color=MUTED, lw=0.8, ls=":")
    ax.set_xlabel("Active subcarriers $M$"); ax.set_ylabel("Memory (MB)")
    ax.set_title("Memory", loc="left", color=INK); ax.legend(frameon=False, loc="upper left")
    for a in axes:
        a.set_xticks([256, 512, 1024, 2048, 4096, 8192], ["256", "512", "1024", "2048", "4096", "8192"])
        a.minorticks_off()
    fig.suptitle("Operation counts (not MATLAB wall-clock time)", x=0.01, ha="left", fontsize=10, color=INK)
    fig.tight_layout(); save(fig, "fig9_complexity")


# ---------------------------------------------------------------------------
# Background figures (redrawn from the published model and gate1/gate2 data)
# ---------------------------------------------------------------------------
def gamma(f, f0, xr, xi):
    """Single-pole element response of Iudice et al. (Table I rates in Hz)."""
    xr, xi = 2 * np.pi * xr, 2 * np.pi * xi
    return -1 + 2 * xr / (1j * 2 * np.pi * (f - f0) + xr + xi)


def fig_bg_dispersion():
    fc = 3.594e9
    f = fc + np.linspace(-49.14e6, 49.14e6, 801)
    g0 = gamma(f, 3.471e9, 97.5e6, 8.9e6); g1 = gamma(f, 3.710e9, 128.4e6, 9.0e6)
    x = (f - fc) / 1e6
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.4))
    ax = axes[0]
    ax.plot(x, 20 * np.log10(abs(g0)), color=OURS, label="State 1 ($c_0$, 4 V)")
    ax.plot(x, 20 * np.log10(abs(g1)), color=BEM, label="State 2 ($c_1$, 19.25 V)")
    ax.set_title("Reflection strength", loc="left", color=INK); ax.set_ylabel("$|\\Gamma|$ (dB)")
    ax.legend(frameon=False, fontsize=7, loc="lower center")
    ax = axes[1]
    ax.plot(x, np.degrees(np.unwrap(np.angle(g0))), color=OURS)
    ax.plot(x, np.degrees(np.unwrap(np.angle(g1))), color=BEM)
    ax.set_title("Reflection phase", loc="left", color=INK); ax.set_ylabel("Phase (deg)")
    ax = axes[2]
    ax.plot(x, np.mod(np.degrees(np.angle(g0 / g1)), 360), color=INK)
    ax.axhline(180, color=MUTED, lw=0.8, ls=":")
    ax.set_title("Phase difference of the states", loc="left", color=INK); ax.set_ylabel("Degrees")
    for a in axes: a.set_xlabel("Offset from 3.594 GHz (MHz)")
    fig.tight_layout(); save(fig, "figA_element_dispersion")


def fig_bg_harmonics():
    d = load("gate1")["coupling"]
    h = np.array(d["h"]); dis = np.array(d["disp_db"]); idl = np.array(d["ideal_db"])
    fig, ax = plt.subplots(figsize=(4.6, 2.7))
    w = 0.38
    ax.bar(h - w/2, np.maximum(idl, -60) + 60, w, bottom=-60, color="#b9b8b2", label="Ideal ±1 switch")
    ax.bar(h + w/2, np.maximum(dis, -60) + 60, w, bottom=-60, color=OURS, label="Real element (dispersive)")
    ax.set_ylim(-60, 0); ax.set_xticks(h[::2])
    ax.set_xlabel("Harmonic index $h$ (shift in subcarriers)")
    ax.set_ylabel("Harmonic strength (dB)")
    ax.text(0.19, dis[h == 0][0] + 1.5, f"{dis[h == 0][0]:.0f}", ha="center", fontsize=6.5, color=INK)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=7)
    ax.set_title("Where one subcarrier's energy lands (2-slot pattern)", loc="left", color=INK)
    fig.tight_layout(); save(fig, "figB_harmonic_leakage")


def fig_bg_energy():
    g = load("gate2"); snr = g["snr_db"]
    cols = {"500": TOEP, "625": BEM, "750": OURS, "875": "#7a7975"}
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8), gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    for dkey, col in cols.items():
        r = g[f"S1_D{dkey}"]; lab = f"duty {int(dkey)/1000:g}"
        ax.plot(snr, r["rate_mmse"], color=col, marker="o", ms=4, mec="white", label=f"Joint MMSE, {lab}")
        ax.plot(snr, r["rate_tin"], color=col, ls=(0, (4, 2)), lw=1.5)
    ax.plot(snr, g["S1_D750"]["rate_mf_bound"], color=INK, ls=":", lw=1.3, label="Upper bound (all energy)")
    ax.set_xlabel("SNR (dB)"); ax.set_ylabel("Rate (bit/s/Hz)")
    ax.set_title("Solid: joint equaliser.  Dashed: one-tap", loc="left", color=INK)
    ax.legend(frameon=False, fontsize=6.8, loc="upper left")
    ax = axes[1]
    duties = [0.625, 0.75, 0.875]
    meas = [g[f"S1_D{int(d*1000)}"]["energy_gain_db"] for d in duties]
    pred = [-10 * np.log10((2 * d - 1) ** 2) for d in duties]
    x = np.arange(len(duties))
    ax.bar(x - 0.18, pred, 0.34, color="#b9b8b2", label="Predicted $1/(2D-1)^2$")
    ax.bar(x + 0.18, meas, 0.34, color=OURS, label="Measured")
    for xi, m in zip(x, meas): ax.text(xi + 0.18, m + 0.3, f"{m:.1f}", ha="center", fontsize=7, color=INK)
    ax.set_xticks(x, [f"{d:g}" for d in duties]); ax.set_xlabel("Duty cycle $D$")
    ax.set_ylabel("Energy-collection gain (dB)")
    ax.set_title("Theory vs simulation", loc="left", color=INK)
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout(); save(fig, "figC_energy_collection")


if __name__ == "__main__":
    fig_estimators()
    fig_model_knowledge()
    fig_fullband()
    fig_fair()
    fig_delay()
    fig_system()
    fig_curves()
    fig_complexity()
    fig_bg_dispersion()
    fig_bg_harmonics()
    fig_bg_energy()
