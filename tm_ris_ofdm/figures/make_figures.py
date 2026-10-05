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


if __name__ == "__main__":
    fig_estimators()
    fig_model_knowledge()
    fig_fullband()
    fig_fair()
    fig_delay()
