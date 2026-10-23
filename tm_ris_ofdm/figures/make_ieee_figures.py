"""Figures for the short IEEE paper (paper/ieee/), drawn from saved results/*.json and closed-form expressions.

No simulation runs here. Sized for IEEEtran: single column 3.5 in, double column 7.16 in.
Run:  .venv/bin/python -I tm_ris_ofdm/figures/make_ieee_figures.py
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
OUT = HERE.parent / "paper" / "ieee"

# Reference categorical slots 1-4 (validated: CVD and normal-vision pass; contrast WARN -> markers + legend + tables)
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, MUTED, GRID, NEUTRAL = "#0b0b0b", "#52514e", "#e4e3df", "#8a8985"

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.titlesize": 8,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.4,
    "lines.markersize": 4.5, "savefig.dpi": 300, "figure.facecolor": "white",
})


def load(name):
    return json.loads((RES / f"{name}.json").read_text())


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------------------
# Fig. 2: why the problem is hard (closed form): element dispersion and harmonic tail
# ---------------------------------------------------------------------------
def fig_mechanism():
    f = np.linspace(3.544e9, 3.644e9, 401)
    def gamma(f0, xr, xi):                   # Iudice et al., OpenRIS cell, xi given as xi/2pi in Hz
        return -1 + 2 * xr / (1j * (f - f0) + xr + xi)
    g1, g2 = gamma(3.471e9, 97.5e6, 8.9e6), gamma(3.710e9, 128.4e6, 9.0e6)
    K, s = 8, np.array([1] * 6 + [2] * 2)
    h = np.arange(-20000, 20001)
    def t(st):
        k = np.where(s == st)[0]
        return np.exp(-1j * np.pi * h / K) / K * np.sinc(h / K) * np.exp(-2j * np.pi * np.outer(h, k) / K).sum(1)
    e = np.abs(t(1) - t(2)) ** 2
    H = np.unique(np.round(np.logspace(0, np.log10(1000), 60)).astype(int))
    out = np.array([1 - e[np.abs(h) <= x].sum() for x in H])

    fig, ax = plt.subplots(1, 2, figsize=(3.5, 1.55))
    fg = (f - 3.594e9) / 1e6
    ax[0].plot(fg, np.degrees(np.angle(g1)), color=C1, label=r"$\Gamma_1$")
    ax[0].plot(fg, np.degrees(np.angle(g2)), color=C2, ls="--", label=r"$\Gamma_2$")
    ax[0].set_ylim(-150, 150)
    ax[0].set_xlabel("$f-f_c$ (MHz)"); ax[0].set_ylabel(r"$\angle\Gamma_s(f)$ (deg)")
    ax[0].legend(frameon=False, loc="best", handlelength=1.6)
    ax[0].set_title("(a) element response", loc="left")
    ax[1].loglog(H, 100 * out, color=C1)
    for x in (12, 40, 128):
        y = 100 * (1 - e[np.abs(h) <= x].sum())
        ax[1].plot([x], [y], "o", color=C1, mec="white", mew=0.6)
        ax[1].annotate(f"{y:.1f}%", (x, y), textcoords="offset points", xytext=(3, 3), fontsize=6.5, color=INK)
    ax[1].set_xlabel("harmonics kept, $|h|\\leq H$"); ax[1].set_ylabel("energy outside (%)")
    ax[1].set_title("(b) harmonic tail", loc="left")
    fig.tight_layout(w_pad=0.6)
    save(fig, "fig_mechanism")


# ---------------------------------------------------------------------------
# Fig. 3: fair comparison (week10_fair_baselines.m): loss to perfect CSI vs pilot symbols
# ---------------------------------------------------------------------------
def fig_fair():
    d = load("week10_fair_baselines")
    nps, ests = d["Np"], d["estimators"]
    style = [("ours", "Proposed (timing only)", C1, "o", "-"),
             ("ideal", "Ideal-element LMMSE", C3, "^", "-"),
             ("bemall", "CE-BEM, all harmonics", C2, "s", "-"),
             ("bem128", r"CE-BEM, $|h|\leq128$", C2, "s", ":"),
             ("omp", "OMP", C4, "D", "-")]
    fig, axes = plt.subplots(1, 4, figsize=(7.16, 1.75), sharey=True)
    panels = [("phys30", 0, "(a) shared path, 20 dB"), ("phys30", 1, "(b) shared path, 30 dB"),
              ("TDLC100", 0, "(c) TDL-C 100 ns, 20 dB"), ("TDLC100", 1, "(d) TDL-C 100 ns, 30 dB")]
    for ax, (ch, isn, title) in zip(axes, panels):
        r = np.array(d[ch]["rate_trials"], float)[:, :, 1, :, isn, 0]          # robust-residual eq., known G
        rp = np.array(d[ch]["rate_perfect_trials"], float)[:, isn]
        for key, lab, col, mk, ls in style:
            g = rp[:, None] - r[:, ests.index(key), :]
            m = np.nanmean(g, 0); n = np.sum(np.isfinite(g), 0)
            ci = 1.96 * np.nanstd(g, 0, ddof=1) / np.sqrt(np.maximum(n, 1))
            ok = np.isfinite(m)
            ax.errorbar(np.array(nps)[ok], np.maximum(m[ok], 1e-3), yerr=ci[ok], color=col, marker=mk, ls=ls,
                        label=lab, capsize=1.5, mec="white", mew=0.5, elinewidth=0.8)
        ax.axhline(0.3, color=MUTED, lw=0.7, ls="--")
        ax.set_xscale("log", base=2); ax.set_yscale("log")
        ax.set_xticks(nps, [str(n) for n in nps]); ax.set_ylim(3e-3, 12)
        ax.set_xlabel("pilot symbols $N_p$"); ax.set_title(title, loc="left")
    axes[0].set_ylabel("loss to perfect CSI\n(bit/s/Hz)")
    axes[3].text(33, 0.3, "0.3", fontsize=6.5, color=MUTED, va="center")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.08, 1, 1), w_pad=0.4)
    save(fig, "fig_fair")


# ---------------------------------------------------------------------------
# Fig. 4: coded BLER (week11_coded_stats.m), 2000 blocks per point
# ---------------------------------------------------------------------------
def fig_bler():
    d = load("week11_coded_stats")
    snr, names, b = d["snr_db"], d["names"], np.array(d["bler"], float)
    style = [("ours", "Proposed", C1, "o", "-"), ("perfect", "Perfect CSI", INK, None, "--"),
             ("bem12", r"CE-BEM $|h|\leq12$", C2, "s", "-"), ("bem40_robust", r"CE-BEM $|h|\leq40$, robust", C2, "s", ":"),
             ("onetap", "One-tap", NEUTRAL, "v", "-")]
    fig, ax = plt.subplots(figsize=(3.5, 2.0))
    for key, lab, col, mk, ls in style:
        y = b[names.index(key)]
        ok = y > 0
        ax.semilogy(np.array(snr)[ok], y[ok], color=col, marker=mk, ls=ls, label=lab, mec="white", mew=0.5)
    ax.axhline(0.1, color=MUTED, lw=0.7, ls="--")
    ax.set_xlabel("SNR (dB)"); ax.set_ylabel("block error rate"); ax.set_ylim(1e-3, 1.5)
    fig.legend(*ax.get_legend_handles_labels(), frameon=False, loc="lower center", fontsize=6.5, ncol=3,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.16, 1, 1))
    save(fig, "fig_bler")


if __name__ == "__main__":
    fig_mechanism()
    fig_fair()
    fig_bler()
