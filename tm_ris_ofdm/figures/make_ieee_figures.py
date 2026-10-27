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


# ---------------------------------------------------------------------------
# Standard-format results (week12_standard_curves.m): families share a hue, variants differ by dash/marker
#   element-model LMMSE family (green): ideal +-1 (solid), exact element (dashed), 20 MHz off (dotted)
#   CE-BEM family (orange): all harmonics with Np = 32 (solid), |h|<=128 with ours' Np (dotted)
# ---------------------------------------------------------------------------
STD = [("ours", "Proposed (timing only)", C1, "o", "-", None),
       ("genie", "Exact-element LMMSE", C3, "^", "--", "white"),
       ("mm20", "Element model 20 MHz off", C3, "v", ":", None),
       ("ideal", "Ideal-element LMMSE", C3, "^", "-", None),
       ("bemall_Np32", r"CE-BEM, all $h$, $N_p=32$", C2, "s", "-", None),
       ("bem128", r"CE-BEM, $|h|\leq128$", C2, "s", ":", "white"),
       ("omp", "OMP", C4, "D", "-", None),
       ("onetap", "One-tap", NEUTRAL, "x", "-", None)]


def _std_plot(ax, x, ys, key, lab, col, mk, ls, mfc, logy=True):
    kw = dict(color=col, marker=mk, ls=ls, label=lab, mec=col if mfc == "white" else "white", mew=0.7)
    if mfc == "white":
        kw["mfc"] = "white"
    ax.plot(x, ys, **kw)


def fig_std_snr():
    s = load("week12_standard_curves")["S1"]["phys30"]; e = s["estimators"]; snr = np.array(s["snr_db"])
    nm = np.array(s["nmse"], float); rt = np.array(s["rate"], float); rp = np.array(s["rate_perfect"], float)
    be = np.array(s["ber"], float); bp = np.array(s["ber_perfect"], float)
    fig, axes = plt.subplots(1, 3, figsize=(7.16, 1.7))
    for key, lab, col, mk, ls, mfc in STD:
        i = e.index(key)
        if key != "onetap":
            _std_plot(axes[0], snr, 10 * np.log10(np.nanmean(nm[:, i, :], 0)), key, lab, col, mk, ls, mfc)
        _std_plot(axes[1], snr, np.maximum(np.nanmean(rp - rt[:, i, :], 0), 3e-3), key, lab, col, mk, ls, mfc)
        y = np.nanmean(be[:, i, :], 0); ok = y > 0
        if key not in ("bem128", "onetap"):
            _std_plot(axes[2], snr[ok], y[ok], key, lab, col, mk, ls, mfc)
    y = np.nanmean(bp, 0); ok = y > 0
    axes[2].plot(snr[ok], y[ok], color=INK, ls="--", lw=1.0, label="Perfect CSI")
    axes[0].set_ylabel("NMSE (dB)"); axes[0].set_title("(a) channel estimate", loc="left")
    axes[1].set_yscale("log"); axes[1].set_ylim(5e-3, 15); axes[1].axhline(0.3, color=MUTED, lw=0.7, ls="--")
    axes[1].set_ylabel("loss to perfect CSI\n(bit/s/Hz)"); axes[1].set_title("(b) robust-MMSE rate loss", loc="left")
    axes[2].set_yscale("log"); axes[2].set_ylim(1e-5, 0.6)
    axes[2].set_ylabel("uncoded 16-QAM BER"); axes[2].set_title("(c) detection", loc="left")
    for ax in axes:
        ax.set_xlabel("SNR (dB)"); ax.set_xticks(snr)
    h, l = axes[1].get_legend_handles_labels(); h2, l2 = axes[2].get_legend_handles_labels()
    fig.legend(h + h2[-1:], l + l2[-1:], loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.04),
               fontsize=6.5)
    fig.tight_layout(rect=(0, 0.15, 1, 1), w_pad=0.5)
    save(fig, "fig_std_snr")


W2, H2 = 3.45, 2.35            # one graph at full column width (two per row in a figure*)
W3, H3 = 2.33, 1.95            # one graph at a third of the text width


def save_legend(handles, labels, name, ncol, width=7.0):
    """Legend strip as its own file, shared by the separately saved graphs of one figure."""
    fig = plt.figure(figsize=(width, 0.1 + 0.17 * np.ceil(len(labels) / ncol)))
    fig.legend(handles, labels, loc="center", ncol=ncol, frameon=False, fontsize=7.5)
    save(fig, name)


def fig_snr_tang(smax=30):
    """Estimators vs SNR in the format of Tang et al. Figs. 10-11 (week12_snr_tang): one file per channel and metric."""
    R = load("week12_snr_tang")
    metrics = [("nmse", "channel NMSE", None), ("loss", "loss to perfect CSI (bit/s/Hz)", None),
               ("ber4", "uncoded BER, QPSK", None), ("ber16", "uncoded BER, 16-QAM", None)]
    for c in ("phys30", "TDLC100"):
        fine = RES / f"week12_snr_tang_fine_{c}.json"          # Step 12: 2 dB grid over 10-30 dB, when complete
        s = json.loads(fine.read_text())[c] if fine.exists() else None
        s = s if s is not None and len(s["nmse"]) == 30 else R[c]           # partial file while running: use Step 11
        e = s["estimators"]; snr = np.array(s["snr_db"]); keep = snr <= smax; snr = snr[keep]
        nm = np.array(s["nmse"], float)[:, :, keep]; rt = np.array(s["rate"], float)[:, :, keep]
        rp = np.array(s["rate_perfect"], float)[:, keep]
        be = np.array(s["ber"], float)[:, :, keep]; bp = np.array(s["ber_perfect"], float)[:, keep]
        for mname, ylab, ylim in metrics:
            fig, ax = plt.subplots(figsize=(W2, H2))
            for key, lab, col, mk, ls, mfc in STD:
                if key not in e:
                    continue
                i = e.index(key)
                if mname == "nmse":
                    y = np.nanmean(nm[:, i], 0)
                elif mname == "loss":
                    y = np.maximum(np.nanmean(rp - rt[:, i], 0), 3e-3)
                else:
                    y = np.nanmean(be[:, i, :, 0 if mname == "ber4" else 1], 0)
                ok = y > 0
                _std_plot(ax, snr[ok], y[ok], key, lab, col, mk, ls, mfc)
            if mname.startswith("ber"):
                y = np.nanmean(bp[:, :, 0 if mname == "ber4" else 1], 0); ok = y > 0
                ax.plot(snr[ok], y[ok], color=INK, ls="--", lw=1.1, label="Perfect CSI")
            if mname == "loss":
                ax.axhline(0.3, color=MUTED, lw=0.8, ls="--")
            ax.set_yscale("log"); ax.margins(x=0.02, y=0.06); ax.grid(True, which="major", lw=0.4)
            ax.set_xlabel("SNR (dB)"); ax.set_ylabel(ylab); ax.set_xlim(snr[0] - 0.5, snr[-1] + 0.5)
            ax.set_xticks(np.arange(10, smax + 1, 5))
            fig.tight_layout()
            if mname == "ber4" and c == "phys30":
                save_legend(*ax.get_legend_handles_labels(), "fig_snr_legend", ncol=4)
            save(fig, f"fig_snr_{c}_{mname}")


def fig_std_alpha():
    s = load("week12_standard_curves")["S2"]; e = s["estimators"]; a = np.array(s["alpha"])
    rt = np.array(s["rate"], float); rp = np.array(s["rate_perfect"], float)
    for k in range(2):
        fig, ax = plt.subplots(figsize=(W2, H2))
        for key, lab, col, mk, ls, mfc in STD:
            if key not in e or key in ("bem128", "onetap"):
                continue
            g = rp[:, :, k] - rt[:, e.index(key), :, k]
            m = g.mean(0); h = 1.96 * g.std(0, ddof=1) / np.sqrt(g.shape[0])
            kw = dict(color=col, marker=mk, ls=ls, label=lab, mec="white", mew=0.5, capsize=2, elinewidth=0.8)
            if mfc == "white":
                kw.update(mfc="white", mec=col)
            ax.errorbar(a, np.maximum(m, 1e-3), yerr=h, **kw)
        ax.axhline(0.3, color=MUTED, lw=0.8, ls="--")
        ax.set_yscale("log"); ax.set_ylim(0.02, 4); ax.set_xticks([0, 0.25, 0.5, 0.75, 1], ["0", "0.25", "0.5", "0.75", "1"])
        ax.set_xlabel(r"element non-ideality $\alpha$"); ax.set_ylabel("loss to perfect CSI (bit/s/Hz)")
        fig.tight_layout()
        if k == 0:
            save_legend(*ax.get_legend_handles_labels(), "fig_alpha_legend", ncol=4)
        save(fig, f"fig_alpha_{s['snr_db'][k]}dB")


def fig_std_pilots():
    """NMSE vs pilot symbols in the style of Tang et al. Figs. 6-8: linear MSE on a logarithmic axis."""
    s = load("week12_standard_curves")["S3"]; e = s["estimators"]; nps = np.array(s["Np"])
    nm = np.array(s["nmse"], float)
    fig, ax = plt.subplots(figsize=(3.5, 2.2))
    rows = [r for r in STD if r[0] in e] + [("bemall", r"CE-BEM, all $h$", C2, "s", "-", None)]
    rows = [r for r in rows if r[0] != "bemall_Np32"]
    for key, lab, col, mk, ls, mfc in rows:
        if key not in e:
            continue
        y = np.nanmean(nm[:, e.index(key), :], 0); ok = np.isfinite(y)
        _std_plot(ax, nps[ok], y[ok], key, lab, col, mk, ls, mfc)
    ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.set_xticks(nps, [str(n) for n in nps])
    ax.set_ylim(1e-5, 2); ax.grid(True, which="both", lw=0.4)
    ax.set_xlabel("pilot symbols $N_p$"); ax.set_ylabel("channel NMSE")
    fig.legend(*ax.get_legend_handles_labels(), loc="lower center", ncol=3, frameon=False, fontsize=6.3,
               bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.17, 1, 1))
    save(fig, "fig_std_pilots")


# ---------------------------------------------------------------------------
# Element response in the format of Iudice et al. Figs. 2-4 (single-pole model, Table I parameters)
# ---------------------------------------------------------------------------
F0 = np.array([3.471e9, 3.710e9]); XR = 2 * np.pi * np.array([97.5e6, 128.4e6]); XI = 2 * np.pi * np.array([8.9e6, 9.0e6])
PHI0 = np.radians(-10.7)                      # common renormalization phase of Iudice Table I (Fig. 3)


def gam(f, k):
    return -1 + 2 * XR[k] / (1j * 2 * np.pi * (f - F0[k]) + XR[k] + XI[k])


def fig_element():
    """Three separate graphs as in Iudice et al. Figs. 2-4: magnitude, phase, phase difference."""
    f = np.linspace(3.3e9, 3.8e9, 1001); fg = f / 1e9
    g = [np.exp(1j * PHI0) * gam(f, k) for k in (0, 1)]
    lab = [r"$c_0$ (4 V), state 1", r"$c_1$ (19.25 V), state 2"]

    def axis():
        fig, ax = plt.subplots(figsize=(W3, H3))
        ax.axvspan(3.544, 3.644, color=GRID, alpha=0.9, lw=0, zorder=0)
        ax.axvspan(3.594 - 0.00384, 3.594 + 0.00384, color=NEUTRAL, alpha=0.45, lw=0, zorder=0)
        ax.set_xlabel("frequency (GHz)"); ax.set_xlim(3.3, 3.8); ax.set_xticks([3.3, 3.4, 3.5, 3.6, 3.7, 3.8])
        return fig, ax
    fig, ax = axis()
    for k, (col, ls) in enumerate(((C1, "-"), (C2, "--"))):
        ax.plot(fg, 20 * np.log10(np.abs(g[k])), color=col, ls=ls, label=lab[k])
    ax.set_ylabel(r"$|\Gamma|$ (dB)"); ax.set_ylim(-2.5, 0.05); ax.legend(frameon=False, loc="lower left", fontsize=6.5)
    fig.tight_layout(); save(fig, "fig_element_mag")
    fig, ax = axis()
    for k, (col, ls) in enumerate(((C1, "-"), (C2, "--"))):
        ax.plot(fg, np.degrees(np.angle(g[k])), color=col, ls=ls, label=lab[k])
    ax.set_ylabel(r"$\angle\Gamma$ (deg)"); ax.set_yticks([-150, -100, -50, 0, 50, 100, 150])
    ax.legend(frameon=False, loc="lower left", fontsize=6.5)
    fig.tight_layout(); save(fig, "fig_element_phase")
    fig, ax = axis()
    ax.plot(fg, np.degrees(np.angle(g[0] / g[1])) % 360, color=INK, lw=1.2, label="OpenRIS cell")
    ax.axhline(180, color=C3, ls=":", lw=1.2, label=r"ideal $180^\circ$")
    ax.set_ylabel("phase difference (deg)"); ax.legend(frameon=False, loc="lower left", fontsize=6.5)
    fig.tight_layout(); save(fig, "fig_element_dphase")


# ---------------------------------------------------------------------------
# Coupling-matrix structure: reproduction of Iudice et al. Fig. 5 next to this paper's configuration
# ---------------------------------------------------------------------------
def fig_coupling():
    d = load("fig_coupling_data")
    for k, name in (("iudice_K2", "fig_coupling_K2"), ("ours_K8", "fig_coupling_K8")):
        Hm = np.array(d[k]["re"]) + 1j * np.array(d[k]["im"])
        fig, ax = plt.subplots(figsize=(W3, H3))
        im = ax.imshow(20 * np.log10(np.abs(Hm) + 1e-12), cmap="viridis", vmin=-60, vmax=0, origin="upper",
                       interpolation="nearest")
        ax.set_xlabel("input subcarrier $m$"); ax.set_ylabel(r"output tone $\bar m$")
        ax.set_xticks([0, 20, 40, 60]); ax.set_yticks([0, 20, 40, 60]); ax.grid(False)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.set_label(r"$|H_{\bar m m}|$ (dB)", fontsize=7)
        fig.tight_layout(); save(fig, name)
    K, h = 8, np.arange(-20000, 20001)
    sq = np.array([1] * 6 + [2] * 2)

    def t(st):
        k = np.where(sq == st)[0]
        return np.exp(-1j * np.pi * h / K) / K * np.sinc(h / K) * np.exp(-2j * np.pi * np.outer(h, k) / K).sum(1)
    e = np.abs(t(1) - t(2)) ** 2
    H = np.unique(np.round(np.logspace(0, np.log10(1000), 60)).astype(int))
    out = np.array([1 - e[np.abs(h) <= x].sum() for x in H])
    fig, ax = plt.subplots(figsize=(W3, H3))
    ax.loglog(H, 100 * out, color=C1)
    for x in (12, 40, 128):
        y = 100 * (1 - e[np.abs(h) <= x].sum())
        ax.plot([x], [y], "o", color=C1, mec="white", mew=0.6)
        ax.annotate(f"{y:.1f}%", (x, y), textcoords="offset points", xytext=(3, 3), fontsize=6.5, color=INK)
    ax.set_xlabel(r"harmonics kept, $|h|\leq H$"); ax.set_ylabel("energy outside (%)")
    fig.tight_layout(); save(fig, "fig_coupling_tail")


if __name__ == "__main__":
    import sys
    if "struct" in sys.argv:
        fig_element(); fig_coupling(); sys.exit()
    if "tang" in sys.argv:
        fig_snr_tang(); sys.exit()
    if "std" in sys.argv:
        fig_snr_tang(); fig_std_alpha(); fig_std_pilots(); sys.exit()
    fig_mechanism()
    fig_fair()
    fig_bler()
