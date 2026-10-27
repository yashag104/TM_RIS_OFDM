# Reproducing and checking the results yourself

This guide explains how to check every number and figure in the IEEE paper (`paper/ieee/main.tex`) on your own machine. There are three levels, from minutes to days:

1. **Level 1 (about 10 minutes):** rebuild the figures and print every number from the saved results. No simulation.
2. **Level 2 (minutes to hours per experiment):** rerun one experiment in MATLAB and compare it with the saved results.
3. **Level 3:** read the pre-registration log and the claims ledger, to check that each claim was tested as it was registered.

---

## 0. What you need

| Tool | Version used | Used for |
|---|---|---|
| MATLAB | R2026a, base MATLAB only (no toolboxes) | all simulations (`matlab/*.m`) |
| Python | 3.12, with numpy 2.5, scipy 1.18, matplotlib 3.11 and pymupdf | figures and analysis (`figures/*.py`) |
| LaTeX | any TeX Live, MiKTeX or Overleaf | the paper (`paper/ieee/main.tex`) |

The Python environment is the WSL virtual environment at `/mnt/c/Users/hp/RIS/.venv`. Run Python from a WSL (Ubuntu) terminal.

If you need to recreate the environment:

```bash
cd /mnt/c/Users/hp/RIS
python3 -m venv .venv
.venv/bin/pip install numpy scipy matplotlib pymupdf
```

---

## Level 1: rebuild figures and numbers from the saved results (no simulation)

Every experiment saved its raw per-channel results as JSON in `results/`, and these files are committed to git. The commands below read only those files. All commands run from `/mnt/c/Users/hp/RIS` in WSL.

```bash
# Figures 2 and 6 (mechanism, coded BLER) and the old comparison figure
.venv/bin/python -I tm_ris_ofdm/figures/make_ieee_figures.py

# Figures 3, 4 and 5 (SNR curves, non-ideality continuum, NMSE vs pilots)
.venv/bin/python -I tm_ris_ofdm/figures/make_ieee_figures.py std

# Table IV (held-out pilot overhead); also rewrites the table rows inside main.tex
.venv/bin/python -I tm_ris_ofdm/figures/analyse_heldout.py

# Every number quoted from the week-12 runs: Figs. 3-5, the learned-baseline table,
# the low-SNR test and the pre-registered checks
.venv/bin/python -I tm_ris_ofdm/figures/analyse_week12.py
```

Then compile the paper. Run it twice so the references resolve:

```bash
cd tm_ris_ofdm/paper/ieee
pdflatex main.tex
pdflatex main.tex
```

On Overleaf instead: upload the whole `paper/ieee/` folder (`main.tex`, `IEEEtran.cls`, the `algorithm*.sty` files and every `fig_*.pdf`) and compile.

**What to check.** The figures in `paper/ieee/` should match the PDF. Every number in the analysis printout should match the text and tables. Many older experiments store their summary numbers inside the JSON itself, for example `"m1_eps0p05_20dB": [mean, ci]`, so you can open the file and read them directly.

---

## Level 2: rerun an experiment in MATLAB

### How to run a script
1. Open MATLAB on Windows.
2. Set the Current Folder to `C:\Users\hp\RIS\tm_ris_ofdm\matlab`.
3. Type the script name in the Command Window, for example `week12_standard_curves`. The script prints progress per channel.

From a terminal instead (this is how the runs were made):

```bash
cd /mnt/c/Users/hp/RIS/tm_ris_ofdm/matlab
"/mnt/c/Program Files/MATLAB/R2026a/bin/matlab.exe" -batch "week12_standard_curves" > run.log 2>&1
```

### Protect the saved results first
Most scripts write their JSON to `results/` and **overwrite** the saved file. Since the saved results are committed to git, do one of the following:

- **Week-12 scripts:** give them a different output name. For example, to run with 1 channel as a quick smoke test:
  `ntrial = 1; outname = 'my_check'; week12_standard_curves`
  This writes `results/my_check.json` and leaves the original alone. The same works for `week12_c6_lowsnr` (also `nsel`, `ntrains`) and `week12_c6b` (also `ntrain`).
- **Other scripts:** rerun them, then compare with `git diff --stat tm_ris_ofdm/results/`. To restore the originals, run `git checkout tm_ris_ofdm/results/`.

### Will I get the same numbers?
- **Seeds:** every script fixes its random seeds with `rng(...)` per channel. The same MATLAB version therefore draws the same channels and gives the same results.
- **Small differences:** multithreaded BLAS can change the last digit or so of iterative results. These differences are far inside the reported 95% intervals.
- **Timing results** (equaliser wall-clock, H2) depend on the machine and on what else is running. Run them only on an idle machine.

### Which script produces which part of the paper

| Paper item | MATLAB script | Result file | Analysis or figure script | Approx. run time |
|---|---|---|---|---|
| Fig. 2 (element response, harmonic tail) | none (closed form) | n/a | `make_ieee_figures.py` | seconds |
| Figs. 4–5 (NMSE / loss / BER vs SNR, 10–30 dB in 2 dB steps, Tang format; one file per graph) | `week12_snr_tang.m` (Step 11: 10:5:40; Step 12: `snrs=10:2:30`, `icsel=1` / `2`) | `week12_snr_tang.json`, `week12_snr_tang_fine_{phys30,TDLC100}.json` | `make_ieee_figures.py tang`, `analyse_snr_tang.py` | about 4.5 h (4.7 min per channel) |
| Fig. 4 (continuum), Fig. 5 (NMSE vs pilots) | `week12_standard_curves.m` | `week12_standard_curves.json` | `make_ieee_figures.py std`, `analyse_week12.py` | about 3 h (6.4 min per channel) |
| Fig. 6 (coded BLER) | `week11_coded_stats.m` | `week11_coded_stats.json` | `make_ieee_figures.py` | hours |
| Comparison text: pilot counts, paired gains over the best model-based baseline | `week10_fair_baselines.m` | `week10_fair_baselines.json` | summary fields inside the JSON | hours |
| Learned-baseline table (C6, shift, self-trained) | `week12_c6_lowsnr.m` (part 2), `week12_c6b.m` | `week12_c6_lowsnr.json`, `week12_c6b.json` | `analyse_week12.py` | about 3–4 h and 1.5–2 h |
| Low-SNR rank test (FAIL) | `week12_c6_lowsnr.m` (part 1) | same | `analyse_week12.py` | included above |
| Robustness table (model violations) | `week9_mismatch.m`, `week10_parametric.m`, `week11_cdl_rank.m` | same names `.json` | fields inside the JSON | hours |
| Timing table | `week9_timing_physics.m`, `week9_resid_vs_evidence.m`, `week9_sync_window.m`, `week9_timing_landscape.m` | same names `.json` | fields inside the JSON | hours |
| Full-band SER | `week11_fullband30.m` | `week11_fullband30.json` | fields inside the JSON | hours |
| Overhead table IV (held-out) | `week11_overhead_heldout.m`; `week11_overhead_heldout_part.m` with `icsel = 2, 3, 4` for TDL-A/C/D | `week11_overhead_heldout*.json` | `analyse_heldout.py` | about 8 h with the three parts in parallel |
| Equaliser timing (break-even about 12 symbols) | `week11_eq_amortised.m` (idle machine) | `week11_eq_amortised.json` | printed by the script | 25 min |
| Hardware (CFO, phase noise, rise time) | `week6_hardware_fix.m`, `week6_phasenoise.m`, `week11_rise_small.m` | same names `.json` | fields inside the JSON | hours |
| Emission (ACLR 36.8 / 52.7 dB) | none (analytic) | `k3_oob.json` | `figures/k3_emission.py` (recomputes and compares; ACLR within 0.01 dB) | seconds |
| Cost of time modulation | `week11_motivation.m` | `week11_motivation.json` | fields inside the JSON | about 1 h |

For any claim not listed here, `results/CLAIMS_LEDGER.md` names the script and result file behind it.

### Core functions, if you want to read the method
- `coupling_matrix.m`: the true channel through the time-modulated RIS (from Iudice et al.).
- `element_response.m`, `params_openris.m`: the OpenRIS element model.
- `est_general.m`: **our estimator** (Algorithm 1, timing only, low rank).
- `est_bem_cg.m`: the CE-BEM baseline.
- `dpss_basis.m`: the delay prior.
- `phys_channel.m`, `tdl38901.m`: the channel generators.

The ideal-element, exact-element, OMP and learned-prior estimators are local functions at the end of the scripts that use them.

### Suggested first exercise (about 15 minutes)
1. In MATLAB, run `ntrial = 1; outname = 'my_check'; week12_standard_curves`. This is one channel per part.
2. In the log, find the line `S2 ch 1 | 20 dB gap (ideal, genie, ours) at a = 0: ...`.
3. Compare it with the first channel in `results/week12_standard_curves.json`. The numbers should be identical, because the seeds are the same.

---

## Level 3: check that the claims were tested honestly

- `results/PREREG_LOG.md`: every test with its pass/fail rule, the date and the git commit, written **before** the run. Use `git log -p tm_ris_ofdm/results/PREREG_LOG.md` to confirm that each entry was committed before its result file. Compare with `git log -- tm_ris_ofdm/results/<file>.json`.
- `results/CLAIMS_LEDGER.md`: every claim with its status (OK, FIXED, NARROWED, WITHDRAWN, FAIL) and its evidence.
- `results/SUMMARY_*.md`: one summary per workstream, including the negative results: the learned prior beating ours, the low-SNR rank FAIL, the 120 km/h limit and the CG speed correction.
- `results/SOFTENED_CLAIMS.md`: claims that were weakened or withdrawn after review.
