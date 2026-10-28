# Summary STD: standard-format curves with intermediate references (week12_standard_curves)

Pre-registered on 2026-10-24 (PREREG_LOG Step 8). The run used 30 new channels per point.
- **Data:** `results/week12_standard_curves.json`.
- **Figures:** `paper/ieee/fig_std_snr`, `fig_std_alpha` and `fig_std_pilots`, made with `figures/make_ieee_figures.py std`.

## Why this was run
The paper's figures showed ours next to perfect CSI and every baseline at an error floor, with nothing in between.

This run adds the formats used in the literature (NMSE, rate and BER against SNR, and NMSE against pilots). It also adds references that sit between ours and the baselines:
- an exact-element LMMSE, which knows the true Γ_s(f);
- the same estimator with an element model 5, 10 or 20 MHz off;
- CE-BEM with 16 times the pilots;
- an element non-ideality continuum.

## S1: against SNR (shared path, Np = 2), loss to perfect CSI in bit/s/Hz

| SNR (dB) | 0 | 5 | 10 | 15 | 20 | 25 | 30 |
|---|---|---|---|---|---|---|---|
| Ours | 0.02 | 0.02 | 0.03 | 0.06 | 0.08 | 0.12 | 0.17 |
| Exact-element LMMSE | 0.00 | 0.01 | 0.02 | 0.04 | 0.10 | 0.21 | 0.36 |
| Model 20 MHz off | 0.00 | 0.01 | 0.02 | 0.05 | 0.13 | 0.28 | 0.56 |
| Ideal-element LMMSE | 0.01 | 0.02 | 0.07 | 0.19 | 0.50 | 1.19 | 2.32 |
| CE-BEM, all harmonics, Np = 32 | 0.21 | 0.21 | 0.30 | 0.37 | 0.44 | 0.50 | 0.57 |
| OMP | 0.14 | 0.26 | 0.69 | 2.10 | 4.10 | 5.82 | 7.50 |
| CE-BEM, abs(h) ≤ 128, Np = 2 | 0.91 | 1.86 | 3.20 | 4.77 | 6.43 | 8.12 | 9.81 |
| One-tap | 0.79 | 1.71 | 3.02 | 4.54 | 6.15 | 7.79 | 9.44 |

**NMSE (dB) at 20 dB:**

| Estimator | NMSE |
|---|---|
| Ours | −31.4 |
| Exact-element LMMSE | −26.7 |
| Model 20 MHz off | −26.3 |
| Ideal-element LMMSE | −22.6 |
| CE-BEM, Np = 32 | −22.5 |
| OMP | −13.0 |
| CE-BEM, Np = 2 | −0.3 |

**Uncoded 16-QAM BER at 25 dB** (shared path):
- perfect CSI 3.3e-5;
- ours 7.5e-5;
- exact-element LMMSE 8.8e-4;
- ideal-element LMMSE 6.0e-4;
- CE-BEM with 32 pilots 3.2e-4.

At 30 dB, neither perfect CSI nor ours made an error in 3.1e5 bits.

**TDL-C 100 ns, Np = 4**, loss at 20 / 30 dB:

| Estimator | 20 dB | 30 dB |
|---|---|---|
| Ours | 0.20 | 0.74 |
| Exact-element LMMSE | 0.26 | 1.02 |
| Ideal-element LMMSE | 0.63 | 2.55 |
| CE-BEM, Np = 32 | 1.21 | 1.64 |

**Honest finding: there is a crossover.**
- At low SNR, the estimators that know the element model are slightly better than ours (crossover 17–22 dB depending on the channel set; Step 12 adds the 22 dB shared-path case). On the shared path the loss is 0.00–0.04 for them against 0.02–0.06 for ours. On TDL-C at 0 dB the NMSE is −13.7 dB for them against −7.2 dB for ours.
- From 20 dB upward ours is better, and the margin grows with SNR.
- **Why:** ours fits a nonlinear (ALS) low-rank model. At very low SNR this costs a little against a linear LMMSE that has the right structure. At high SNR, the low rank (r = 8 against 256 unknowns) removes more noise than any element knowledge can.

## S2: element non-ideality continuum (shared path, Np = 2)
Here H_α = (1−α)·H_ideal + α·H_OpenRIS. The table gives the loss to perfect CSI in bit/s/Hz.

| α | 0 | 0.1 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|---|
| Ideal-element LMMSE, 20 dB | 0.053 | 0.058 | 0.098 | 0.303 | 0.474 | 0.471 |
| Ours, 20 dB | 0.078 | 0.078 | 0.078 | 0.076 | 0.072 | 0.073 |
| Exact-element LMMSE, 20 dB | 0.053 | 0.053 | 0.053 | 0.053 | 0.052 | 0.053 |
| Ideal-element LMMSE, 30 dB | 0.21 | 0.26 | 0.60 | 1.70 | 2.28 | 2.27 |
| Ours, 30 dB | 0.18 | 0.18 | 0.18 | 0.16 | 0.15 | 0.15 |

Together with the BER panel, this shows a continuum rather than two extremes:
- with ideal elements, the literature's estimator works as well as ours;
- as the elements become realistic it degrades steadily, while ours stays flat.

## S3: NMSE (dB) against pilots at 20 dB (shared path)

| Np | 1 | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|---|
| Ours | −29.3 | −32.1 | −36.0 | −39.3 | −42.6 | −45.4 |
| Exact-element LMMSE | −26.1 | −28.8 | −31.7 | −34.5 | −37.3 | −40.2 |
| Model 10 MHz off | −26.1 | −28.7 | −31.4 | −34.0 | −36.6 | −38.8 |
| Ideal-element LMMSE | −22.4 | −23.2 | −23.9 | −24.2 | −24.5 | −24.6 |
| CE-BEM, abs(h) ≤ 128 | 0.0 | −0.3 | −2.8 | −13.6 | −20.2 | −24.2 |
| CE-BEM, all harmonics | n/a | n/a | n/a | −9.3 | −17.9 | −23.0 |
| OMP | −11.3 | −12.9 | −15.0 | −15.9 | −16.3 | −16.4 |

## Pre-registered checks
- **(a) Baseline sanity: passes, but the check is weaker than intended.** At α = 0 the ideal-element LMMSE equals the exact-element LMMSE, because at α = 0 they are the same estimator by construction. So the check holds trivially (pre-registration oversight, noted here). The informative comparison is the ideal-element LMMSE against ours at α = 0: 0.053 vs 0.078 at 20 dB and 0.21 vs 0.18 at 30 dB. These are within the 95% intervals, so the baseline is as good as ours when its assumption holds.
- **(b) Ours against the exact-element LMMSE at α = 1, Np = 2, 20 dB: PASS.** The losses are 0.073 and 0.053, a difference of 0.02, which is under 0.1. At 30 dB ours is better (0.15 vs 0.20).
- **(c) Crossover α:** the ideal-element LMMSE loss first exceeds 0.3 bit/s/Hz at α = 0.5 at 20 dB (0.303) and at α = 0.25 at 30 dB (0.60).
- **(d) Mismatch ordering: holds.** At 20 dB the NMSE orders as exact ≤ 5 MHz ≤ 10 MHz ≤ 20 MHz off, on both channels. The effect is small (≤ 0.4 dB at 20 dB and ≤ 1.4 dB at 30 dB on the shared path). An approximately known element model is nearly as good as an exact one; assuming an ideal switch is what fails.
