# Summary C: baselines, ablation, parametric fit, rise time (2026-10-22)

All runs pre-registered in PREREG_LOG.md; 30 channels per point; paired t-intervals.

## C1/C4 (+E2, E3, E6): fair baselines — `week10_fair_baselines.json`
Rules: every estimator gets naive MMSE, robust MMSE (error power from the pilot residual) and robust MMSE (oracle
error power); CE-BEM with no harmonic cap (matrix-free LMMSE `est_bem_cg`, 7.3e-9 from the dense solver);
the Toeplitz straw man is replaced by the ideal-element model-aware LMMSE, H = (T1 - T2) o (F C F^T), same prior.
Channels: shared path 30 ns and TDL-C 100 ns; Np = 2, 4, 8, 16, 32; 20/30 dB.

Loss to perfect knowledge (robust-residual equaliser, known G), Np = 2 / 8 / 32:

| | shared 20 dB | shared 30 dB | TDL-C 20 dB | TDL-C 30 dB |
|---|---|---|---|---|
| ideal-element LMMSE | 0.48 / 0.44 / 0.43 | 2.27 / 2.16 / 2.14 | 0.85 / 0.50 / 0.41 | 3.05 / 2.12 / 1.95 |
| CE-BEM all harmonics | - / 3.42 / 0.43 | - / 5.10 / 0.55 | - / 5.58 / 1.17 | - / 9.04 / 1.55 |
| CE-BEM |h|<=128 | 6.55 / 2.34 / 0.46 | 9.89 / 3.60 / 1.51 | 6.05 / 5.31 / 1.02 | 9.30 / 8.83 / 2.14 |
| OMP | 4.16 / 3.13 / 2.92 | 7.49 / 6.67 / 6.44 | 4.13 / 2.90 / 2.59 | 7.61 / 6.12 / 5.80 |
| ours | 0.07 / 0.01 / 0.00 | 0.11 / 0.02 / 0.01 | 0.52 / 0.08 / 0.02 | 1.80 / 0.18 / 0.03 |

- E3: minimal Np for loss <= 0.3 at 20 dB: ours 2 (shared) / 4 (TDL-C); every baseline > 32.
- C4: ours minus best baseline per channel (min over 6 baselines x 3 equalisers): shared -0.41 +-0.02 / -2.16 +-0.04,
  TDL-C -0.30 +-0.06 / -1.14 +-0.25 at Np = 2 (20/30 dB); every Np excludes zero. PASS.
- E6: estimated-gain detector adds <= 0.15 to every loss; ranking unchanged. PASS.
- C1 sanity (oracle-robust CE-BEM not worse at 30 dB): FAIL on the shared-path channel for the truncated sets
  (|h| <= 12, 40; drop up to 0.28). Interpretation: a scalar error power cannot represent truncation error
  concentrated on the missing harmonics. PASS on TDL-C and for |h| <= 128 / all harmonics. Reported in the paper.
- New finding: the ideal-element approximation is fair at 20 dB (0.4-0.5 loss) and fails at 30 dB (2.0-2.3,
  flat in Np: model error). This is the quantitative case for not needing the element model.
- Runs launched before the fallback removal: conservative for ours.

## C2 ablation — `week10_ablation.json`
Linear first stage alone gives most of the gain (0.104 / 0.229 at Np 2); low rank r = 8: 0.065 / 0.108; +1 DD 0.019 / 0.029.

## C3 parametric element fit — `week10_parametric.json`
Six-parameter element fit not better than timing-only (+0.08 +-0.05 worse at 30 dB) and fails with element spread
(+0.73 +-0.43). Not adopted.

## C5 rise time — `week10_risetime.json`
CE-BEM's failure depends on sharp switching edges: with a 0.3-slot ramp CE-BEM |h|<=40, Np = 8 reaches 8.2 vs ours
8.56 (Np 2); at ramps <= 0.05 slot our advantage exceeds 4 bit/s/Hz. Banded equalisers do not recover even with smooth
edges (explanation in the paper corrected).

## Not done
C6 (learning-based estimator baseline): stated as a limitation.
