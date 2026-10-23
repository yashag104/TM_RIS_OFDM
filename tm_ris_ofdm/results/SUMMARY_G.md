# Summary G: theory checks and ablation of the pipeline (2026-10-22)

## G1 effective degrees of freedom — `week11_theory.json` (30 channels, G = 1, M = 256)
dof = tr(hat matrix) of the linear first stage: 19 / 22 / 24 / 27 (20 dB) and 23 / 27 / 30 / 33 (30 dB) for
Np = 1 / 2 / 4 / 8, against 256 / 512 / 1024 / 2048 observations and 2L^2 = 512 nominal unknowns. Why two pilot
symbols suffice: the DPSS prior concentrates the problem on ~20-33 effective parameters.

## G2 Bayesian MSE bound vs achieved NMSE (same LMMSE, dB)
| Np | bound 20 | achieved 20 | bound 30 | achieved 30 |
|---|---|---|---|---|
| 1 | -29.0 | -23.5 | -33.3 | -30.8 |
| 2 | -32.4 | -26.4 | -39.4 | -34.3 |
| 4 | -35.6 | -29.3 | -43.9 | -37.8 |
| 8 | -38.6 | -32.1 | -47.5 | -41.0 |
Achieved error is 2.5-6.5 dB above the bound: the bound holds under the assumed DPSS prior, which is not the true
channel distribution. Reported as prior mismatch, not as estimator suboptimality.

## G3 ablation / G3b fallback removal — `week10_ablation.json`, `week11_nofallback.json`
Linear stage alone does most of the work; low rank adds 0.04-0.12; one DD pass more. The pilot-fit fallback was an
overfitting trap on TDL-C (chosen in 37-57% of channels at 30 dB, worse gap); removed after a fresh-seed test: never
worse, up to -0.25 +-0.17 bit/s/Hz better. Algorithm 1 in the paper has no fallback.

## G4: not done (stated as a limitation).
