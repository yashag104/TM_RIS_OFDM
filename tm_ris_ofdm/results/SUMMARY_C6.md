# Summary C6: learning-based baseline, distribution shift, self-trained prior, low-SNR rank

The work was pre-registered in two steps:
- **Step 9 (2026-10-25):** `week12_c6_lowsnr.json`, the learned prior and the low-SNR rank.
- **Step 10 (2026-10-25):** `week12_c6b.json`, distribution shift and the self-trained prior.

The Step 9 outcome was logged before Step 10 ran.

## Baseline
**Learned-prior LMMSE.** Build a PCA basis of vec(H) from N_train channels. Keep the components with eigenvalue at least 1e-6 of the largest. Use a zero-mean Gaussian prior with the sample eigenvalues and estimate with the LMMSE from the pilots.

This is the Bayes-optimal linear estimator for the learned statistics, and an MSE-trained network approaches it for near-Gaussian channels. It learns the element non-ideality from data.

## Step 9: matched training (N_train = 100, 300, 1000), 30 test channels
The learned prior beats ours at every SNR, every Np, on both channels, already from N_train = 100.

| Channel | Learned (1000) loss | Ours loss | Paired difference (ours worse by) |
|---|---|---|---|
| Shared path, 20 dB | 0.006 | 0.070 | +0.064 ± 0.009 |
| TDL-C, 20 dB | 0.074 | 0.285 | +0.211 ± 0.056 |
| TDL-C, 30 dB | 0.293 | 0.997 | +0.704 ± 0.258 |

**Why the learned prior does so well.** The simulated families are low-dimensional: 6 PCA components on the shared path and 50–52 on TDL-C.
- At M = 256 the band is 7.7 MHz. With a 30 ns delay spread and K = 10 dB, the cascade is nearly frequency-flat, and channels correlate at 0.92 on average.
- TDL uses fixed normalised tap delays.

## Step 10: shift and self-training (N_train = 300), loss to perfect CSI in bit/s/Hz at 20 / 30 dB

| Test | Ours | Learned, matched | Learned, shifted (trained on) | Self-trained (ours, Np = 8) |
|---|---|---|---|---|
| Shared 30 ns, K 10 dB | 0.071 / 0.143 | 0.007 / 0.031 | n/a | 0.029 / 0.118 |
| Shared 100 ns, K 3 dB | 0.124 / 0.468 | 0.054 / 0.247 | 0.346 / 1.416 (shared 30 ns) | n/a |
| TDL-C 100 ns | 0.185 / 0.668 | 0.058 / 0.213 | n/a | 0.078 / 0.362 |
| TDL-A 100 ns | 0.178 / 0.726 | 0.047 / 0.184 | 0.121 / 0.611 (TDL-C) | n/a |
| TDL-D 100 ns | 0.483 / 1.676 | 0.043 / 0.204 | 0.445 / 2.047 (TDL-C) | n/a |

**Readout (a), shift: mixed.** These are paired differences of the shifted prior against ours at 20 dB; positive means ours is better.

| Shift | 20 dB | 30 dB |
|---|---|---|
| Shared 30 → 100 ns | +0.22 ± 0.17 (ours better) | +0.95 ± 0.46 |
| TDL-C → TDL-D | +0.04 ± 0.07 (tie) | +0.37 ± 0.20 (ours better) |
| TDL-C → TDL-A | −0.06 ± 0.03 (learned better) | tie |

An offline prior breaks when the delay spread or the line-of-sight strength changes, and it transfers between similar non-line-of-sight profiles.

**Readout (b), self-trained prior: PASS, adopted as an optional mode.**
- At 20 dB it beats ours by 0.042 ± 0.012 on the shared path and 0.107 ± 0.027 on TDL-C.
- At 30 dB it gains 0.025 ± 0.034 (not significant) on the shared path and 0.305 ± 0.163 on TDL-C.
- It needs no true channels and no element model: 300 channels measured with 8 pilot symbols by our own estimator.
- Because it is collected in place, it matches the deployment distribution by construction.

## Step 9: low-SNR rank, FAIL
- **Rule A (rank per SNR from 20 selection channels, 30 held-out):** it improves 0–15 dB by +0.008–0.015 on the shared path and +0.034–0.050 on TDL-C, with CIs above zero. It is worse by 0.034 ± 0.035 at 20 dB on TDL-C, which breaks the "never worse by > 0.01" condition.
- **Rule B (BIC):** it always picks rank 1 and is far worse at high SNR.
- **Outcome:** the frozen rank stays, and the low-SNR crossover remains a stated limit.

## Consequences for the claims
- **Narrowed claim:** "beats every baseline we could build" is withdrawn. The new claim is that ours is the best receiver needing neither an element model nor training data.
- **Learned prior:** an LMMSE with a prior learned from about 100 or more channels of the deployment distribution is better. Trained on a different distribution, it can be much worse.
- **Self-trained prior:** this two-stage receiver (ours collects estimates, then a learned prior) keeps the no-model and no-ground-truth property and closes most of the gap.
- **Caveat for the paper:** our simulated channel families are low-dimensional. A real deployment is richer, which favours ours over offline-learned priors, but this has not been measured.
