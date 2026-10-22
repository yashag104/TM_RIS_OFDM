# Summary, review step 1 (workstreams A2 + B1-B4), 2026-10-22

Tested: physical clock/timing models for claim C-C, with per-hop bulk delays (tau_b 50-300 ns, tau_u 30-300 ns),
receiver synced to the first cascaded arrival; 30 shared-path channels per point, paired t-intervals.
Scripts: week9_timing_physics (T1, T2), week9_timing_landscape (T3), week9_prior_sensitivity (T4).

| test | criterion (pre-registered) | outcome |
|---|---|---|
| T1 identities | rel. error < 1e-10 | **PASS**: max 6.9e-15. RIS locked to the incident signal: bulk delays cancel exactly after sync. RIS locked to the BS clock: the BS-RIS delay is exactly a clock offset d = -tau_b/T_slot. |
| T2 mode (i) | nominal and evidence within 0.3 of matched | **PASS**: evidence - matched +0.001 +-0.007 (20 dB), +0.040 +-0.082 (30 dB); with 1 ns per-edge jitter +0.015 / +0.086 +-0.06 |
| T2 mode (ii) | evidence within 0.3 of known offset | **PASS**: +0.011 +-0.019 (20 dB), +0.033 +-0.066 (30 dB). Without a search the loss is 1.62 / 4.23 bit/s/Hz, so the search is needed when the RIS follows the BS clock. |
| T3 (B4 overfitting) | held-out minimiser within 0.005 slot of evidence maximiser in >= 90% | **PASS**: 93% at 20 and 30 dB |
| T4 (B3 prior reference) | loss vs known within 0.5 dB for all 4 variants | **FAIL**: t0 = +100 ns loses 6.68 bit/s/Hz (20 dB). t0 = -100 ns +0.31 dB, tau_max x0.5 -0.27 dB, x1.5 +0.07 dB pass. |

Findings that change the paper:
1. "Residual search cannot recover the offset" is FALSE as stated. With the full estimator on a fine grid (0.01 slot),
   in-sample residual search is not significantly worse than the evidence: gap difference +0.021 +-0.032 (20 dB),
   +0.092 +-0.140 (30 dB); 1/30 channels with error > 0.01 slot vs 0/30 for the evidence. The week-7 failure
   came from the coarse procedure (8-iteration estimator on a 0.025 grid, refinement only within +-0.01), not
   from the criterion. The evidence's real advantages: smoother landscape (works on the coarse grid), uses only
   the linear first stage (no ALS per candidate, much cheaper), and no failures in 30 channels.
2. The overfitting explanation (B4) is supported: the held-out residual and the evidence pick the same offset in
   93% of channels and give the same rate (difference -0.001 +-0.029 at 30 dB).
3. "Shifted responses are representable, projection error ~1e-4" holds only for small shifts: the advanced input
   hop has projection error 1.1e-5 at 50 ns, 2.7e-4 at 100 ns, 1.4e-2 at 200 ns, 0.82 at 500 ns (L = 16,
   tau_max = 0.3 us). The near-flatness of the residual is a small-offset (< ~100 ns) effect.
4. The prior window must contain the sync reference (start at or before 0). A window starting after 0 cannot
   hold both hops after sync (their first arrivals sum to 0) and fails completely (+20 dB). This is a
   configuration rule, but it is a real fragility and must be stated.
5. Physically, the offset is either ~0 (RIS locked to the incident signal: no search needed; the search then
   costs +0.04 +-0.08 bit/s/Hz at 30 dB, not significant) or equal to the BS-RIS delay (search needed and works).

Claims that change: C-C is revised (see CLAIMS_LEDGER.md). Per the stopping rule, work stops here for review.

## Step 1b (fresh seeds), 2026-10-22
F1 (fix for the T4 failure): receiver sync error +-150 ns. Loss vs reference (perfect sync, known offset):
default window +0.22/+1.16, fixed margin window +0.21/+0.79, window start chosen by the evidence jointly with the
offset +0.050 +-0.046 / +0.046 +-0.078 bit/s/Hz (20/30 dB). PASS for the evidence-chosen window (criterion:
mean <= 0.165 and upper CI <= 0.33 at 30 dB); fixed margin FAILS. The paper's Algorithm 2 now chooses the window
start jointly with the offset. F2 (residual vs evidence regimes) pending.
