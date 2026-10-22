# Pilot-Based Estimation and Equalisation of Harmonic-Coupled OFDM Channels through Dispersive Time-Modulated RIS

MATLAB R2026a, base MATLAB only (no Communications Toolbox). Independent of the
Python FL/RIS codebase in the parent folder.

Literature map and gap: `RESEARCH_NOTES.md`.

## Status (2026-10-01)

| Gate | What | Result |
|---|---|---|
| 0 | Novelty: nobody equalises or estimates TM-RIS harmonic coupling at a legitimate receiver | Passed (searched; see notes) |
| 1 | Reproduce Iudice et al. (arXiv 2609.18360) and validate the closed form against an independent time-domain simulation | **Passed** |
| 2 | Joint MMSE vs treat-interference-as-noise (TIN), perfect H; pre-registered energy-collection prediction | **Passed** |
| 3 | Structured pilot estimator within ~1 dB of perfect H with Np <= 2 | **Passed** (physics-aware LMMSE + decision-directed); one marginal miss in the non-physical channel (see "Gate 3, closed") |

### Gate 1 (`matlab/gate1_reproduce.m`, `results/gate1.json`)
* Dispersion of the OpenRIS cell over 98 MHz. Differential phase 172–184 deg (paper:
  180 +/- 6.5), imbalance +/-0.57 dB (< 0.75), per-state phase variation 46.5/49.2 deg
  (~47/50), ripple 0.62/0.50 dB (~0.6). Matches.
* K = 2 coupling: odd harmonics on the sinc(h/2) envelope; even harmonics exactly zero;
  the DC harmonic appears only through dispersion, at -19 to -27 dB relative to h = 1.
* Closed-form coupling matrix vs direct time-domain simulation (filter -> gate -> FFT,
  Q = 4 elements, K = 4 random sequences, multipath on both hops): relative error
  0.27% median, 0.49% max (residual of sampled vs continuous gates).

### Gate 2 (`matlab/gate2_equalization.m`, `results/gate2.json`)
Q = 16 elements, K = 8 slots, 256 subcarriers (central 128 scored), 30 trials,
dispersive cell, 30 ns RMS multipath per hop.

| Scenario | duty D | energy-collection gain, measured | predicted (ideal) | rate @ 20 dB: TIN / MMSE / MF bound |
|---|---|---|---|---|
| S1 same sequence | 0.5   | 23.4 dB | inf (b0 = 0) | 0.01 / 6.22 / 6.54 |
| S1 | 0.625 | 11.6 | 12.0 | 0.10 / 6.14 / 6.45 |
| S1 | 0.75  | 5.9  | 6.0  | 0.42 / 6.19 / 6.52 |
| S1 | 0.875 | 2.4  | 2.5  | 1.19 / 6.11 / 6.43 |
| S2 random shifts | 0.5 | 26.4 | – | 0.01 / 5.59 / 6.61 |
| S2 | 0.75 | 6.5 | – | 0.62 / 5.13 / 6.62 |

* The measured gain matches the closed-form prediction 1/(2D-1)^2 to within 0.1–0.4 dB.
* MMSE tracks the MF bound within ~0.3 bit/s/Hz for S1. For random space-time shifts
  (S2) the gap grows to ~1.5 bit/s/Hz, because H is worse conditioned. That motivates
  sequence design for conditioning (week 3).
* With binary 50% duty (the paper's own example), TIN gets essentially nothing.
  Joint equalisation recovers about 95% of the MF-bound rate.

Caveat for the paper: the TIN gap depends on how much ICI reaches the user. A user in
a TM-IRS beam *designed* to keep h = 0 dominant will see a smaller gain than these
generic-position numbers. This must be shown explicitly, not hidden.

### Gate 3, first attempt: FAILED (`matlab/gate3_estimation.m`, `results/gate3.json`)
Model-agnostic estimator H = sum_h Shift_h diag(F c_h) with +/-8 or +/-12 harmonics and
delay-limited gains. NMSE floors at -14 to -16 dB whatever the pilot count, and the
rate after MMSE is only 2–3 bit/s/Hz against 5–6 with perfect H. Cause: rectangular
switching leaves a sinc^2 ~ 1/h^2 harmonic tail, so ~3–5% of the energy lies outside
any small harmonic set. Baselines at 20 dB: conventional one-tap 0.4, Toeplitz ~1.4–2.6,
banded unstructured LS ~0–0.6 bit/s/Hz.

### Gate 3, physics-aware estimator: PASSED for SNR <= 20 dB (`gate3b_physics.m`, `gate3c_tune.m`)
The legitimate receiver knows the control sequence and the element-response family,
so all harmonics are tied by a known matrix B:
S1: H = B o (Fd C Fs^T), with C the Ld x Ls two-hop delay matrix (64 unknowns).
Truth keeps all in-band harmonics.

| SNR | Np = 1 | Np = 2 | Np = 4 | perfect H |
|---|---|---|---|---|
| 10 dB | 3.21 | 3.21 | 3.26 | 3.31 |
| 20 dB | 6.21 | 6.30 | 6.23* | 6.45 |
| 30 dB | 7.9 | 8.3–8.6 | 9.56* | 9.6–9.8 |

(rate after MMSE, bit/s/Hz. Rows 20/30 dB at Np = 1, 2 use the tuned ridge (osr 2,
relative ridge 1e-3..1e-2); * = gate3b, untuned ridge.)

* Pre-registered criterion (Np <= 2, within ~1 dB of perfect) is met at 10 and 20 dB,
  but NOT at 30 dB, where 4 pilot symbols are needed. Residual NMSE floor ~ -30 dB from
  off-grid delays in the delay basis (DPSS basis is the planned fix).
* **Dispersion-aware B matters.** B from the ideal +/-1 model floors at -24.5 dB NMSE.
  The exact model reaches -52 dB. At 30 dB / Np = 4: 9.56 vs 7.34 bit/s/Hz.
* S2 (per-element random shifts, independent per-element channels): 1024 unknowns, so it
  needs Np >= 8 (5.0 vs 5.1 bit/s/Hz at 20 dB). The independent-per-element channel is
  Iudice's general form but unphysical for a compact surface. Co-located elements share
  path delays, which should cut the unknowns to ~64 + Q phases (week 3).

### Gate 3, closed (`gate3d_fix.m`, `gate3d_lmmse.m`, `gate3e_dd.m`)
One calibrated setting (DPSS basis, tau_max = 0.3 us, L = 16, LMMSE with DPSS-eigenvalue
prior), chosen on a calibration seed and evaluated on fresh seeds. Per-element sequences
use a rank-1 (a_q s_q^T) ALS estimator. Channels: Iudice's general independent
per-element TDL, and a physically consistent model (shared paths, geometric phases,
K = 10 dB).

Loss vs perfect-H MMSE, dB, Np = 2 pilot symbols (rate-equivalent):

| case | 10 dB | 20 dB | 30 dB pilots only | 30 dB + decision-directed (2 iters, 6 data symbols) |
|---|---|---|---|---|
| iid S1 | 0.06 | 0.45 | 1.06 | 0.3 |
| iid G4 | 0.29 | 0.58 | 1.87 | **1.2** |
| iid S2 | 0.30 | 0.56 | 1.31 | 0.7 |
| phys S1 | 0.04 | 0.16 | 0.47 | – |
| phys G4 | 0.29 | 0.57 | 1.17 | 0.5 |
| phys S2 | 0.38 | 0.56 | 1.19 | 0.5 |

Verdict: criterion met everywhere for the physical channel. The single miss is iid G4 at 30 dB
(1.2 dB, limited by 64-QAM decision errors at SER 2–4%) under the non-physical
independent-per-element channel. Reported as is.

## Week 3 results

### (a) Fair case: TM designed for the legitimate user (`week3_designed_user.m`)
Per-element polarity (sequence or its complement) co-phases h = 0 at the user;
random per-element shifts keep the signal scrambled elsewhere. Physical channel, Np = 2 estimate.

| duty | user | one-tap SIR | TIN rate 10/20/30 dB | joint MMSE (estimated) 10/20/30 dB |
|---|---|---|---|---|
| 0.875 | designed | +10.0 dB | 2.51 / 3.34 / 3.47 | 3.13 / 6.18 / 9.37 |
| 0.75  | designed | +3.5 dB | 1.44 / 1.70 / 1.74 | 2.77 / 5.63 / 8.80 |
| 0.625 | designed | -2.3 dB | 0.62 / 0.70 / 0.71 | 2.24 / 4.53 / 7.38 |
| 0.75  | generic  | -12.3 dB | 0.17 / 0.19 / 0.19 | 2.24 / 4.77 / 7.66 |

Even in the most TIN-friendly configuration, TIN saturates at its SIR (~3.5 bit/s/Hz).
The gain is modest at 10 dB (+0.6 bit/s/Hz, ~2 dB) and large at 20–30 dB.
Report it that way.

### (b) Uncoded BER, 16-QAM (`week3_banded_ber.m`)
* One-tap: BER 0.39–0.45, unusable without joint equalisation.
* Full MMSE with the estimated channel (Np = 2 + 2 DD iterations) matches perfect-H BER
  (S1: 2.9e-3 at 16 dB, 2e-5 at 20 dB, 0 at 24 dB).
* **Banded sliding-window MMSE fails.** With W = 16 it floors at 2.5e-3 (S1) / 3e-2 (S2): the
  1/h^2 harmonic tail reaches far outside any short window. Negative result, worth stating.
* S2 (random shifts) has a high BER even with perfect H (7.7e-3 at 28 dB), because of conditioning.

### (c) Exact fast operator + CG equaliser (`fast_coupling.m`, `week3_fast_eq.m`)
* Key identity: B_q(mbar, m) = sum_s T_{q,s}(mbar - m) Gamma_s(f_m). Each element's coupling
  is an exact sum of 2 Toeplitz matrices times the dispersive state responses, so H x and
  H' y cost O(Q M log M) by FFT, with no M x M matrix needed.
* Verified against the explicit matrix: relative error 7e-16.
* 20 CG iterations give the same SER as the direct full MMSE (perfect and estimated H).
* Honest timing: in MATLAB the dense direct solve is still faster in wall-clock up to
  M = 1024 (0.15 s vs 0.54 s). The advantage is asymptotic in operations and memory
  (M = 3276 dense H is 170 MB). The paper must state flops/memory, not "faster".

### (d) Channel-aware shift design (`week3_seq_design.m`): DROPPED from the paper (2026-10-06)
Reason: 8 trials only, and the gain is mostly the known space-time-coding harmonic steering. Kept here for the record.

Greedy per-element shift choice that maximises the legit user's MMSE rate:
5.58 -> 8.85 bit/s/Hz (identical shifts: 5.52). The eavesdropper's one-tap SIR stays scrambled
(-6.5 -> -8.2 dB). Most of the gain is harmonic *beamforming* (shifts set harmonic phases,
so power adds coherently at the user, at fixed Tx power and noise). That is the known
space-time-coding harmonic-steering idea (Zhang/Cui 2018). Only the joint MMSE-rate /
dispersive design is new. Supporting result, not a headline.

### (e) Element-model mismatch, and a model-free estimator (`week3_mismatch.m`, `week3_modelfree.m`)
Truth uses Table I. The receiver's B uses shifted resonances (df0) and +/-10% decay rates.
Np = 2 + 2 DD iterations. Loss vs perfect H, dB at 20 / 30 dB:

| receiver knowledge | S1 | S2 (rank-1) |
|---|---|---|
| exact model | 0.05 / 0.12 | 0.24 / 0.51 |
| df0 = 5 MHz | 0.53 / 2.59 | 0.24 / 0.49 |
| df0 = 10 MHz | 0.66 / 3.36 | 0.24 / 0.50 |
| df0 = 20 MHz | 0.55 / 2.91 | 0.24 / 0.51 |
| ideal +/-1 (frequency-flat literature model) | 1.41 / 6.57 | 0.24 / 0.53 |
| **model-free per-state (timing only)** | **0.08 / 0.11** | not run yet |

* The physics-aware S1 estimator passes the stated criterion (< 1 dB at 20 dB for df0 <= 10 MHz)
  but is fragile at 30 dB.
* **Model-free estimator.** B_q(mbar, m) = sum_s T_{q,s}(mbar - m) Gamma_s(f_m), and T_{q,s}
  depends only on the control timing. Gamma_s(f_m) multiplies the input side, so it is absorbed
  into a per-state input response: H = sum_s T_s o (F C_s F^T), with 2 L^2 unknowns.
  It needs NO element model, matches the exact-model estimator (0.08 / 0.11 dB), and
  removes the 3–6.5 dB loss of mismatched or frequency-flat models. This replaces
  "dispersion-aware B" as the headline estimator.
* S2 under mismatch is insensitive (losses dominated by estimation noise). A per-(element, state)
  model-free variant for S2 still needs a bilinear form (shared a_q, per-state s_{q,s}); to do.

### (f) Full n78 band, M = 3276, matrix-free (`week3_fullband*.m`)
Estimation regressors and equalisation use only the exact FFT operator (no M x M matrix).
* First run (full delay matrix C, L = 36). SER at 20 dB / 28 dB (S1): estimated 0.045 / 0.057
  vs true-H 0.033 / 0.018. Gap at high SNR. More DD passes made it WORSE (0.24):
  64-QAM decision errors propagate.
* Diagnosis (`week3_fullband_err.m`). CG converges by 25 iterations, so it is not the cause. The model is correct:
  operator NMSE reaches -45..-49 dB at 60 dB SNR. At 28 dB the full-C estimator is
  noise-limited at -23..-29 dB, because the two-hop delay matrix has ~(2NW)^2 effective unknowns,
  which grows quadratically with bandwidth.
* Fix: C is low rank (shared propagation paths: rank <= #paths). A rank-r ALS estimator
  C = U V^T has 2 L r unknowns: -39.0 dB at Np = 2 (rank 6 or 16, same), vs -26.3 dB for full C.
* SER with low rank (S1, L = 30):

| SNR / QAM | one-tap | CG true H | estimate, pilots only | estimate + 1 DD |
|---|---|---|---|---|
| 20 dB / 16-QAM | 0.84 | 0.0327 | 0.034 | 0.034 |
| 28 dB / 64-QAM | 0.96 | 0.0182 | 0.042 | 0.026 |

* S2 (per-element sequences) at full band, shared path-subspace estimator
  (`week3_fullband_s2.m`, a_q = U alpha_q, s_q = V beta_q, block ALS). Pre-registered bar:
  operator NMSE <= -30 dB at 28 dB AND SER(est) <= 1.5 x SER(true) at 20 dB.
  **FAILED.** SER at 20 dB is 0.059–0.062 vs 0.044 (1.3–1.4x, OK), but operator NMSE is -24 dB at
  both 20 and 28 dB (a floor, not noise), and SER at 28 dB is 0.17–0.18 vs 0.054.
  Diagnostic (`week3_fullband_s2_diag.m`, mean-B init + 3 random restarts, best by pilot
  residual): operator NMSE -26.9 dB at 28 dB and -26.5 dB at 60 dB. The range across restarts is
  -23 to -30 dB. SER at 28 dB: 0.044 vs 0.015 true; at 60 dB: 0.0125 vs 0.0028.
  Conclusion: not noise (no change 28 -> 60 dB), and not the model (some restarts reach -29/-30 dB).
  The cause is the non-convex 4-block ALS settling in local minima, and the pilot residual does not
  reliably pick the best restart. **Recorded as a limitation**: per-element sequences
  work at M = 256 but not yet at full band. Possible future fixes (not done): convex/lifted
  formulation, better-conditioned initialisation, more pilots.

### (g) Verification for publication (`week3_modelfree_confirm.m`, `week3_delayspread.m`)
* Headline re-run on a fresh seed (1801), 30 trials, 95% CIs. Loss vs perfect H, dB, at 10/20/30 dB:
  exact model 0.02/0.06/0.16; df0 = 10 MHz 0.06/0.47/2.44; ideal +/-1 0.18/1.42/6.60;
  **timing-only 0.04/0.10/0.21** (CI on rate loss <= +/-0.05 bit/s/Hz). Reproduces (e).
* Delay spread beyond the assumed support (estimator unchanged, tau_max = 0.3 us), timing-only,
  20 trials. Loss dB at 20/30 dB: 30 ns 0.10/0.16; 100 ns 0.12/0.31 (bar <= 1 dB at 20 dB: pass);
  200 ns 0.30/1.29 +/- 0.77.
* Code review: no estimator uses the true channel; calibration was on a separate seed;
  mismatched filters are scored on the true H. The receiver is assumed to know the noise
  variance and the control timing (state in the paper). Lemma 1 is a regrouping of
  Iudice et al. eq. (5) and is credited as such in `paper/letter.tex`.

## Week 4: strengthening for submission

### (a) Comparison with the literature's estimators (`week4_bem_compare.m`)
Same pilots and seeds, physical channel, S1, M = 256, 20 trials, 95% CIs. Rate (bit/s/Hz) of the
MMSE equaliser built from each estimate and scored on the true H:

| estimator | 20 dB, Np = 2 | 20 dB, Np = 8 | 30 dB, Np = 2 | 30 dB, Np = 8 |
|---|---|---|---|---|
| one-tap (conventional) | 0.42 | 0.42 | 0.43 | 0.43 |
| Toeplitz mixing (ideal model, as in Tao et al.) | 2.27 +/- 0.08 | 2.46 +/- 0.10 | 1.42 +/- 0.10 | 1.52 +/- 0.03 |
| CE-BEM, Q = 4 | 1.29 +/- 0.04 | 1.33 +/- 0.03 | 0.57 +/- 0.05 | 0.58 +/- 0.02 |
| CE-BEM, Q = 8 | 1.82 +/- 0.06 | 1.88 +/- 0.04 | 0.81 +/- 0.07 | 0.99 +/- 0.04 |
| CE-BEM, Q = 12 | 2.20 +/- 0.08 | 2.32 +/- 0.04 | 1.08 +/- 0.10 | 1.38 +/- 0.06 |
| **ours (timing-only)** | **6.52 +/- 0.02** | **6.58 +/- 0.02** | **9.68 +/- 0.08** | **9.87 +/- 0.02** |
| perfect H | 6.61 | 6.61 | 9.92 | 9.92 |

CE-BEM (the standard estimator for intra-symbol time variation) was given the same LMMSE/DPSS
prior as ours, so this is a strong version of it. It does not improve with pilots: the 1/h^2
harmonic tail is outside any truncated harmonic set. It also gets *worse* at 30 dB than at 20 dB,
because an MMSE filter that trusts a wrong model hurts more as noise falls.

### (b) Per-element sequences at full band: FIXED via the shift structure (`week4_fullband_s2_shift.m`)
Exact identity: for a cyclic shift k of the base sequence,
t_{q,s}[h] = t_{base,s}[h] e^{-j2 pi h k/K} = diag(D_k) T_s diag(D_k^*), with D_k = e^{-j2 pi m k/K}.
So H = sum_s sum_k (D_k F U) Lam_k (D_k^* ... F V_s) o T_s: U is shared (output side), V_1 and V_2
are per state (input side, absorbing dispersion), and Lam_k (r x r) is per shift value, so there
are at most K = 8 small matrices instead of 16 per-element pairs. Timing-only, no element model.
Given (U, V) the Lam step is a convex linear solve.
Same pre-registered bar as week3_fullband_s2. M = 3276, Np = 2, r = 6, 5 trials:

| | op NMSE | SER one-tap | SER true H | SER estimated | ratio |
|---|---|---|---|---|---|
| 20 dB, 16-QAM | -30.1 dB | 0.863 | 0.0543 | 0.0598 | **1.10** |
| 28 dB, 64-QAM | **-34.4 dB** (trials -32..-36) | 0.965 | 0.0397 | 0.0443 | 1.12 |

**PASSED** (bar: <= -30 dB at 28 dB, ratio <= 1.5 at 20 dB). The earlier limitation (4-block ALS
stalling at -24..-27 dB) is resolved. A bug was caught before the run (the output-side modulation
was applied twice) and fixed.

### (c) Full band, 10 trials with 95% CIs, timing-only vs CE-BEM (`week4_fullband_trials.m`)
M = 3276, S1, physical channel, Np = 2, Nd = 6. Timing-only low-rank estimator with a shared output
factor (H = sum_s T_s o ((F U)(F V_s)^T), r = 8). CE-BEM with Q = 8 and the same prior. SER:

| receiver | 20 dB / 16-QAM | 28 dB / 64-QAM |
|---|---|---|
| one-tap | 0.839 +/- 0.003 | 0.955 +/- 0.002 |
| true H (CG-MMSE) | 0.031 +/- 0.024 | 0.035 +/- 0.021 |
| timing-only, pilots | 0.034 +/- 0.025 (1.09x true) | 0.054 +/- 0.032 (1.54x) |
| timing-only, + 1 DD | 0.032 +/- 0.024 (1.05x) | 0.040 +/- 0.022 (1.14x) |

**Reporting rule (point 4):** always quote the pilots-only ratio next to the DD ratio. At 28 dB / 64-QAM
pilots-only is 1.54x here and was 2.3x in the earlier 3-trial run (`week3_fullband_lr`, DD 1.4x); the
headline must not quote only 1.05-1.14x. More trials are needed to settle it (point 5).
| CE-BEM, pilots / + 1 DD | 0.549 / 0.563 | 0.865 / 0.865 |

The CIs are dominated by channel-to-channel variation. Per trial, the estimate tracks the true-H SER
closely. The one pilots-only outlier (28 dB trial 9: 0.158 vs 0.039) drops to 0.073 after the DD
pass. This makes every full-band result timing-only (no element model).

### (d) Sensitivity to every design choice (`week5_param_sweep.m`, `est_general.m`, fig10)
Answers "why these settings?". One setting is varied at a time around the default (Q = 16, K = 8,
D = 0.75), M = 256, Np = 2, 12 channels per point. Same timing-only estimator everywhere
(`est_general`, which handles any set of G distinct switching patterns; identical and shifted
patterns are special cases). Gap of ours to perfect H, bit/s/Hz at 20 / 30 dB:
* Elements Q = 4, 16, 64, 144: 0.06–0.09 / 0.10–0.20.
* Slots K = 2, 4, 8, 16 (duty 0.5 and 0.75): 0.06–0.09 / 0.11–0.16.
* Duty D = 1/8 ... 7/8 (K = 8): 0.05–0.13 / 0.12–0.23.
* Random patterns, G distinct sequences over 16 elements:
  2 pilots 0.07, 0.12, 0.22, 0.34, 0.32 / 0.10, 0.19, 0.32, 0.78, 1.03 for G = 1, 2, 4, 8, 16;
  4 pilots 0.03, 0.05, 0.08, 0.17, 0.19 / 0.04, 0.08, 0.13, 0.32, 0.38.
CE-BEM stays at about 2.2 / 1.1 and one-tap at 0.01–1.2 bit/s/Hz throughout.
Conclusion: the results do not depend on Q, K or D. **Design rule:** with many distinct per-element
patterns (G >= 8), use 4 pilot symbols. (Rows for different G and Np use independent random
patterns and channels: compare gaps, not absolute rates.)

## Week 5-6: supervisor's improvement list

Backup before this work: `D:\RIS_backup_2026-10-06\` (folder copy + verified git bundle of all branches).

### (1) Coherence time and mobility (`week5_mobility.m`)
Previously assumed: block fading. Both hops are constant over pilots + data (as in Iudice et al.:
static geometry). Justification: switching period 33 us << coherence time; isolates TM leakage from Doppler.
Now tested: only the RIS -> user hop moves. Each path gets Doppler nu = fD cos(phi), both
symbol-to-symbol rotation and intra-symbol ICI. Estimate once (2 pilots), reuse for later symbols.
Gap to perfect H (bit/s/Hz) vs symbols since estimation n:

| speed | Tc | 20 dB, n = 2/10/40/150 | 30 dB, n = 2/10/40/150 | intra-symbol Doppler |
|---|---|---|---|---|
| 3 km/h | 42 ms (1188 sym) | 0.06/0.06/0.06/0.06 | 0.08/0.08/0.08/0.10 | -60 dB, 0 bit |
| 30 km/h | 4.2 ms (119 sym) | 0.07/0.09/0.26/0.54 | 0.23/0.36/0.82/1.27 | -47 dB, 0.008 bit |
| 120 km/h | 1.1 ms (30 sym) | 0.12/0.32/0.52/0.75 | 0.51/1.12/1.49/1.90 | -34 dB, 0.19 bit |

1-dB reuse length: >300 symbols (walking), ~40 (30 km/h, 20 dB) / ~8 (30 dB), ~8 (120 km/h, 20 dB),
under 2 at 120 km/h and 30 dB. Pilot overhead 2/N: <1% walking, ~5% urban at 20 dB (within 5G NR
DMRS density of 7-28% per slot). **Limitation:** vehicular speed at high SNR needs tracking. Natural
fix (not built): the input-side factor V (static BS -> RIS hop) can be reused; only U needs refreshing.

### (12) Direct link, two RIS panels, interference (`week6_directlink.m`)
* Unblocked direct path at -20 ... +10 dB relative to the RIS path. Ours, with NO direct-path term,
  vs perfect (bit/s/Hz, 20/30 dB): -10 dB 6.14/9.31 vs 6.25/9.53; 0 dB 5.28/8.13 vs 5.39/8.49;
  +10 dB 5.95/8.85 vs 6.01/9.01. The explicit direct-path estimator (`est_general_direct.m`) gives
  the same numbers. **Reason:** the state kernels sum to the identity (T_1 + T_2 = I, since the
  element is in one of the two states at all times), so an unswitched path is already inside the
  model (equal input response in both states). No change to the estimator is needed.
  One-tap rises to 2.9-3.1 with a strong direct path but stays far below.
* Two TM-RIS panels (independent channels, different patterns, G = 2): ours 5.54/8.67 vs perfect
  5.64/8.86 (rank 8; rank 16 the same). A receiver that models only one panel collapses to 1.15:
  **the receiver must know every panel's timing.**
* Co-channel OFDM interferer, unknown and treated as noise, SIR 20/10/0 dB: ours within
  0.03-0.17 bit/s/Hz of perfect at every SIR. Interference limits both receivers equally.

### (5) Hardware imperfections (`week6_hardware.m`), loss vs perfect, dB at 20/30 dB
* (a) element-to-element spread (resonance sigma 5-20 MHz, +/-5% losses): 0.20-0.28 / 0.28-0.46. Robust.
* (b) switching transients (linear ramp, rise time as a fraction of a slot) 0.05/0.1/0.2/0.3:
  1.06/2.78/3.73/3.88 at 20 dB, 4.76/10.19/11.37/11.73 at 30 dB. **Weakness:** the receiver assumes
  instant switching. Follow-up running (`week6_hardware_fix.m`): kernels from a calibrated or
  blindly estimated rise time.
* (c) receiver CFO 0.005/0.01/0.02/0.05 subcarrier: 0.26/0.32/0.52/1.70 at 20 dB, 0.72/1.00/2.45/7.50 at 30 dB.
  Requires residual CFO below ~1% of the spacing (standard OFDM synchronisation).
* (d) phase noise: the first model was WRONG (a fresh random phase per OFDM symbol) and its 8-25 dB
  losses are discarded. Redone with a continuous phase across symbols, with/without common-phase
  correction, in `week6_hardware_fix.m`.

### (6) Coded performance (`week6_coded.m`)
Rate-1/2 convolutional code (K = 7, 133/171 octal), soft Viterbi (checked: 0 errors noiseless;
0/506 at Eb/N0 1.9 dB vs 62 uncoded errors), random interleaver, 16-QAM Gray, max-log LLRs from each
receiver's OWN estimate. 506-bit packets, 30 channels x 4 packets per point. BLER vs SNR 4:2:24 dB:
* one-tap: 1.00 at every SNR (no packet ever correct)
* CE-BEM: best 0.31 (12-14 dB), then rises back to 0.97 at 24 dB (an overconfident wrong model)
* ours (2 pilots): 0.93 / 0.28 / 0.025 / 0 at 6 / 8 / 10 / 12 dB, 0 above
* perfect H: 0.93 / 0.21 / 0.008 / 0
Ours is within ~0.3 dB of perfect at 10% BLER. Throughput saturates at the code limit
(1.98 bit/s/Hz per subcarrier) from 12 dB; BEM peaks near 1.4; one-tap delivers 0.

### (10) Estimator complexity (`figures/complexity_estimator.py`, `results/week6_complexity.json`)
Operations per estimate (20 ALS iterations, r = 8, Np = 2), M = 256 / 1024 / 3276:
G = 1: 5.7e8 / 4.1e9 / 3.9e10; G = 16: 1.8e10 / 5.4e10 / 1.8e11. For comparison, one dense MMSE solve
at M = 3276 is 4.7e10 and the CG equaliser 4.2e8 per symbol. **The estimator dominates the cost.**
Re-estimating every ~40 symbols (30 km/h) at the full band is ~30 TOPS as written: not realistic.
Mitigations, now built and measured (`est_general_ws.m`, `week7_complexity_fix.m`,
`results/week7_complexity_fix.json`, Fig. 18). User at 30 km/h, refresh every 20 symbols, 8 blocks,
8 trials, S1, M = 256; gap to perfect (bit/s/Hz) at 20/30 dB, ops per update, MATLAB time:
| update | gap 20/30 dB | ops | vs full | time |
|---|---|---|---|---|
| full re-estimate (20 it, r = 8) | 0.07 / 0.17 | 1.2e9 | 1x | 307 ms |
| full, rank 6 | 0.07 / 0.16 | 7.2e8 | 1.7x | 215 ms |
| warm start, 2 iterations | 0.10 / 0.20 | 1.0e8 | 12x | 31 ms |
| V frozen, U + Lam only, 1 iteration | 0.10 / 0.40 | 1.1e7 | 105x | 8 ms |
| same, 1 pilot symbol | 0.13 / 0.49 | 6.0e6 | 197x | 5 ms |
**Warm start is the safe default: 12x fewer operations for +0.03 bit at both SNRs.** Freezing V is ~100x
cheaper and as good at 20 dB, but costs +0.23 bit at 30 dB: the factorisation does not split the
static and moving hops exactly, so V must be refreshed now and then (e.g. U-only updates between
occasional warm starts). Rank 6 is free. Applied to the full band, warm start brings the ~30 TOPS above to
roughly 2.5 TOPS; this is an extrapolation of the iteration count, not a full-band measurement.
The earlier "~50x" estimate was optimistic for the safe option (12x) and pessimistic for U-only (105x).

### (8) More baselines (`week6_baselines_pilots.m`), rate bit/s/Hz at 20/30 dB, 10 channels
ours 6.48/9.65 | CE-BEM 2.17/1.07 | OMP (harmonic x delay dictionary, stops at noise level) 2.32/1.19 |
unstructured banded LS with **32** orthogonal pilot symbols (16x our overhead) 0.85/0.54 | perfect 6.58/9.89.
Sparse recovery stalls like BEM (a truncated harmonic set). Brute-force LS is worst even with 16x
the pilots. Deep-learning estimators were not run: they need training data, and the comparison of
interest here is model structure. Listed as future work.

### (11) Pilot design, same number of pilot REs (512) and same pilot energy, rate at 20/30 dB
2 full pilot symbols 6.49/9.71 | comb-2 over 4 symbols, gaps EMPTY (NR DMRS-like) 6.52/9.79 |
comb-4 over 8 symbols, gaps empty 6.51/9.77 | comb-2 with DATA in the gaps 4.79/5.10 |
same + 1 decision-directed pass 6.36/8.82.
**Design rule:** standard comb pilots work as well as full pilot symbols if the neighbouring REs are
left empty (allowed in NR). Data next to pilots leaks into them through the switching; one DD pass
recovers most of it at 20 dB but not at 30 dB.

### (2) Timing knowledge (`week6_timing.m`)
Loss vs perfect (dB at 20/30 dB) with the receiver using the NOMINAL timing:
* common offset 0.02/0.05/0.1/0.2 slot: 4.26/9.86/13.55/15.07 at 20 dB, 11.3/20.5/24.1/25.3 at 30 dB.
  **Very sensitive:** 0.02 slot = 83 ns.
* random edge jitter rms 0.01/0.02/0.05 slot (new every symbol): 2.57/5.44/9.78 at 20 dB, 9.9/14.3/19.7 at 30 dB.
* **Blind offset recovery** (scan the offset, keep the smallest pilot residual; no truth):
  true 0.05/0.20/0.40 -> estimated 0.050/0.200/0.400; loss blind vs known 0.18/0.19, 0.24/0.27,
  0.15/0.16 at 20 dB (0.28/0.27, 0.41/0.39, 0.35/0.36 at 30 dB). Bar (0.2 slot within 0.3 dB of
  known): PASS. Caveat: these offsets lay ON the scan grid. Off-grid offsets with golden-section
  refinement, and finer jitter levels (0.1-0.5% of a slot), are running in `week6_timing_offgrid.m`.
* Conclusion so far: constant clock offsets are recovered blindly; jitter cannot be estimated
  away and becomes a hardware requirement (well below ~0.5% of a slot, ~20 ns).

### (5b) Hardware follow-up (`week6_hardware_fix.m`)
* Switching transients **fixed**: kernels from a calibrated or blindly estimated rise time.
  Loss (20/30 dB), nominal -> calibrated -> blind: rise 0.05: 0.75/4.78 -> 0.30/0.41 -> 0.43/0.42;
  0.1: 1.96/10.2 -> 0.30/0.61 -> 0.26/0.61; 0.2: 3.13/11.5 -> 0.23/0.50 -> 0.23/0.50;
  0.3: 3.83/11.7 -> 0.23/0.36 -> 0.21/0.37. Rise time estimated within 0.03 (exact from 0.1).
* Phase noise with a continuous phase and CPE removal still shows 2.2/8.1 dB at 100 Hz, but that
  comparison was UNFAIR: its reference knew each symbol's random phase-noise ICI. Rerun with a fair
  reference (true channel, unknown phase noise) in `week6_phasenoise.m`.

### (5c) Phase noise, fair comparison (`week6_phasenoise.m`)
Continuous Wiener phase over pilots and data, common phase removed per symbol. Rate (bit/s/Hz) at
20 / 30 dB, no phase noise 6.63 / 9.94:
| linewidth | genie channel (knows H0, not the phase noise) | ours | ours vs genie channel |
|---|---|---|---|
| 10 Hz | 6.57 / 9.44 | 6.45 / 9.30 | 0.37 / 0.42 dB |
| 100 Hz | 6.31 / 8.31 | 6.15 / 7.99 | 0.52 / 0.99 dB |
| 1 kHz | 4.72 / 5.24 | 4.45 / 4.79 | 0.90 / 1.37 dB |
| 10 kHz | 1.96 / 2.00 | 1.59 / 1.54 | 1.20 / 1.38 dB |
Most of the loss is the generic OFDM phase-noise ICI that any receiver suffers. Ours adds 0.4-1.4 dB
because the pilots are corrupted too. Bar (within 0.5 dB of genie channel at <= 1 kHz, 20 dB):
**FAIL** at 1 kHz (0.90 dB); met at 10 Hz, borderline at 100 Hz. Requirement for the paper: equivalent
linewidth below ~100 Hz at 30 kHz spacing, or standard phase tracking (NR PT-RS). A receiver
knowing the phase noise itself loses nothing (phase noise is unitary), but that bound is unreachable.

### (4) 3GPP TDL channels, larger surfaces, near field (`week6_channels.m`, `tdl38901.m`)
Gap to perfect (bit/s/Hz) at 20/30 dB, rank 8, 2 pilots; first number with the 0.3 us window, second
with a CP-wide (2.34 us) window:
* TDL-A: 30 ns 0.08/0.19 (CP 0.72/1.40); **100 ns 0.54/2.49; 300 ns 2.13/5.39** (CP 1.27/4.39)
* TDL-C: 30 ns 0.07/0.26; **100 ns 0.87/3.30; 300 ns 1.66/4.66** (CP 1.12/3.92)
* TDL-D: 30 ns 0.14/0.79; **100 ns 1.39/3.83; 300 ns 1.64/4.45** (CP 1.67/3.67)
* Q = 64/256/1024 (TDL-C, 30 ns): 0.06-0.11 / 0.17-0.25. Near field (Q = 1024, TDL-D LoS from exact
  distances, user at 3/10/30 m): 0.09-0.18 / 0.38-0.60. Fine.
**Weakness found:** standard 3GPP channels with realistic delay spread (100-300 ns) break the default
estimator. They have 14-24 taps from different directions (our own model had 6 paths), so rank 8
and a 0.3 us window are too small. **Correction:** the earlier "robust to 100 ns" result (week3_delayspread)
used our 6-path model and was too optimistic for standard channels.
Fix under test (`week6_channels_fix.m`): window = min(CP, 10 x DS), rank 8/16/32, 2/4/8 pilots.

### (2b) Timing follow-up (`week6_timing_offgrid.m`)
* Blind offset recovery with OFF-grid offsets (uniform in +/-0.5 slot; coarse scan then golden-section
  refinement of the pilot residual): offset error rms 0.0033 slot (~14 ns). Loss blind/known:
  0.20/0.19 dB at 20 dB (**bar met**), but **3.69/0.45 dB at 30 dB**: the residual 14 ns error is
  too large at high SNR. Open issue: needs a finer timing estimate (final refinement with the full
  estimator, or averaging over more pilots).
* Jitter requirement, loss at 20/30 dB: 4 ns rms (0.1% slot) 0.24/0.94; 8 ns 0.40/2.05;
  21 ns 0.97/5.91. **Requirement: switching jitter below ~20 ns for 20 dB operation and ~4 ns for
  30 dB.** FPGA controllers typically achieve sub-ns jitter, but it must be stated.

### (2b, 2c) Strongest honest baselines (`week7_strong_baselines.m`), 20 trials, rate at 20/30 dB
**2b, BEM with more harmonics and pilots** (same LMMSE/DPSS prior). Np = 2 / 8 / 32:
* |h| <= 12: 2.13/2.26/2.33 | 1.09/1.26/1.41
* |h| <= 40: 2.06/4.20/4.44 | 1.82/3.34/3.82
* |h| <= 128: 0.15/4.25/6.08 | 0.07/6.06/8.15
* ALL |h| <= 255 (4 trials only, ~8700 unknowns): 0.13/2.97/6.18 | 0.08/4.19/9.37
* ours, Np = 2: 6.51 | 9.69; perfect 6.59 | 9.89
**Correction:** the earlier statement "BEM does not improve with more pilots" holds only for the
truncated BEM. Full-harmonic BEM with 32 pilot symbols reaches 6.18 / 9.37, close to ours.
**Correct claim:** ours with 2 pilots beats full-harmonic BEM with 32 (a 16x pilot saving); at equal
overhead (Np = 2) full BEM collapses to 0.13.
**2c, ICI equalisers from the high-mobility OFDM literature**, on our estimate (oracle H in brackets):
banded MMSE (Schniter/Rugini-style) window 4 / 16 / 32: 3.23 (3.24) / 4.82 (4.83) / 5.45 (5.47) at 20 dB,
3.35 / 5.26 / 6.20 at 30 dB: limited by the harmonic tail **even with the true channel**.
LSQR/CG (Hrycak-style) 10 iterations 6.51 (6.59) | 9.66 (9.83); full MMSE 6.51 (6.59) | 9.69 (9.89).
The iterative matrix-free equaliser matches full MMSE in 10 iterations.

### (3) Precise blind timing (`week7_timing_refine.m`), 12 trials, off-grid offsets
Offset error rms (slot) coarse / two-stage / + data-aided, and loss (dB) incl. known timing:
* 20 dB: 0.0045 / 0.0043 / 0.0029 (~12 ns); loss 0.57 / 0.52 / **0.20** (known 0.25). PASS.
* 30 dB: 0.021 / 0.018 / 0.019 (~80 ns); loss 4.80 / 4.32 / 3.78 (known 0.41). **Bar (within 0.5 dB
  of known at 30 dB) FAILED.**
Diagnosis: the offset error is WORSE at 30 dB than at 20 dB, a near-ambiguity rather than noise. A timing
shift multiplies the channel by phase ramps that look like extra propagation delay, so the
delay-domain channel factors partly absorb a wrong offset and the pilot-residual surface is flat
near the truth. **Limitation:** blind timing recovery is adequate up to ~20 dB; at 30 dB the RIS and
receiver need timing alignment within ~10 ns (shared clock or calibration).

### (4b) Rich 3GPP channels: window, rank, pilots (`week6_channels_fix.m`), 6 trials
Window = min(CP, 10 x DS) (1.0 us at 100 ns, 2.34 us at 300 ns). Gap (bit/s/Hz) at 20 dB (30 dB),
Np = 2 / 4 / 8; rank 8, 16 or 32 makes little difference (ranges over rank shown):
* TDL-A 100 ns: 0.69-0.76 (2.5-2.8) / 0.21-0.32 (0.8-0.9) / **0.08-0.09 (0.24-0.27)**
* TDL-C 100 ns: 0.47-0.50 (1.7-1.9) / 0.17-0.24 (0.6-0.8) / **0.09-0.10 (0.18-0.26)**
* TDL-D 100 ns: 0.65-0.69 (2.1-2.3) / 0.31-0.37 (1.3-1.4) / **0.23-0.24 (0.85-0.95)**
* 300 ns (A/C/D): 1.1-1.5 (3.4-4.7) / 0.6-1.1 (2.6-3.3) / 0.35-0.92 (1.4-2.3)
Bar (< 0.3 at 20 dB for every model at 100 ns): **PASS with 8 pilot symbols**. 300 ns: open limitation.
The pilot count, not the rank, is the lever: a delay window matched to the spread means more unknowns.
Overhead: 8 pilots per coherence block, e.g. ~3% walking (>300 symbols), ~20% urban (~40 symbols);
comb pilots with empty gaps (item 11) let these spread over time.

## Week 8: attacking the five remaining limitations (pre-registered 2026-10-21)

Bars fixed BEFORE running. A failed bar stays failed and is reported.

| # | Limitation | Idea | Bar (pass if ...) |
|---|---|---|---|
| L5 | Fast equaliser only shown by op count | time matrix-free CG vs dense MMSE in MATLAB | CG is faster in wall-clock at M = 3276 with the same SER (within CI) |
| L4 | Full-band update cost extrapolated | measure cold vs warm-start updates at M = 3276, moving user | warm start >= 5x faster in wall-clock and SER within 1.2x of cold |
| L2 | Timing must be known to ~10 ns at 30 dB | a clock offset d equals delaying the output hop and advancing the input hop by d slots (D T D* o (FU)(FV)^T = T o (D FU)(D* FV)^T), so a two-sided delay window with a guard absorbs it | unknown off-grid offset up to +-0.02 slot (83 ns), NO timing search: loss within 1 dB of known timing at 30 dB |
| L2b (added 2026-10-21 after L2 failed, before running) | same | choose the offset by the Bayesian evidence of the linear first-stage model (penalises fits that need improbable spread factors), instead of the pilot residual, which is nearly flat in the offset | off-grid offsets uniform in +-0.5 slot (as week7), 30 dB: evidence-blind within 1 dB of known timing |
| L1 | 300 ns delay spread unsolved | learn the static input-side factor V (BS -> RIS hop, full inner dimension r = Q) once from several frames; per frame estimate only U (linear, convex) | TDL-A/C/D at 300 ns, 20 dB, 8 pilots: gap <= 0.3 bit/s/Hz for all three; and lower than the single-frame estimator at 30 dB |
| L1b (added 2026-10-21 after the first L1 model showed learned-W = single-frame, before running) | same | learn W from the row space of ALL frames' stacked linear estimates (SVD), then ALS, with the U prior covariance also learned over frames; Nf = 10 and 30 frames (40 cut to 30 for run time, before running) | same bar as L1: 300 ns, 20 dB, 8 pilots, gap <= 0.3 for A, C, D, and below single at 30 dB |
| L3b (added 2026-10-21 after L3 failed, before running) | same | proper Kalman filter on vec(U): covariance propagated, process noise P0w kron(Ru, diag lam) with the same physics-set P0w, W and Ru from learn_input_side, measurement noise per symbol = max(s2, residual power of that symbol with the predicted channel and its decisions) | same bar as L3 |
| L1c (added 2026-10-21, before running) | 300 ns | L1b learning with Nf = 100 frames (stacked-SVD init from the first 30, then 6 ALS iterations over all 100); the BS -> RIS hop is static for minutes, i.e. thousands of frames | same bar as L1: 300 ns, 20 dB, 8 pilots, gap <= 0.3 for A, C and D, and below single at 30 dB |
| L5b (added 2026-10-21, before running) | equaliser speed | batched shared-kernel operator (week8_eq_walltime2) as the method; CG stopped when the normal-equation residual falls below 1e-2 of its start (max 50 iterations), instead of a fixed count | at M = 3276 (12 symbols, 20 dB, 16-QAM): CG faster than dense in wall-clock, and SER no worse than 1.1x dense; same at M = 2048 reported |
| L3 | 120 km/h outdates the estimate in a few symbols | with V fixed, track U symbol by symbol from decisions (Kalman-style ridge to the previous U) | 120 km/h: gap <= 0.25 bit/s/Hz averaged over symbols 1-40 after the pilots at 20 dB, and <= 0.5 over 1-20 at 30 dB |


### L5 result: equaliser wall-clock (`week8_eq_walltime.m`, `week8_eq_walltime2.m`)
Pre-registered run (fast_coupling as used everywhere, 25 CG iterations, 12 symbols, same channel):
M = 512 / 1024 / 2048 / 3276: dense 0.09 / 0.33 / 1.49 / 5.56 s, CG 2.25 / 4.57 / 10.95 / 18.68 s;
SER identical (3276: 0.0652 vs 0.0648). **Bar FAILED: CG is 3.4x slower at M = 3276 in MATLAB.**
Post-hoc (labelled as such): with identical sequences all elements share one kernel per state, so the
operator can do one batched FFT per state (states summed in frequency, one inverse FFT). With this
operator: M = 3276 dense 4.52 s vs CG-10 2.88 s (1.6x faster, SER 0.0088 vs 0.0090) and CG-25 7.41 s;
M = 6552 (dense not run, 690 MB per matrix): CG-10 6.6 s, CG-25 14.3 s, but here CG-10 is no longer
converged (SER 0.153 vs 0.125). Honest statement for the paper: the crossover in MATLAB is near the
full n78 band and needs the batched operator and ~10 iterations; the robust advantage is memory
(O(QM) vs O(M^2)) and the operation count, not MATLAB seconds. The 6552 point also shows the CG
iteration count must grow with bandwidth.

### L2 result: guarded delay window (`week8_timing_guard.m`), 12 trials, loss dB at 20/30 dB
| case | plain | guard | known |
|---|---|---|---|
| (a) common offset <= 83 ns, guard 100 ns | 0.56 / 1.76 | 0.92 / 2.32 | 0.21 / 0.39 |
| (b) common offset <= 208 ns, guard 250 ns | 4.02 / 12.28 | 3.20 / 9.21 | 0.19 / 0.28 |
| (c) per-element skew <= 42 ns, guard 50 ns | 0.57 / 1.77 | 0.39 / 1.26 | 0 (gap to perfect H) |
**Bar FAILED** (a: guard is worse than doing nothing). The identity behind the idea is exact
(checked: 2.5e-15), and a diagnostic showed the shifted factors ARE representable even in the plain
basis (projection error 1e-4). So representation was never the problem: the DPSS prior shrinks the
improbable shifted factors, and a wider window spreads the prior over more dimensions (more noise).
The same fact explains week7's 30 dB failure: the pilot residual is nearly flat in the offset.
This led to L2b (evidence-based timing, bar fixed before running).

### L2b result: blind timing by Bayesian evidence (`week8_timing_evidence.m`), 12 trials
Same seed, channels and offsets as week7_timing_refine (its residual numbers are reproduced exactly).
Off-grid offsets uniform in +-0.5 slot, Np = 2:
| | offset error rms (slot) | loss resid / evid / known (dB) |
|---|---|---|
| 20 dB | resid 0.0045, evid 0.0011 | 0.57 / **0.25** / 0.25 |
| 30 dB | resid 0.0210, evid 0.0012 (5 ns) | 4.80 / **0.40** / 0.41 |
**Bar PASSED** (30 dB: evidence-blind 0.40 dB vs 0.41 dB known). The receiver no longer needs the RIS
timing offset: it finds it from 2 pilot symbols to ~5 ns. The 30 dB "~10 ns alignment" limitation is
removed. What remains: random per-edge jitter changes the pattern shape, not just its delay, so it
is not absorbed; it must stay below ~4 ns at 30 dB (cheap clocks have ps-level jitter). A fixed
per-element skew would need 16 offsets; not tested.

### L1 result: two-timescale estimator, first version (`week8_twoscale_300ns.m`, `est_twoscale.m`)
TDL-A 300 ns, 6 trials, gap (bit/s/Hz) for Np = 2/4/8:
| | single frame | learned W (10 frames) | true W (genie) |
|---|---|---|---|
| 20 dB | 1.53 / 0.77 / 0.33 | 1.38 / 0.77 / 0.37 | 0.57 / 0.29 / 0.11 |
| 30 dB | 4.54 / 3.42 / 1.91 | 4.41 / 3.03 / 1.51 | 1.89 / 0.84 / 0.40 |
**Bar FAILED** (TDL-A at 20 dB, 8 pilots: 0.37 > 0.3; the bar needs all three models), so the run was
stopped after TDL-A to free the machine (`results/week8_twoscale_300ns.json` transcribed from the log).
The genie column shows the idea is right: with the true static side, 300 ns is solved (0.11 at 20 dB,
0.40 at 30 dB). The weak link was learning W: it started from ONE frame's poor estimate. This led to
L1b (all-frame initialisation by SVD of the stacked linear estimates, plus a learned U prior), bar
unchanged and fixed before running. Development of L1b used one separate smoke-test channel
(TDL-A, 30 dB: single 1.80, stacked-SVD init 0.61, + learned U prior 0.48, genie 0.13).

### L3 result: decision-directed tracking at 120 km/h (`week8_tracking_120.m`), 8 trials
Average gap (bit/s/Hz) over symbols 1-40 at 20 dB / 1-20 at 30 dB:
reuse single-frame 0.33 / 1.13 | reuse two-timescale 0.28 / 1.11 | **track 0.52 / 1.30**.
**Bar FAILED**: tracking is worse than reusing the pilot estimate. Two flaws: (i) the update ignored the
uncertainty of the previous estimate (a ridge to the previous U, not a Kalman filter), so decision
errors (16/64-QAM) bias it and are never corrected; (ii) it used the first, weak learning of W.
This led to L3b (proper Kalman filter), bar unchanged and fixed before running.

### L3b result: Kalman tracking (`week8_tracking_120b.m`), same channels as L3, 8 trials
Average gap (bit/s/Hz) over symbols 1-40 at 20 dB / 1-20 at 30 dB:
reuse single-frame 0.33 / 1.13 | reuse two-timescale (L1b learning) **0.24** / 1.02 | Kalman 0.46 / 0.95.
**Bar FAILED.** Gap at symbols 1/10/20/40, 20 dB: Kalman 0.15 / 0.17 / 0.47 / 0.75 (decision-directed
drift after ~10 symbols); 30 dB: even symbol 1 is at 0.44-0.52 for every receiver, because the
intra-symbol Doppler ICI at 120 km/h (-34 dB, week5) is not in our model, so <= 0.5 at 30 dB is out of
reach without modelling it. Stopped after two principled attempts. Remaining limitation, stated as
such: at 120 km/h pilots must recur every ~10 symbols at 20 dB (as 5G's additional DMRS do), and 30 dB
needs a Doppler-ICI-aware model. Side note (not the pre-registered receiver, not counted as a pass):
simply reusing the two-timescale estimate gives 0.24 at 20 dB.

### L4 result: warm start at the full band (`week8_fullband_warm.m`), M = 3276, 4 trials x 5 updates
30 km/h, refresh every 40 symbols, Np = 2, 16-QAM, 20 dB, data 20 symbols after the pilots:
time per update cold 12.5 s, **warm 1.4 s (8.6x faster)**; SER cold 0.3976, warm 0.3986 (ratio 1.00).
**Bar PASSED** (>= 5x faster, SER within 1.2x). The ~2.5 TOPS extrapolation of week 7 is now a
measurement: ~9x fewer seconds per update at the full band.
Caveat found here (applies to week5 too): the absolute SER (0.40 vs 0.029 with the current true channel)
is set by staleness, not by the estimator. Diagnostic (1 trial): data 2 symbols after the pilots, SER cold
0.056 / true 0.043 (1.3x, as week4); 20 symbols after, cold 0.31, and even the TRUE channel of the pilot
time gives 0.32. At 30 km/h the dominant path rotates by up to ~0.45 rad in 20 symbols. The rate metric
of week5 (MMSE SINR with |diag(WH)|^2 as signal) is blind to a common phase rotation, i.e. it assumes
the receiver tracks the common phase (as 5G phase-tracking RS do). Week5's "~40 symbols at 30 km/h"
holds with common-phase tracking; an uncoded detector without it needs fresh pilots far more often.

### L1b/L1c results: two-timescale learning at 300 ns (`week8_twoscale_300ns_b.m`, `_c.m`)
L1b (stopped after TDL-A, verdict already fixed): TDL-A, 8 pilots, gap at 20/30 dB: single 0.40/1.84,
10 frames 0.36/1.08, 30 frames 0.31/0.77, genie 0.15/0.38. **L1b bar FAILED** (0.31 > 0.3).
L1c (100 frames, 6 trials), gap at 20 dB / 30 dB, 8 pilots:
| | single frame | learned W, 100 frames | true W (genie) |
|---|---|---|---|
| TDL-A 300 ns | 0.36 / 1.91 | **0.23** / 0.85 | 0.12 / 0.37 |
| TDL-C 300 ns | 0.46 / 1.84 | **0.26** / 0.99 | 0.12 / 0.32 |
| TDL-D 300 ns | 0.85 / 2.74 | 0.64 / 2.89 | 0.56 / 1.61 |
**L1c bar FAILED on TDL-D** (A and C pass at 20 dB and halve the 30 dB gap). TDL-D fails even with
the TRUE static side (0.56), so it is not the learning: at 300 ns TDL-D has taps up to 12.5 x 300 ns
= 3.8 us, beyond the 2.34 us normal CP and outside any CP-limited delay window (taps beyond the CP
carry ~ -20 dB). A real system would use an extended CP there. Statement for the paper: with the
static BS -> RIS hop learned over ~100 frames (well within its coherence of minutes), 300 ns delay
spread is handled for TDL-A and TDL-C with 8 pilots; channels whose taps exceed the CP remain out of
scope. The learning principle itself is known for conventional RIS (RESEARCH_NOTES, fourth sweep).

### L5b result: batched operator + stopping rule (`week8_eq_walltime3.m`), fresh seed
| M | dense | CG (iterations used) | SER dense / CG |
|---|---|---|---|
| 2048 | 0.88 s | 0.83 s (8) | 0.0128 / 0.0151 |
| 3276 | 2.95 s | **1.12 s (5)** | 0.0001 / 0.0001 |
| 6552 | not run (memory) | 1.85 s (5) | - / 0.0001 |
**Bar PASSED** at M = 3276 (2.6x faster, same SER). Caveats: at 3276 both SERs are ~1e-4, so the SER
check is weak; at 2048 (reported, not in the bar) CG's SER is 1.18x dense, so the 1e-2 stopping rule
is slightly loose. With the batched operator CG time grows ~linearly in M (1.1 s -> 1.9 s from 3276
to 6552) while dense grows ~M^3 and does not fit at 6552.

### Pilot overhead: net throughput vs pilots, refresh interval, speed, channel (`week8_pilot_overhead.m`, Fig. 19)
Net rate = (N - Np)/N x mean rate over the data symbols of a refresh block of N symbols with Np pilot
symbols in front; N and Np optimised per case (N in 4..320, Np in 1..8 ours, 2..32 CE-BEM, CE-BEM's
best harmonic range per case). Channels: 6-path Rician (30 ns, 8 trials) and 3GPP TDL-A/C/D at 100 ns
(per-tap Doppler, 4 trials). CE-BEM on TDL limited to |h| <= 40 (|h| <= 128 needs ~1 GB per estimate).
Three ways the receiver handles the stale estimate: per-subcarrier genie (optimistic), ONE common phase/
gain corrected per symbol as 5G phase-tracking RS do (realistic, main numbers), no correction
(pessimistic). The first run used only the genie and the no-correction metrics and mislabelled the
genie as "phase tracking"; it was rerun with the realistic metric (same seed; the old metrics reproduce
exactly). First run kept as `results/week8_pilot_overhead_run1.json`.

| channel | speed | SNR | ours, common-phase tracking: net (Np/N, overhead) | ours: no tracking ... genie | CE-BEM best (same tracking) | one-tap | perfect |
|---|---|---|---|---|---|---|---|
| 6-path Rician, 30 ns | 3 | 20 | 6.28 (2/120, 2%) | 5.98 ... 6.33 | 5.14 (32/320, 10%) | 0.38 (8/320, 2%) | 6.49 |
| 6-path Rician, 30 ns | 3 | 30 | 9.29 (2/80, 2%) | 8.43 ... 9.47 | 6.87 (16/240, 7%) | 0.42 (8/320, 2%) | 9.79 |
| 6-path Rician, 30 ns | 30 | 20 | 6.13 (1/40, 2%) | 5.38 ... 6.53 | 4.83 (16/240, 7%) | 0.39 (8/320, 2%) | 6.60 |
| 6-path Rician, 30 ns | 30 | 30 | 8.93 (1/30, 3%) | 7.21 ... 9.40 | 6.02 (16/240, 7%) | 0.42 (8/320, 2%) | 9.90 |
| 6-path Rician, 30 ns | 120 | 20 | 5.59 (1/160, 1%) | 3.70 ... 6.23 | 3.43 (8/320, 2%) | 0.40 (8/320, 2%) | 6.49 |
| 6-path Rician, 30 ns | 120 | 30 | 7.55 (1/10, 10%) | 4.40 ... 8.56 | 2.96 (16/320, 5%) | 0.43 (8/320, 2%) | 9.78 |
| TDL-A, 100 ns | 3 | 20 | 5.05 (4/60, 7%) | 4.91 ... 5.34 | 3.44 (16/160, 10%) | 0.33 (8/240, 3%) | 5.91 |
| TDL-A, 100 ns | 3 | 30 | 6.61 (4/40, 10%) | 6.50 ... 7.03 | 3.25 (16/120, 13%) | 0.39 (8/240, 3%) | 8.94 |
| TDL-A, 100 ns | 30 | 20 | 4.70 (1/10, 10%) | 4.47 ... 5.29 | 2.25 (8/40, 20%) | 0.34 (8/80, 10%) | 6.31 |
| TDL-A, 100 ns | 30 | 30 | 5.89 (2/15, 13%) | 5.41 ... 6.96 | 1.93 (16/60, 27%) | 0.40 (8/80, 10%) | 9.54 |
| TDL-A, 100 ns | 120 | 20 | 4.26 (1/8, 12%) | 3.65 ... 4.97 | 1.53 (16/240, 7%) | 0.34 (8/320, 2%) | 6.28 |
| TDL-A, 100 ns | 120 | 30 | 5.09 (1/6, 17%) | 4.34 ... 6.08 | 1.24 (16/240, 7%) | 0.38 (8/320, 2%) | 9.51 |
| TDL-C, 100 ns | 3 | 20 | 5.71 (4/80, 5%) | 5.55 ... 6.02 | 3.70 (16/240, 7%) | 0.38 (8/320, 2%) | 6.34 |
| TDL-C, 100 ns | 3 | 30 | 8.05 (4/60, 7%) | 7.61 ... 8.60 | 3.33 (16/240, 7%) | 0.42 (8/320, 2%) | 9.58 |
| TDL-C, 100 ns | 30 | 20 | 4.81 (1/10, 10%) | 4.69 ... 5.20 | 2.10 (8/40, 20%) | 0.33 (8/240, 3%) | 5.95 |
| TDL-C, 100 ns | 30 | 30 | 6.03 (1/8, 12%) | 5.80 ... 6.73 | 1.73 (8/30, 27%) | 0.37 (8/240, 3%) | 9.13 |
| TDL-C, 100 ns | 120 | 20 | 4.34 (1/8, 12%) | 3.65 ... 4.95 | 1.66 (16/320, 5%) | 0.33 (8/320, 2%) | 6.07 |
| TDL-C, 100 ns | 120 | 30 | 5.37 (1/6, 17%) | 4.39 ... 6.39 | 1.39 (16/320, 5%) | 0.40 (8/320, 2%) | 9.30 |
| TDL-D, 100 ns | 3 | 20 | 5.58 (4/80, 5%) | 5.22 ... 5.77 | 3.63 (16/160, 10%) | 0.39 (8/240, 3%) | 6.25 |
| TDL-D, 100 ns | 3 | 30 | 8.06 (2/60, 3%) | 7.46 ... 8.51 | 3.26 (16/160, 10%) | 0.40 (8/240, 3%) | 9.51 |
| TDL-D, 100 ns | 30 | 20 | 5.06 (1/20, 5%) | 4.53 ... 5.60 | 2.88 (32/320, 10%) | 0.37 (8/320, 2%) | 6.49 |
| TDL-D, 100 ns | 30 | 30 | 7.06 (1/15, 7%) | 5.96 ... 7.73 | 2.55 (16/320, 5%) | 0.42 (8/320, 2%) | 9.75 |
| TDL-D, 100 ns | 120 | 20 | 3.91 (1/8, 12%) | 3.10 ... 4.74 | 1.64 (16/160, 10%) | 0.32 (8/240, 3%) | 6.13 |
| TDL-D, 100 ns | 120 | 30 | 4.94 (1/8, 12%) | 3.66 ... 5.51 | 1.36 (16/160, 10%) | 0.34 (8/160, 5%) | 9.32 |

Findings:
* Ours beats the best CE-BEM after overhead in all 24 cases, by 1.1-4.8 bit/s/Hz (median 2.7); one-tap
  stays at ~0.4.
* Overhead needed (ours, realistic tracking, 20 dB): walking 2-7%; 30 km/h 2% (6-path) to 5-10% (TDL);
  120 km/h 12% on TDL (1 pilot every 8 symbols). The 6-path model needs less at speed because its
  dominant LoS path makes the Doppler mostly a common phase, which tracking removes; the NLoS TDL
  numbers are the realistic ones.
* Longer delay spread needs more pilots per refresh only when the channel is slow (TDL walking: 4
  pilots); when it moves, 1 pilot refreshed often beats many pilots refreshed rarely.
* Without any phase tracking the 120 km/h overhead rises to ~25% (1 pilot every 4 symbols).
* Replaces the earlier rough estimate ("~3% walking, ~20% urban"): walking 2-7%, 30 km/h 2-10%.

## Run
```
matlab -batch "cd('C:\Users\hp\RIS\tm_ris_ofdm\matlab'); gate1_reproduce"
matlab -batch "cd('C:\Users\hp\RIS\tm_ris_ofdm\matlab'); gate2_equalization"
```

## Plan
* Week 2: done (Gate 3 above).
* Week 3: done (results above).
* Remaining before writing: full-band scale (M = 3276) run using the fast operator;
  more trials for (d); robustness to element-model mismatch (receiver's Table I
  parameters perturbed by +/-5%, since the physics-aware estimator depends on them);
  coded BER (optional).
* Weeks 4–5: letter draft.
