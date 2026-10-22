# Pre-registration log

Each entry is written BEFORE the corresponding run. Format: date | git hash | experiment | frozen
hyperparameters | criterion. Earlier pre-registrations (weeks 6-8) are in README.md.

## 2026-10-22 | git 75f31e0 | Step 1: workstreams A2 + B1-B4 (physical timing, bulk delays, prior sensitivity)
Frozen estimator (unchanged from the paper): S1, K=8, Q=16, M=256, DPSS tau_max=0.3 us, L=16, r=8, 20 ALS
iterations, ridge 1e-3; evidence search grid -0.5:0.025:0.5 slot + 14 golden steps; Np=2.
Channels: shared-path model (1 LoS + 5 NLoS per hop, 30 ns rms, K=10 dB), 30 independent draws per point,
seeds shared across receivers, paired differences with t-based 95% CIs.

* T1 (physics check, no criterion on performance): with per-hop bulk delays tau_b ~ U(50,300) ns,
  tau_u ~ U(30,300) ns and receiver sync to the first cascaded arrival,
  (i) RIS locked to the incident signal  ==> synced H equals the no-bulk H (rel. error < 1e-10);
  (ii) RIS locked to the BS clock        ==> synced H equals the no-bulk H with clock offset d = -tau_b/T_slot.
  If either identity fails, the timing claims are re-examined before anything else.
* T2 (A2/B1, performance), 20 and 30 dB, receivers: nominal timing (no search), evidence search (Alg. 2),
  known offset, perfect CSI.
  PASS mode (i): nominal and evidence each within 0.3 bit/s/Hz of the matched (no-bulk) result.
  PASS mode (ii): evidence within 0.3 bit/s/Hz of the known-offset receiver.
  Mode (i) also with switching jitter 10 ps and 1 ns rms (per edge, fixed over the frame).
* T3 (B2/B4, single channel + 30 channels): curves of in-sample pilot residual, HELD-OUT residual (2 extra
  pilot symbols not used for estimation) and evidence vs d over +-0.5 slot; projection error of the
  shifted responses vs |offset| up to 2 us for L=16. Hypothesis (B4): at a wrong d the in-sample residual
  is near its minimum but the held-out residual is not. PASS if, over 30 channels, the held-out residual
  minimiser is within 0.005 slot of the evidence maximiser in >= 90% of channels.
* T4 (B3, prior sensitivity), 30 dB, mode (ii): DPSS window shifted to [t0, t0+tau_max] with
  t0 in {-100, +100} ns, tau_max scaled 0.5x and 1.5x. Report rate loss and d_hat bias/rms.
  PASS if the rate loss vs known offset stays within 0.5 dB for all four variants.

### Outcomes (recorded 2026-10-22 after the runs)
T1 PASS (6.9e-15). T2 PASS (both modes, both SNRs). T3 PASS (93%). T4 FAIL (t0 = +100 ns: +6.68 bit/s/Hz loss);
other three variants pass. Details: SUMMARY_B.md.

## 2026-10-22 | git 6d95ea3 | Step 1b: corrections of C-C (fresh seeds, decided before running)
Frozen estimator as in step 1. Mode (ii) physical timing (RIS locked to the BS clock), 30 NEW channels.

* F1 (fix for T4: prior window vs sync reference). A window starting after the sync reference is equivalent to
  the receiver syncing late. Model a receiver sync error delta_s ~ U(-150, +150) ns (first-arrival detection
  error). Receivers, all with the evidence search:
    (a) default window [0, 0.3 us]                     (paper as is)
    (b) margin window [-150 ns, 0.3 us]                (margin = assumed max sync error; L = 16)
    (c) evidence-selected window start t0 in {-150, -75, 0, +75} ns, jointly with d (max evidence)
  Reference: perfect sync, known offset, default window (paired, same noise).
  PASS for a variant: mean paired loss vs reference <= 0.165 bit/s/Hz (0.5 dB) at 30 dB AND upper t-CI bound
  <= 0.33 bit/s/Hz. 20 dB reported too. The paper adopts a passing variant; if none passes, T4 stays a limitation.
* F2 (where does the evidence beat the residual?). Paired gap difference residual - evidence, t-CI, in four
  regimes, ALL reported: (R1) week-7 procedure (8-iteration estimator, 0.025 grid, refine +-0.025) vs evidence
  on the same grid; (R2) fine-grid residual (0.01, 20 iterations) vs evidence, Np = 1; (R3) same with TDL-C
  100 ns, window 1 us, r = 16, Np = 8; (R4) wall-clock per search, fine-grid residual vs evidence (M = 256).
  The paper claims an evidence advantage only in regimes whose CI excludes zero; otherwise "equivalent, cheaper".

## 2026-10-22 | git 5ac0dde | Step 2: workstream A1 + A3 (independent ground truth, model mismatch)
Frozen estimator unchanged (tau_max 0.3 us, L 16, r 8, 20 it, Np 2, pilots only). Shared-path channel, 30 channels,
20/30 dB, paired t-CIs against the matched result (same channel, nominal truth, same noise).
* A1a: received pilots/data from the sampled-gating time-domain simulator (td_simulate, Nfft = 4096) at
  Q = 16, K = 8, multipath on both hops; true H for scoring obtained by probing the simulator column by column.
* A1b: element memory: add the exact switching transient of the single-pole element (state-space solution:
  at each edge the internal state carries over and decays with the new state's pole). Report its relative
  size and the loss.
* A3 (truth built per slot, H = sum_k T_k o (Aout_k Ain_k^T), validated against coupling_matrix):
  (m1) state-dependent output response: a_d,q -> a_d,q (1 + eps eta_{q,s}), eps in {0.05, 0.1, 0.2};
  (m2) mutual coupling depending on BOTH elements' states, kappa in {0.02, 0.05, 0.08}, for identical
       patterns (G = 1) and for 4 random per-element patterns (G = 4);
  (m3) frequency-dependent array phase (beam squint, k = 2 pi f / c), Q = 16 and Q = 256;
  (m4) per-element resonance spread 10 MHz (kept from week 6) for completeness.
  Predictions written before running: state-independent coupling is absorbed exactly; joint-state coupling is
  absorbed for G = 1 and not for G = 4.
PASS (per mismatch at its middle level, and A1a/A1b): loss vs matched within 0.3 bit/s/Hz at 20 and 30 dB
(mean paired difference). Every level is reported; failures are diagnosed, not removed.
