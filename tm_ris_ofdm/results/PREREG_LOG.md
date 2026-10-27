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

## 2026-10-22 | git 6c66544 | Step 3: workstream C1 + C4 + E2 + E3 + E6 (fair baselines, re-derived pilot counts)
Channels: shared-path 30 ns (estimator window 0.3 us, L 16, r 8; CE-BEM two-hop window 0.6 us, L2 17) and
3GPP TDL-C 100 ns (cascade 1.73 us < CP; ours window 1 us, L 18, r 16; CE-BEM window 2 us, L2 28).
30 channels, 20/30 dB, Np in {2, 4, 8, 16, 32}, pilots only, frozen hyperparameters (listed in the script).
Estimators: one-tap; ideal-element model-aware LMMSE (strongest "ideal +-1 element" estimator, replaces the
25-coefficient Toeplitz straw man); CE-BEM |h| <= 12, 40, 128 and ALL harmonics (matrix-free CG LMMSE,
identical estimator to the dense one; all-harmonic at Np >= 8); OMP (harmonic x delay); ours (est_general).
Equalisers for EVERY estimator: naive MMSE; robust MMSE with sigma_e^2 from the pilot residual,
W = (Hh'Hh + (s2 + se2) I) \ Hh', se2 = max(||Y - Hh S||^2/(M Np) - s2, 0); robust with the ORACLE error power
(upper bound for any scalar regularisation). Rates scored on the true H, with the detector knowing G_mm
(metric 13) and with the detector using its own estimated gain diag(W Hh) (E6).
* C1 sanity (before any baseline claim): with the oracle-robust equaliser CE-BEM's rate must not fall from 20 to
  30 dB; if it does, the baseline is debugged before reporting.
* E3: minimal Np in the grid with mean gap <= 0.3 bit/s/Hz at 20 dB (robust-residual equaliser), per method and
  channel; "> 32" if never reached.
* E2: pilot ratio Np_min(best CE-BEM) / Np_min(ours).
* C4 / claim: paired difference ours - best baseline (best over harmonic sets and equalisers, i.e. favouring the
  baseline) per Np and SNR, t-CI. The claim "beats CE-BEM" survives only where the CI excludes zero.

## 2026-10-22 | git cdbc657 | Step 4: workstreams C2 (+G3), C3, C5
30 channels per point, 20/30 dB, pilots only unless stated, paired t-CIs. Frozen hyperparameters as in step 3.
* C2 ablation (shared-path 30 ns and TDL-C 100 ns; Np in {2, 4, 8}): (a) linear first stage only; (b) ALS low
  rank r in {2, 4, 6, 8, 16 (= Q)} WITHOUT the fallback; (c) full pipeline (ALS r = 8 / 16 with fallback);
  (d) full + one DD pass (6 data symbols, 16-QAM at 20 dB, 64-QAM at 30 dB). G3: fraction of channels where the
  fallback is chosen; pilot residual per ALS iteration; 5 random restarts at r = 8 (spread of final gap).
  No pass/fail: the ablation reports what each stage buys. If (a) gets most of the gain, the paper says so.
* C3 parametric element fit: the receiver knows the single-pole family but starts from parameters off by
  N(0, 10 MHz) in resonance and N(0, 10%) in decay rates; it fits the 6 parameters by maximising the evidence
  of the model-aware linear LMMSE (fminsearch), then uses that model. Compared with timing-only, oracle model,
  and unfitted wrong model, on (t1) the exact single-pole truth and (t2) truth with 10 MHz element spread.
  Reported as is; if the parametric fit beats timing-only, the paper's "no model needed" is reframed as
  robustness (equal accuracy without assuming a model family), not accuracy.
* C5 finite rise time rho in {0.01, 0.05, 0.1, 0.3} slot (linear ramp, both truth and our calibrated kernels):
  energy outside |h| <= 12; banded MMSE (+-4, +-16, +-32) with the TRUE channel; CE-BEM |h| <= 12, 40 with Np = 8;
  ours (Np = 2). Reported as is; the "banded equalisers fail" claim is restated with its dependence on rho.

## 2026-10-22 | git 2d53670 | Step 5: workstreams F1, F2, I1-I3, E4 (decided before running)
* F1 (week11_fullband30): full band M = 3276, 30 fresh channels, same receivers as week4. Report paired SER
  ratios (ours/true) with clustered bootstrap 95% CIs (resample channels, 10^4 resamples). No pass bar; the
  paper quotes the ratio and its CI, pilots-only and +1 DD side by side.
* F2 (week11_coded_stats): 200 channels x 10 blocks = 2000 blocks per SNR in 6..14 dB. Crossing SNR of 10%
  BLER (log-linear interpolation) with clustered bootstrap CI, for ours, perfect CSI, CE-BEM (12, 40, 40 robust).
  The paper quotes the gap ours - perfect with its CI.
* I1-I3 (week11_motivation): reported as is. Pre-stated expectation: an optimised static RIS gives the
  legitimate user at least as much rate as a TM-RIS; the paper's motivation is then that time modulation is
  used for security and our receiver reduces its cost to the legitimate user. I2: secrecy rate with
  eavesdropper without / with the pattern.
* E4 (week11_overhead_heldout): (N, Np) chosen per method and case by maximum mean net rate on the 20-channel
  selection set; reported on the 20-channel held-out set: mean net rate, t-CI, ratio to perfect CSI, 20 and
  30 dB, common-phase metric (main), all CE-BEM with robust equaliser and no harmonic cap.

## 2026-10-22 | git 5272ea6 | Step 6: D2, D3, G1, G2, K2, K3, H2, H3 (decided before running; descriptive, no pass bars
except where stated)
* D2 week11_cdl_rank: CDL-A/C/D (cluster angles, 20 rays) vs TDL + random directions at 100 ns; ours Np 2/4/8.
  Paired CDL - TDL gap differences reported; the paper uses CDL as the main 3GPP model if it differs.
* D3: effective rank (99% / 99.9% energy) of the true element-domain matrix vs Q 16..1024, narrow / 98 MHz
  bandwidth, far field without and with squint, near field 3 m / 10 m. The low-rank claim (Cor. 5) is restated
  with the measured ranks; where rank approaches Q, the paper says the low-rank property is weak.
* G1/G2 week11_theory: effective dof vs observations; BMSE bound vs achieved NMSE for the linear first stage.
* K2 week11_rise_small: rise times 2-83 ns, nominal vs calibrated kernels.
* K3 (computed, results/k3_oob.json, analytic): ACLR of the switched reflection, ideal edges and rise times.
* H2/H3 week11_eq_amortised (idle machine): amortised dense vs CG wall-clock, break-even N; CG tolerance needed
  to match the dense SER at each M.

## 2026-10-22 | git 5ae54bb | Step 7: G3b, removing the pilot-fit fallback (decided before running)
The ablation (week10_ablation) showed on TDL-C that the fallback (keep the linear estimate if it fits the
pilots better) is chosen in 37-57% of channels at 30 dB and makes the gap worse (Np 4: 0.77 vs 0.59). Validation on
FRESH channels (new seeds): shared-path 30 ns, TDL-A 100 ns, TDL-C 100 ns; Np 2/4/8; 30 channels; ranks frozen
(8 shared-path, 16 TDL); estimator otherwise unchanged.
PASS (adopt "no fallback" in the paper's Algorithm 1): mean gap without fallback <= mean gap with fallback + 0.02
bit/s/Hz in every case, AND the paired t-CI excludes zero in favour of "no fallback" in at least one TDL case at
30 dB. Otherwise the fallback stays and the ablation finding is reported as such.

## 2026-10-22 — H2/H3 rerun (deviation note)
The first week11_eq_amortised run started while week11_overhead_heldout was still running (the idle check in
idle_run.sh failed silently). Its output is kept as week11_eq_amortised_loaded.json and is NOT used for wall-clock
claims (H2). Iteration counts, residuals and SER (H3) do not depend on load and are reported from it. The identical
script is rerun unchanged on an idle machine (idle_run2.sh) for H2. Git 6dbd32e.

## 2026-10-23 — E4 parallelisation (no change to the test)
week11_overhead_heldout seeds every (channel, set) separately (rng(11400 + 100*ic + 10*iset)), so channels can be run
in separate processes with identical draws. TDL-C (ic 3) and TDL-D (ic 4) are run by week11_overhead_heldout_part.m
into week11_overhead_heldout_TDLC100.json / _TDLD100.json; the sequential job is stopped (by PID) after it writes
TDL-A. Shared-path and TDL-A come from week11_overhead_heldout.json. Selection rule unchanged (figures/, Python).

## 2026-10-24 | git 1b0f4eb | Step 8: standard-format curves with intermediate references (week12_standard_curves)
**Motivation.** The paper's figures show ours next to perfect CSI and every baseline at an error floor, with nothing in between. That is hard for a reader to compare with the literature, and it invites the question whether the baselines are weak. This step adds the formats used by RIS and doubly-selective channel-estimation papers, plus references that fill the middle.

**Fixed settings.** All hyperparameters are frozen as in week10_fair_baselines:
- channels: shared path 30 ns and TDL-C 100 ns;
- M = 256 and Q = 16;
- identical patterns, K = 8, D = 0.75;
- robust MMSE with the error power estimated from the residual;
- 30 new channels per point, with new seeds (rng 12100 + 100*part + channel).

**New references:**
- **genie:** a model-aware LMMSE that knows the true element responses Γ_s(f) exactly. H = (T1 diag Γ1 + T2 diag Γ2) ∘ (F C Fᵀ).
- **mmδ:** the same estimator with both resonances shifted by δ = 5, 10 or 20 MHz, a mismatched element model.
- **CE-BEM, all harmonics, Np = 32:** 16 times the pilots of ours, as a reference.

**Parts:**
- **S1:** NMSE, robust-MMSE rate and uncoded 16-QAM BER versus SNR from 0 to 30 dB in 5 dB steps. Np = 2 on the shared path and Np = 4 on TDL-C. Uncoded BER uses 20 data symbols per channel, sliced after unbiased MMSE.
- **S2:** element non-ideality continuum, H_α = (1−α)H_ideal + αH_OpenRIS, with α in {0, 0.1, 0.25, 0.5, 0.75, 1}, at 20 and 30 dB. Shared path, Np = 2.
- **S3:** NMSE versus Np in {1, 2, 4, 8, 16, 32} at 20 dB, shared path.

**Pre-registered checks.** None of these change existing claims. Each is reported whichever way it goes:
- **(a) Baseline sanity:** at α = 0 the ideal-element LMMSE rate gap is within 0.05 bit/s/Hz of genie. If it fails, the baseline implementation is suspect, and that is reported.
- **(b) Ours versus genie:** ours is within 0.1 bit/s/Hz of genie at α = 1, Np = 2 and 20 dB.
- **(c) Crossover:** report the smallest α at which the ideal-element LMMSE gap exceeds 0.3 bit/s/Hz at 20 dB.
- **(d) Mismatch ordering:** at 20 dB the mismatched models order as genie ≤ mm5 ≤ mm10 ≤ mm20 in NMSE. Report the observed order.

## 2026-10-25 | git ffc462b | Step 9: learning-based baseline (C6) and low-SNR rank selection (week12_c6_lowsnr)

**C6, learned-prior LMMSE.** Training channels come from the same generator as the test channels: the physical OpenRIS element, the shared path and TDL-C, with separate seeds (rng 13000+). The steps are:
1. Form the PCA basis B of vec(H) from N_train channels, with N_train in {100, 300, 1000}.
2. Keep the components with eigenvalue at least 1e-6 of the largest.
3. Use a zero-mean Gaussian prior with the sample eigenvalues.
4. Estimate the channel as the LMMSE from the pilots: y = (Sᵀ ⊗ I) B c + n.

The baseline therefore learns the element dispersion from data. It is the Bayes-optimal linear estimator for the learned statistics, which an MSE-trained network approaches for near-Gaussian channels.

- **Test channels:** 30 new channels (rng 13300+).
- **Parts:** SNR 0–30 dB with Np = 2 on the shared path and Np = 4 on TDL-C; and Np in {1, 2, 4, 8, 16, 32} at 20 dB on the shared path.
- **Comparison:** ours on the same channels, with paired t-intervals.
- **Reporting:** whichever way it goes. If the learned LMMSE beats ours, that is reported, together with its cost: N_train recorded channels from the deployed surface.

**Low-SNR rank.** The receiver knows σ² (A6), so a rank chosen per SNR is legitimate when selected on separate channels.
- **Candidate ranks:** {1, 2, 4, 8} on the shared path and {2, 4, 8, 16} on TDL-C.
- **Grid:** SNR 0:5:30, with Np = 2 on the shared path and Np = 4 on TDL-C.
- **Selection rule A:** pick r*(SNR) by maximum mean rate on 20 selection channels (rng 13100+), and report on 30 held-out channels (rng 13200+). The frozen rank is 8 on the shared path and 16 on TDL-C.
- **Selection rule B (pilots only, per channel):** BIC = ‖Y − ĤS‖²/σ² + 3Lr·ln(M·Np). This is reported on the held-out channels.
- **PASS:** the adopted rule (A, or B if B is no worse than A within 0.01) improves on the frozen rank at 0–15 dB by at least 0.01 bit/s/Hz with the paired CI above zero at one or more SNR points, AND is never worse by more than 0.01 at any SNR. Otherwise the frozen rank stays and the low-SNR crossover stands as reported.

### Step 9 outcome (recorded before Step 10 runs)
- **Low-SNR rank: FAIL.** Rule A improves 0–15 dB (shared path +0.008–0.015, TDL-C +0.034–0.050, CIs above zero). It is worse by 0.034 at 20 dB on TDL-C, which breaks the "never worse by > 0.01" condition. Rule B (BIC) always picks rank 1 and is far worse at high SNR. The frozen rank stays.
- **C6: learned-prior LMMSE beats ours** at every SNR and Np, on both channels, from N_train = 100.
  - At 20 dB the learned LMMSE loses 0.006 against ours 0.070 on the shared path, and 0.074 against 0.285 on TDL-C.
  - The test families are low-dimensional: 6 PCA components on the shared path and 52 on TDL-C.

## 2026-10-25 | git afe82e1 | Step 10: C6 under distribution shift, and a self-trained prior (week12_c6b)
**Training sets.** N_train = 300 for each set, 30 test channels per point, SNR {0, 10, 20, 30}. Np = 2 on the shared path and 4 on TDL; ours frozen as before.

| Family | Training sets |
|---|---|
| Shared path | true H from phys30 (rms 30 ns, K = 10 dB); true H from phys100K3 (rms 100 ns, K = 3 dB); and the self-trained set |
| TDL (100 ns) | true H from TDL-C, TDL-A and TDL-D; and the self-trained set |

The **self-trained** set uses ours' estimates of the training channels (phys30 and TDL-C respectively) from 8 pilot symbols at 20 dB. It needs no true channels and no element model.

**Tests:**

| Test family | Matched learned prior | Shifted learned prior | Self-trained |
|---|---|---|---|
| phys30 | phys30 | phys100K3 is not used here | yes |
| phys100K3 | phys100K3 | phys30 | no |
| TDL-C | TDL-C | none | yes |
| TDL-A | TDL-A | TDL-C | no |
| TDL-D | TDL-D | TDL-C | no |

Ours runs on every test.

**Pre-registered readouts.** Paired t-intervals, reported whichever way they go.
- **(a) Shift:** for each shifted pair, is ours ≤ the shifted learned prior at 20 dB?
- **(b) Self-trained:** the self-trained prior beats ours at 20 dB with the paired CI above zero, on both phys30 and TDL-C. If so, it is ADOPTED as an optional receiver mode ("ours + self-trained prior"), with the cost stated: N_train channels measured with 8 pilot symbols.

## 2026-10-27 | git b889849 | Step 11: comparison in the format of Tang et al. 2007, Figs. 10–11 (week12_snr_tang)
**Purpose.** The paper's SNR figure spans 0–30 dB with 16-QAM. Tang et al. (CE-BEM, our main baseline) report channel MSE and BER from 10 to 40 dB with QPSK. This run reproduces that format so the two can be compared directly. No hypothesis is tested, and no claim depends on it; results are reported as they come.

**Settings.** Hyperparameters are frozen as in week12_standard_curves.
- **SNR:** 10:5:40 dB.
- **Modulation:** uncoded BER for QPSK and for 16-QAM after unbiased MMSE, 100 data symbols per channel.
- **Channels:** 30 new channels (rng 17100+).
  - shared path: Np = 2;
  - TDL-C 100 ns: Np = 4.
- **Estimators:** one-tap, ideal-element LMMSE, exact-element LMMSE, CE-BEM with abs(h) ≤ 128 at ours' Np, CE-BEM with all harmonics at Np = 32, OMP, and ours.
- **Reference:** perfect CSI.
