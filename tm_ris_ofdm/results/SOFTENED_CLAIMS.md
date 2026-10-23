# Claims softened, corrected or removed during the review (2026-10-22/23)

| Old wording | New wording / status | Evidence |
|---|---|---|
| "the first receiver for TM-RIS OFDM" | "to our knowledge, the first pilot-based receiver for OFDM through a dispersive TM-RIS"; closest works named | SUMMARY_J |
| "within 0.1 bit/s/Hz with two pilot symbols" (general) | two pilot symbols on the shared-path channel; four to eight on 3GPP TDL-C 100 ns (within 0.2 at 20 dB) | week10_fair_baselines |
| "16x fewer pilots than CE-BEM" | ours with 2 exceeds full-harmonic CE-BEM with 32; no baseline within 0.3 with <= 32 | week10_fair_baselines |
| Toeplitz ideal-element baseline (2.27/1.42) | replaced as strongest baseline by the ideal-element LMMSE (loss 0.4-0.5 @20, 2.0-2.3 @30 dB); old row labelled as the one-sided form | week10_fair_baselines |
| "residual search cannot recover the offset" | false with a fine grid and Np >= 2; evidence needed with Np = 1, 5-58x cheaper otherwise | week9_resid_vs_evidence |
| offset rms "5 ns" | physical offsets (BS-clock-locked RIS, offset = -tau_b) | week9_timing_physics |
| projection error ~1e-4 | only for shifts below ~100 ns | week9_timing_landscape |
| delay window fixed | window start chosen jointly by evidence; must start at or before the sync reference | week9_sync_window |
| two-timescale learning, TDL 300 ns pass | WITHDRAWN (cascade exceeds the CP) | d1_cp_check |
| TDL-D 100 ns "weakly violates the CP" | valid (-56 dB beyond CP) | d1_cp_check |
| coded gap 0.4 dB | 0.14 dB [0.09, 0.20] | week11_coded_stats |
| full-band SER ratio 1.09x / 1.54x | 1.10 [1.06, 1.15] / 1.37 [1.18, 1.86] | week11_fullband30 |
| beats CE-BEM by 1.5-4 in overhead | (pending held-out E4) earlier 1.1-4.8, median 2.7, selection on test channels | week11_overhead_heldout |
| net 3.9-5.7 bit/s/Hz on TDL | (pending held-out E4) | week11_overhead_heldout |
| CG 2.6x faster than dense | limited to N = 12 symbols per estimate; amortised comparison rerun pending (H2) | week11_eq_amortised |
| CG tolerance 1e-2 | 3e-3 (1e-2 gives 1.15x SER at M = 3276) | week11_eq_amortised_loaded |
| banded equalisers fail because of the 1/h^2 tail | fail also with smooth edges: off-diagonal energy outside the window caps SINR | week10_risetime |
| CE-BEM failure is intrinsic | depends on sharp switching edges (0.3-slot ramp: 8.2 vs 8.56) | week10_risetime |
| TM-RIS helps the user (implicit motivation) | time modulation costs the user rate; ours recovers it to within 0.1-0.2 of perfect; a static optimised RIS gives more for one user; motivation is security | week11_motivation |
| designed-user table "proposed" | labelled model-aware rank-1 estimator (not rerun with timing-only) | week3_designed_user |
| pilot-fit fallback in Algorithm 1 | removed (overfitting trap) | week11_nofallback |
| low rank (general) | strong for large Q (rank 16 of 1024), weak at Q = 16 (5 of 16) | week11_cdl_rank |
| achieved error "near optimal" | 2.5-6.5 dB above the Bayesian bound under the assumed prior (prior mismatch) | week11_theory |
| parametric element fit would do better | not better (+0.08 worse) and fails with element spread | week10_parametric |
