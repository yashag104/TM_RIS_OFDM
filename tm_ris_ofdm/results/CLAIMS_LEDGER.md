# Claims ledger (started 2026-10-22, git 75f31e0)

Status: OK = survives current evidence; WEAK = holds but evidence thin / unfair baseline; FIXED = number
corrected; WITHDRAWN = no longer claimed; OPEN = test pending (workstream in brackets).

| # | Claim | Old number | Evidence so far | Status | Corrected number / note |
|---|---|---|---|---|---|
| C-A | Timing only, no element model | loss 0.08/0.11 dB @20/30 | week3_modelfree(_confirm), 30 ch, Np=2 + 2 DD | WEAK | uses 2 DD iterations; generator = model (inverse crime) [A1-A4]; pilots-only rerun pending |
| C-B | 2 pilots within ~0.1 bit/s/Hz of perfect | 6.52 vs 6.61 / 9.68 vs 9.92 | week4_bem_compare, 20 ch | OPEN | gap 0.09 / 0.24; generator delays match the prior [A2]; fair-equaliser baselines [C1] |
| C-C | Offset absorbed blindly by evidence | rms 5 ns; loss 0.40 vs 0.41 dB @30 | week8_timing_evidence, 12 ch, free uniform offset | FIXED wording / OPEN | "5 ns" is reference-dependent (Prop. 2); claim = no rate loss. Physical offsets, prior sensitivity, overfitting test [B1-B4] |
| C-D1 | Two-timescale learning | TDL-A/C 300 ns pass | week8_twoscale_300ns_c, 6 ch | WITHDRAWN | 300 ns cascades exceed the CP; ground truth omits ISI |
| C-D2 | Warm start | 12x ops / 8.6x time | week7_complexity_fix, week8_fullband_warm | OK | n = 8 / 4 channels [H1] |
| C-D3 | CG faster than dense | 2.6x at M=3276 | week8_eq_walltime3 | WEAK | 12 symbols incl. factorisation; amortisation untested [H2]; tolerance issue at M=2048 [H3] |
| C-E1 | Pilot saving vs CE-BEM | "16x" | week7_strong_baselines | FIXED | ours (2) exceeds full BEM (32); BEM never within 0.3 bit/s/Hz; Table I ">32". Fair equaliser pending [C1, C4] |
| C-E2 | Coded within ~0.5 dB | 0.4 dB crossing | week6_coded, 120 blocks/pt | WEAK | 3 vs 1 block errors; re-measure with >= 2000 blocks [F2] |
| C-E3 | Net 3.9-5.7 bit/s/Hz on TDL | TDL 100 ns, 4 ch | week8_pilot_overhead | WEAK | n = 4; (N, Np) chosen on the test channels; TDL-D 100 ns weakly violates CP [E4, D1] |
| C-E4 | Beats best CE-BEM in all 24 overhead cases | margin 1.1-4.8 | same | WEAK | BEM capped at |h|<=40 on TDL; unfair equaliser [C1, C4] |
| T-X | Designed-user table "proposed" | 9.37 @30 dB | week3_designed_user | FIXED label | model-aware rank-1 estimator, not timing-only; rerun pending |
| T-VI | Element-knowledge table | | week3_mismatch | FIXED label | 2 DD iterations; pilots-only pending |
