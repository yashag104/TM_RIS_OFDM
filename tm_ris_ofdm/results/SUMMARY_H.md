# Summary H: equaliser cost (H2 amortised wall-clock, H3 tolerance)

Status: H2 and H3 done (2026-10-24). Source: `results/week11_eq_amortised.json`, from a clean rerun on an idle machine
(4 threads, MATLAB R2026a, 20 dB, 3 repetitions). The first run overlapped another job, so it is kept as
`week11_eq_amortised_loaded.json` and not used (deviation logged in PREREG_LOG).

## Wall-clock time (s) for N symbols per estimate, CG at tolerance 3e-3
"Dense" includes forming and factorising the matrix once per estimate.

| M | dense, N = 1 / 12 / 40 / 120 | CG, N = 1 / 12 / 40 / 120 | CG iterations | SER dense / CG |
|---|---|---|---|---|
| 1024 | 0.14 / 0.14 / 0.14 / 0.16 | 0.03 / 0.41 / 1.24 / 3.65 | 8 | 0.0001 / 0.0001 |
| 2048 | 0.87 / 0.88 / 0.91 / 1.01 | 0.16 / 1.68 / 4.66 / 13.5 | 16 | 0.0312 / 0.0313 |
| 3276 | 2.51 / 2.54 / 2.58 / 2.78 | 0.23 / 2.51 / 8.20 / 23.8 | 20 | 0.189 / 0.193 |
| 6552 | not run (687 MB matrix; pre-registered as CG only) | 0.47 / 4.90 / 17.1 / 51.2 | 19 | n/a / 0.132 |

## Findings
- **Break-even:** about 12 symbols per estimate at M = 3276, 6 at M = 2048, and 4 at M = 1024.
  - This is far below the operation-count estimate of about 140 symbols, because dense LAPACK runs at level-3 BLAS speed while CG is memory-bound.
- **When CG is faster:** for the 6–15-symbol refresh blocks selected at 30 and 120 km/h (SUMMARY_E).
- **When dense is faster:** for slow channels (40–120 symbols per estimate), by 3–9x.
- **CG's lasting advantage:** memory linear in M. It handles M = 6552 without the dense matrix.
- **H3, the tolerance:** at 3e-3, CG's SER is within 2% of dense at every M. At 1e-2 it is 1.15x dense at M = 3276. The paper uses 3e-3.
- **Correction:** the old claim "CG 2.6x faster (1.12 vs 2.95 s)" used tolerance 1e-2 and did not amortise the factorisation. It is withdrawn and replaced by the break-even statement in the IEEE paper.
