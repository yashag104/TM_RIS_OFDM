# Summary D: channel models, CP validity, spatial model, rank (2026-10-22)

## D1 CP validity — `d1_cp_check.json` (computed)
Criterion: cascaded BS-RIS-user energy beyond the CP (2.34 us, element memory 13.8 ns) below -40 dB.
TDL-A/C/D at 29 and 100 ns: valid (TDL-D 100 ns: -56 dB, the earlier "weak violation" wording was wrong).
All 300 ns cases: invalid (-25 to -17 dB). Consequence: every 300 ns result (two-timescale learning) WITHDRAWN.

## D2 CDL vs TDL — `week11_cdl_rank.json` (30 channels, Q = 16, 100 ns)
Gap of ours to perfect CSI (TDL / CDL), Np 2 / 4 / 8, 20 dB:
- A: 0.37/0.44, 0.19/0.21, 0.085/0.087
- C: 0.55/0.50, 0.22/0.21, 0.088/0.084
- D: 0.73/0.71, 0.42/0.47, 0.22/0.24
CDL with cluster angles behaves like TDL + random directions (difference <= 0.07 at Np 2, <= 0.03 at Np 8). The
paper keeps TDL as the main 3GPP model and reports CDL as the spatially consistent check.

## D3 rank of the element-domain matrix
Effective rank (99% energy), Q = 16 / 64 / 256 / 1024:
narrow far 3.5/4.9/5.3/5.5; wide (98 MHz) far 5.0/8.9/11.9/14.2; wide + squint 5.0/8.8/12.9/16.0;
near field 3 m 3.2/4.9/6.1/7.5; 10 m 3.5/4.6/6.6/7.4.
Low rank is strong for large Q (16 of 1024) and weak at Q = 16 (5 of 16). Corollary on low rank restated accordingly.
