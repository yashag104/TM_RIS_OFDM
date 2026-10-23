# Summary K: hardware realism (2026-10-22)

## K1: not done (over-the-air / measured element). Stated as the main limitation.

## K2 realistic rise times — `week11_rise_small.json`
Rise times 2-42 ns: nominal kernels are enough, no calibration needed. 83 ns: calibrated kernels needed
(0.14 bit/s/Hz at 30 dB without).

## K3 out-of-band emission — `k3_oob.json` (analytic)
ACLR of the switched reflection (ideal +-1, D = 0.75, K = 8, 100 MHz): 36.8 dB with ideal edges, below the 45 dB
BS requirement; 42.3 dB at rise 0.01 slot and 52.7 dB at 0.05 slot. Trade-off in the paper: smoother edges lower
emission but shrink the advantage over truncated estimators (C5).
