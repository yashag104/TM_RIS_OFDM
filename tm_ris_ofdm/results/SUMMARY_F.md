# Summary F: full band and coded statistics (2026-10-22)

## F1 full band, 30 channels — `week11_fullband30.json`
M = 3276. SER ratio ours / true channel (clustered bootstrap 95% CI, 10^4 resamples):
pilots only 1.10 [1.06, 1.15] at 20 dB / 16-QAM, 1.37 [1.18, 1.86] at 28 dB / 64-QAM; with one DD pass 1.05 / 1.12.
CE-BEM 18x / 32x. Replaces the earlier 1.09x / 1.54x from 10 channels.

## F2 coded BLER — `week11_coded_stats.json`
200 channels x 10 blocks = 2000 blocks per SNR. 10% BLER crossing gap ours - perfect: 0.14 dB [0.09, 0.20]
(earlier 0.4 dB from 120 blocks). One-tap never decodes; CE-BEM |h| <= 12 floors at 0.34; |h| <= 40 fails with 2 pilots.

## F2 (timing regimes, workstream B numbering in the plan) — see SUMMARY_B.md.
