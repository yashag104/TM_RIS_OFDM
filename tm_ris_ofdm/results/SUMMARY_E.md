# Summary E: pilot overhead and mobility (held-out protocol)

**Status:** E4 done (2026-10-24).

## Sources
- **Data:** `week11_overhead_heldout.json` (shared path) and `week11_overhead_heldout_{TDLA100,TDLC100,TDLD100}.json`.
- **Analysis:** `figures/analyse_heldout.py`, which writes `week11_overhead_heldout_summary.json`.
- **Tables generated from it:**
  - `paper/journal/tab_overhead_heldout.tex`;
  - IEEE Table IV (`paper/ieee/main.tex`, rows between the `%ROWS` markers).

## Protocol (pre-registered)
- **Selection:** pick the pilot count Np and block length N (and the CE-BEM harmonic range) that maximise the mean net rate over 20 selection channels.
- **Reporting:** report on 20 disjoint held-out channels, as mean ± 95% t-interval and paired differences.
- **Receivers:** common-phase tracking; robust MMSE for every estimated receiver; CE-BEM allowed up to all harmonics.
- **Run history:** the TDL runs were killed by a shutdown on 2026-10-23 and rerun from the same seeds, so no numbers changed because of the rerun.

## Results (bit/s/Hz, held-out)

| Channel | km/h | 20 dB: ours / BEM / diff | 30 dB: ours / BEM / diff | Perfect CSI (20/30 dB) |
|---|---|---|---|---|
| Shared | 3 | 6.24 / 5.24 / 1.01±0.14 | 8.99 / 7.34 / 1.65±0.34 | 6.71 / 10.02 |
| Shared | 30 | 5.97 / 4.40 / 1.57±0.30 | 8.56 / 5.42 / 3.14±0.40 | 6.61 / 9.90 |
| Shared | 120 | 5.26 / 3.15 / 2.11±0.30 | 7.30 / 3.38 / 3.92±0.53 | 6.46 / 9.74 |
| TDL-A | 3 | 5.38 / 3.71 / 1.67±0.17 | 7.37 / 4.88 / 2.49±0.21 | 6.11 / 9.30 |
| TDL-A | 30 | 4.56 / 1.60 / 2.97±0.13 | 5.51 / 1.68 / 3.84±0.27 | 6.06 / 9.25 |
| TDL-A | 120 | 3.81 / 1.08 / 2.73±0.28 | 4.35 / 1.07 / 3.27±0.58 | 6.25 / 9.45 |
| TDL-C | 3 | 5.34 / 3.61 / 1.74±0.16 | 7.23 / 4.81 / 2.41±0.27 | 6.11 / 9.34 |
| TDL-C | 30 | 4.53 / 1.65 / 2.88±0.29 | 5.53 / 1.81 / 3.72±0.42 | 6.11 / 9.28 |
| TDL-C | 120 | 3.88 / 0.75 / 3.13±0.31 | 4.61 / 0.74 / 3.88±0.54 | 5.86 / 9.03 |
| TDL-D | 3 | 5.17 / 3.61 / 1.56±0.17 | 6.80 / 4.75 / 2.05±0.41 | 6.11 / 9.30 |
| TDL-D | 30 | 4.65 / 2.17 / 2.48±0.25 | 5.85 / 2.39 / 3.46±0.42 | 6.13 / 9.33 |
| TDL-D | 120 | 3.92 / 1.15 / 2.77±0.31 | 4.68 / 1.15 / 3.53±0.69 | 6.12 / 9.30 |

## Findings
- **Ours beats the best CE-BEM in all 24 cases**, by 1.01–3.92 bit/s/Hz. Every interval is above zero.
- **CE-BEM's selected settings:** 16–32 pilot symbols, with |h| ≤ 128 or all harmonics.
- **Selection bias is small:** in-sample selection on the held-out channels differs by at most 0.13. The earlier test-channel table was not badly biased.
- **Fraction of the perfect-CSI rate that ours keeps:**
  - shared path: 75–93%;
  - TDL at 3 km/h: 73–88%;
  - TDL at 30 km/h: 60–76%;
  - TDL at 120 km/h: 46–66%.
- **Selected designs:**
  - slow channels: 2 pilot symbols per 60–120 (shared path) or 4 per 40–80 (TDL);
  - at 30 km/h: 1–2 per 8–30;
  - at 120 km/h: 1 per 6 (TDL) or 1 per 10–15 (shared path).
- **Negative result, reported honestly:** high mobility on TDL is a real limit of this pilot-based design. The coherence time is about 30 symbols at 120 km/h. The earlier "3.9–5.7 bit/s/Hz on TDL" claim is replaced by 3.8–7.4 (held-out). The IEEE paper's "75–93%" claim is now restricted to the shared path, and the 120 km/h figure is stated in its conclusion.
