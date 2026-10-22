# Summary I: does the problem matter? (I1-I3), 2026-10-22
Fixed transmit power and noise (sigma^2 set so the optimised static RIS has the reference SNR). Designed TM pattern
(per-element polarity co-phasing h = 0 + random cyclic shifts; 8-11 distinct patterns), 30 channels per point.

20 dB, rate (bit/s/Hz): static optimised 6.66; static random 3.1-3.4; no RIS (direct -10 dB) 3.46;
TM one-tap 0.68 / 1.98 / 3.15 (D = 0.625 / 0.75 / 0.875); TM ours (Np 4) 2.98 / 4.06 / 5.35; TM perfect 3.08 / 4.18 / 5.50.
Ours - one-tap: +2.1..+2.3 (20 dB), +4.9..+5.2 (30 dB), +0.3..+0.7 (10 dB). Ours - perfect: 0.1-0.2.
Ours - static optimised: -1.3 .. -3.7 (20 dB): for a single user time modulation does not pay off.
Eavesdropper (other position, same pilots), 20 dB: one-tap 0.08-0.71; generic CE-BEM |h|<=40 without the pattern
1.8-2.2; our receiver WITH the pattern 2.5-3.0. Secrecy rate vs these: 2.9-4.6 / 1.0-3.6 / 0.5-2.9.
K variation (4, 8, 16 at D = 0.75) changes little.

Consequences for the paper:
* Motivation reframed: time modulation is used for security / harmonic steering / sensing; our receiver removes most
  of its cost to the legitimate user (to within 0.1-0.2 of perfect CSI). It does not make TM-RIS better than an
  optimised static RIS for a single user (1.3-3.7 bit/s/Hz below at 20 dB).
* Security caveat: a generic CE-BEM eavesdropper already recovers much of the rate; knowing the pattern helps it
  further. Pattern secrecy (A1 tension) must be stated; security analysis of TM-IRS designs is outside our scope.
