# Time-modulated RIS + OFDM: literature map and candidate paper (2026-10-01)

## Anchor paper
Iudice, Darsena, Gelli, Galdi, "Physically Consistent Modeling of Dispersive
Time-Modulated RIS for Wideband OFDM", arXiv 2609.18360 (16 Sep 2026).

* LPTV model. Each element switches periodically with period T = T_u (the OFDM
  useful symbol), with K slots, and each state has a single-pole dispersive response
  Gamma_k(f) = -1 + 2 xi_r / (j 2 pi (f - f0) + xi).
* Harmonic h couples transmitted subcarrier m to received subcarrier m+h, weighted by
  the element response at the absolute frequency of m. The coupling matrix is
  non-diagonal and, under dispersion, NOT Toeplitz.
* Generalised CP: T_cp >= delay spread + RIS memory T_gamma (13.8 ns for the OpenRIS cell).
* Validation: OpenRIS varactor cell, n78, 100 MHz, 30 kHz SCS. Table I gives fitted
  f0, xi_r, xi_i for the two states.
* Explicitly out of scope / future work: joint equalisation across coupled
  subcarriers, Tx precoding, guard subcarriers, control sequences that suppress
  harmonics, multi-resonant elements, mutual coupling, non-ideal transients,
  oblique incidence, OTA validation.

## Who uses TM-RIS/TMA with OFDM (all assume IDEAL, frequency-flat switching)
* Xu & Petropulu, TM-IRS for waveform security, arXiv 2310.13210 (ICASSP 2024)
* Tao & Petropulu, Secure TM-IRS via GFlowNets, arXiv 2506.14992. The legitimate user
  treats harmonic ICI as interference.
* Tao, Petropulu, Poor, TM-IRS for ISAC + security, arXiv 2509.05565. Future work:
  "hardware impairments".
* Tao, Xu, Petropulu, "How secure is TMA-enabled OFDM DM?", arXiv 2310.08551 (TWC 2025).
  An eavesdropper recovers symbols by blind ICA, resolving ambiguities with the
  TOEPLITZ structure of the mixing matrix (ideal model).
* Verde, Darsena, Galdi, rapidly time-varying RIS, arXiv 2304.11912. No OFDM.
* TMA antenna literature has studied non-ideal switch rise times and element bandwidth
  as a sideband filter (out-of-band), but not in-band OFDM coupling.

## Gate 0 numbers (Iudice Table I, computed here)
OpenRIS cell across 100 MHz: common phase drift 47.3 deg, differential 180 deg +/- small,
imbalance < 0.6 dB, DC-harmonic leakage -19..-27 dB (ideal: -inf). So dispersion is a
second-order effect on the harmonic magnitudes, but it breaks the Toeplitz structure.

## Second novelty sweep (2026-10-01): related areas that must be cited and differentiated
* BEM (basis-expansion model) estimation of doubly selective OFDM channels
  (CE-/P-BEM; e.g. "A BEM for estimation of time-varying channels in OFDM"; RIS-OFDM
  BEM + CNN, Shao et al., Syst. Eng. Electron. 2025). There the time variation is
  *unknown* (mobility/Doppler) and expanded on a chosen basis. Ours: the time
  variation is the *known* RIS modulation (sequence + element physics), so only the
  static two-hop channel is unknown, with rank-1 structure per element.
* PLC linear periodically time-varying channels (mains-synchronous, 10–20 ms period).
  The period is much longer than the OFDM symbol, so there is no harmonic coupling at
  subcarrier spacing. Different regime.
* TMA retrodirective DM with a legitimate-user pilot (Tao/Xu/Petropulu): the pilot is
  used by the *array* to set weights, not by the receiver to estimate the mixing matrix.
* Space-time-coding metasurface sensing/DoA work: harmonic processing for sensing,
  not communication channel estimation.
Result: still no paper on receiver-side estimation/equalisation of the
TM-RIS harmonic-coupled OFDM channel with known modulation.

## Gap (searched; no paper found)
No paper estimates or equalises the harmonic-coupled OFDM channel of a TM-RIS at a
legitimate receiver with pilots, under either the ideal or the dispersive model.
Every TM-RIS paper treats the coupling as interference. The only "undo the mixing" work
is the blind eavesdropper attack above (ideal, Toeplitz, TMA).

## Third novelty sweep, claim by claim (2026-10-05)
| claim | verdict |
|---|---|
| C1 receiver-side estimation + equalisation for TM-RIS OFDM | novel (nothing found; Iudice et al. list it as open) |
| C2 timing-only estimation (dispersion absorbed on input side) | novel |
| C3 low-rank two-hop delay structure | known idea (low-rank/tensor RIS cascaded-channel estimation); adopt and cite |
| C4 cyclic shift = diagonal modulations | identity known (STC: time delay -> -2 pi m tau / T harmonic phase, Zhang/Cui 2018); its use for estimation is new |
| C5 matrix-free CG equaliser | known technique: Hrycak, Das, Matz, Feichtinger, IEEE TSP 58(11) 2010 (product-convolution operators + LSQR/GMRES, BEM). Our only addition: the TM-RIS coupling is EXACTLY product-convolution |
Paper claims only C1, C2 and the use in C4; C3/C5 are cited as adopted methods.

## Fourth novelty sweep (2026-10-21), five weeks after the anchor paper

Searched arXiv / web for follow-ups and for the new week-8 ideas.

**Scoop check: none found.**
- Anchor authors (Iudice arXiv page): no TM-RIS paper after arXiv:2609.18360.
- arXiv:2610.12773 (Wang, Xiang, Petropulu, 8 Oct 2026), generative (flow-matching) control of
  TM-IRS for directional modulation with OFDM. Transmitter-side pattern design; no receiver,
  estimator or equaliser. Cite as concurrent TM-IRS + OFDM work.
- No paper found on channel estimation or ICI equalisation for OFDM through a time-modulated /
  space-time-coded RIS. Closest: time-coded RIS + OFDM tensor DOA/delay estimation (Sensors 2022,
  doi 10.3390/s22187083), STC-RIS localisation (arXiv:2410.20505), STMM modelling (arXiv:2302.08310).
  None estimates the harmonic coupling for detection.

**Week-8 ideas: what is new and what must be credited.**
- Evidence-based blind timing (L2b). The TOOL is known: Bayesian model selection for OFDM channel
  estimation under residual timing offsets (IEEE Access 2019, eprints.gla.ac.uk/188514); model
  averaging over timing offset (IEEE TWC letter 2010, hub.hku.hk/handle/10722/124742); AIC/BIC for
  joint channel, length and time offset (LATINCOM 2016). Related RIS timing work: blind sync-offset
  estimation for RIS symbiotic radio (Chen, Long, Liang, Schober, ICC 2024, doi
  10.1109/ICC51166.2024.10622298): the RIS backscatters its own symbols with an unknown offset that
  causes ICI, estimated via virtual-subcarrier energy detection; ML timing offsets in distributed
  multi-RIS (arXiv:2205.02426, 2107.12717). What appears NEW: for a TM-RIS, the switching offset is
  exactly an opposite delay shift of the two hops (T_d o (FU)(FV)^T = T o (D FU)(D* FV)^T), which makes
  the pilot residual nearly flat in the offset (why residual search fails at high SNR), and the
  evidence under the DPSS prior resolves it. Claim the insight and the application, not the tool.
- Two-timescale estimation (L1/L1b/L1c). WELL KNOWN for conventional RIS: quasi-static common BS-RIS
  channel estimated rarely, RIS-user often (Hu and Dai, arXiv:1912.07990; cell-free RIS
  arXiv:2207.14182; double-RIS arXiv:2211.16183; near field arXiv:2501.02985). Not a novelty claim.
  Our contribution at most: the same principle for the harmonic-coupled TM-RIS model, with the
  stacked-SVD initialisation of the shared input side. Must cite these.
- Matrix-free CG for non-diagonal OFDM: known (Hrycak et al. 2010, already cited).

**Verdict.** The core claim (timing-only estimator plus receiver for TM-RIS OFDM, no element model)
is still unscooped. Evidence-based timing is a solid secondary contribution if credited as above.
Two-timescale learning is an application of a known principle.
