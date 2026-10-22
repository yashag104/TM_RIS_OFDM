# Summary J: novelty check (content), 2026-10-22

Searches (web/arXiv, 2026-10-21 and 2026-10-22): time-modulated RIS channel estimation; space-time-coding
metasurface OFDM receiver; intra-symbol RIS switching ICI equalisation; OFDM backscatter with known tag
modulation timing; periodically time-varying channel estimation with known gating; RIS timing offsets;
Bayesian evidence for timing. No work found that estimates and equalises the harmonic-coupled OFDM channel
of a time-modulated RIS at the legitimate receiver with known switching timing.

| # | Closest work | What it does | What is new in ours |
|---|---|---|---|
| 1 | Iudice, Darsena, Gelli, Galdi, arXiv:2609.18360 (2026) | Physically consistent dispersive TM-RIS OFDM model; coupling non-diagonal and non-Toeplitz; receiver left open | Receiver: estimation and equalisation; state regrouping of their eq. (5) (credited) used to absorb dispersion |
| 2 | Xu & Petropulu, arXiv:2310.13210 (ICASSP 2024) | TM-IRS for waveform security; legitimate user conventional | We estimate/equalise the coupling instead of designing it away |
| 3 | Tao & Petropulu, arXiv:2506.14992; Tao, Petropulu, Poor, arXiv:2509.05565 | TM-IRS pattern design (GFlowNets), ISAC + security; legitimate user treats ICI as interference | Same as 2 |
| 4 | Wang, Xiang, Petropulu, arXiv:2610.12773 (Oct 2026) | Generative (flow-matching) TM-IRS control for directional modulation with OFDM | Transmitter-side; no receiver |
| 5 | Tao, Xu, Petropulu, arXiv:2310.08551 (TWC 2025) | Eavesdropper recovers TMA/OFDM DM symbols by blind ICA, using the ideal-model Toeplitz structure | Pilot-based, legitimate, dispersive (non-Toeplitz) model; no element model |
| 6 | Tang, Cannizzaro, Leus, Banelli, TSP 2007; Zemen & Mecklenbraeuker, TSP 2005 | Basis-expansion / DPSS estimation of doubly selective OFDM; time variation unknown | Time variation known (switching timing); only propagation + dispersion estimated. Baseline in all tables |
| 7 | Hrycak, Das, Matz, Feichtinger, TSP 2010 | Matrix-free equalisation via product-convolution operators (BEM approximation) | Adopted; for a TM-RIS the product-convolution form is exact |
| 8 | Hu, Dai et al., arXiv:1912.07990 | Two-timescale RIS channel estimation (static BS-RIS hop) | Adopted principle; applied to the harmonic-coupled model (no novelty claim) |
| 9 | Chen, Long, Liang, Schober, ICC 2024 | RIS symbiotic radio: RIS backscatters its own data; unknown sync offset causes ICI; blind offset via virtual subcarriers | Offset = opposite hop delay (exact identity); evidence over the joint (offset, delay-window) choice; legitimate OFDM link with known pattern |
| 10 | TScatter (Liu et al., NSDI 2021); ICI reduction for backscattered OFDM (2024); interference-free OFDM symbiotic radio (arXiv:2408.09192) | Tag switches within OFDM symbols; receiver decodes the tag's data; ICI handled by unused subcarriers or design | Inverse problem: we know the switching and estimate the channel; no unused subcarriers needed |
| 11 | ISAC via STC metasurfaces (Nat. Commun. 2025); STC differential modulation for RIS (arXiv:2508.10244) | Harmonic-amplitude DOA; non-coherent RIS modulation avoiding channel estimation | Coherent estimation and equalisation of the coupled OFDM channel |
| 12 | Bayesian model selection for OFDM timing offsets (IEEE Access 2019); MacKay 1992 | Evidence/model selection for timing and model order | Tool adopted; new: why it is needed here (offset == opposite hop delays, residual near-flat for small offsets) |
| 13 | Cyclic short-time-varying PLC channel estimation | LPTV channel, period >> symbol, no harmonic coupling | Period = symbol, harmonic coupling at the subcarrier spacing |

Wording change: "the first receiver" is softened to "to our knowledge, the first pilot-based receiver for OFDM
through a dispersive time-modulated RIS"; the closest works (5, 9, 10) are named in the introduction.
Novelty that survived the review: (i) timing-only estimation absorbing the element dispersion, validated
against an independent time-domain truth and element transients (A1); (ii) the receiver as a whole; (iii)
the opposite-delay identity and its consequences for timing (B). Not claimed: low rank, DPSS, CG equalisation,
two-timescale learning, evidence as a tool.
