# tm_ris_ofdm

Receiver for OFDM through a time-modulated reconfigurable intelligent surface (TM-RIS):
a channel estimator that needs only the switching timing, and a matrix-free equaliser.
Builds on Iudice, Darsena, Gelli, Galdi, arXiv 2609.18360 (2026).

- `tm_ris_ofdm/`: the project. Start with `tm_ris_ofdm/README.md` (results, corrections,
  limitations) and the learning guide `tm_ris_ofdm/report/tm_ris_learning_guide_local.html`.
  MATLAB code in `matlab/`, results as JSON in `results/`, figures in `figures/`.
- `old_ris_project/`: the earlier RIS project (federated learning / link-level / NoC study),
  kept unchanged for reference. Nothing in `tm_ris_ofdm/` depends on it.

Figures: `.venv/bin/python tm_ris_ofdm/figures/make_figures.py` (numpy, matplotlib).
