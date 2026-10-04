# Simultaneous Circuit Tests for Finite Reversible Models of Quantum Control

## Reproduce exact results

Python 3.10 or later; standard library only. No GPU, Mathematica or network access is needed. Decisions use rational arithmetic in Q(sqrt(2), sqrt(3), sqrt(5), sqrt(7), i). Unresolved signs raise an error. Decimal values are descriptive, not exclusion evidence.

The reconstructed traversal differs from the original archive, so its counts differ. Its mathematical endpoint is the same. The negative endpoint phase convention at t=0 is retained under right multiplication by T; full lifted transitions are checked against physical matrix edges.

## Reproduce plots and the sampling experiment

Install NumPy and Matplotlib in your chosen environment, then run:

```bash
python3 generate_figures.py
```

The simulation uses seed 20261003 and 20,000 replicates per law and shot count. It uses synthetic Bernoulli outcomes from known probabilities; no laboratory observations are represented or fitted.

## Scope and provenance

The finite-frame model, shared-history formulation, matching method and the 192/193 benchmark developed in: Reversible Finite-Frame Dynamics for Rational Quantum Grids (research archive acknowledged in the manuscript). Here, we study the explicit four-frame two-circuit separation, a formulation for arbitrary finite collections, reconstructed verification, exhaustive short-word feasibility witness and finite-shot analysis.