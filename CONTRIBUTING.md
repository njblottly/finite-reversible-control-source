# Contributing to Simultaneous Circuit Tests for Finite Reversible Models of Quantum Control

Thanks for considering a contribution! This repository contains the manuscript, Python verification scripts, exact catalogue data, witness transition tables and figure-generation code used in the paper.

This guide explains how to run the computations and contribute changes while keeping the mathematical claims and exported results reproducible. No quantum processor, laboratory equipment, GPU or Mathematica installation is required.

## Conduct

Keep discussions respectful, constructive and focused on the work. Explain disagreements with evidence, a reproducible example or a mathematical argument. Follow `CODE_OF_CONDUCT.md` if one is included in the repository.

## Repository layout

The verification and plotting scripts live at the repository root and resolve their data paths relative to their own location.

| Path | Purpose |
| --- | --- |
| `certify_catalogue.py` | Exact arithmetic, catalogue construction, matching routines and optional circuit search. |
| `verify_results.py` | Checks the four-frame example, the repeated-circuit horizon and the short-word witness. |
| `generate_figures.py` | Generates the manuscript plots and seeded sampling experiment. |
| `data/catalogue_geometry.json` | Admissible quotient edges and phase shifts. |
| `data/repeated_witness.json` | Repeated-circuit witness, lifted permutation and exhaustive search results. |
| `data/length12_feasible.json` | Witness for all nonempty Hadamard--phase words through length twelve. |
| `data/verification_report.json` | Verification results, execution environment and verifier source hashes. |
| `data/sampling_experiment.json` | Simulation settings, software versions and rejection frequencies. |
| `figures/` | Manuscript figures in PDF format and generated PNG previews. |
| Manuscript `.tex` and bibliography `.bib` files | Paper source and references; preserve the filenames and citation keys used by the current manuscript. |
| `README.md` | Project overview and reproduction instructions. |

## Quick start

After cloning the repository, run the following commands from its root.

### 1. Create a branch and Python environment

Use Python 3.10 or later. The recorded verification run used Python 3.12.14.

```bash
git checkout -b feat/short-description
python3 -m venv .venv
source .venv/bin/activate
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
```

### 2. Verify the exact results

The exact verification scripts require only Python's standard library.

```bash
python verify_results.py --rebuild-geometry
```

This rebuilds the catalogue geometry before checking the certificates. For a run using the supplied geometry, use:

```bash
python verify_results.py
```

A successful run ends with `All certificate checks passed.` The current benchmarks include 174 admissible quotient edges, a largest feasible repeated-block depth of 192, exclusion through depth 193, and verification of all 8,190 nonempty Hadamard--phase words of length at most twelve.

Verification rewrites `data/repeated_witness.json` and `data/verification_report.json`. Rebuilding geometry also rewrites `data/catalogue_geometry.json`. Review these changes before committing.

### 3. Reproduce the figures and sampling experiment

Install the plotting and sampling dependencies:

```bash
python -m pip install numpy matplotlib
python generate_figures.py
```

The script writes PDF figures and PNG previews to `figures/` and updates `data/sampling_experiment.json`. The experiment uses seed `20261003` and 20,000 simulated data sets per generating law and shot count. These are simulated Bernoulli outcomes, not laboratory measurements.

Record the Python, NumPy and Matplotlib versions when reporting reproduced or changed results. A fixed seed alone does not guarantee identical output across different software versions.

## Optional exploratory search

```bash
python certify_catalogue.py --search 12
```

This writes `data/circuit_search.json`. The search uses a two-million-node limit per feasibility call. A run that reaches its limit is unresolved; it is not an exclusion certificate. Keep exploratory findings distinct from the checked results reported in the paper.

## Mathematical and numerical contributions

State which model assumptions your change uses, including calibration, symmetry, stationarity, preparation and gate-accuracy constraints.

- Preserve exact comparisons in certificate decisions. Floating-point output can illustrate a result but must not replace its exact verification.
- Distinguish a witness that fits the tested circuits from a proof that every allowed controller fails.
- Explain why any new pruning rule cannot discard a feasible controller.
- Treat unresolved arithmetic signs and incomplete searches as unresolved, not as successful checks.
- If a change alters a benchmark, explain the cause and update the affected manuscript statements, data and figures together.

Search-node counts can change with traversal order. Timing and platform fields can also differ between machines. Explain meaningful differences rather than requiring incidental outputs to match byte for byte.

## Manuscript, references and figures

Keep existing citation keys stable unless the corresponding citations are updated throughout the manuscript. Supply a verifiable source for new references and distinguish established results from claims introduced here.

Make figure changes in the generation code and regenerate the affected outputs. Check labels, legends, units and captions for readability and overlap. Separate schematic illustrations from computed data in the caption.

For manuscript changes, confirm that the current `.tex`, `.bib` and figure files compile together in Overleaf or a local LaTeX installation. Do not commit temporary compilation files.

## What to commit

Commit source changes, relevant documentation and the curated data or figures needed to reproduce the affected result. Keep mathematical fixes and unrelated formatting changes in separate commits where practical.

Do not commit virtual environments, Python caches, editor files, temporary LaTeX files, large exploratory dumps or unrelated generated outputs. Preserve the provenance of witness files and distinguish newly generated reports from earlier archived runs. Keep any tracked checksum manifest consistent with the files it covers.

For large data, provide a stable archive link and document how to obtain and verify it. Do not silently replace exact witness data with rounded values.

## Issues and pull requests

For a bug report, include the command, relevant input, software versions, expected result and actual output. Provide a small reproducing example where possible.

In a pull request, describe:

1. The problem or proposed improvement.
2. The files and mathematical claims affected.
3. The verification commands you ran and their outcomes.
4. Any regenerated outputs, changed dependencies or remaining limitations.

Run the checks relevant to your change. Changes to exact arithmetic, geometry, matching, controller constraints or witness handling should include the full geometry rebuild and verification. Documentation-only changes do not require rerunning unrelated computations.

Follow the licence files included in the repository, retain attribution for reused material, and disclose substantial computational or AI assistance where it affects the provenance or interpretation of the contribution.