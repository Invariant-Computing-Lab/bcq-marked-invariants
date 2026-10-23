# BCQ Invariant Geometry - Reproducibility Package

This package accompanies the manuscript **"Quotient-Invariant Geometry for the Brane-Constituent Quark Hypothesis: Algebraic Reduction and Synthetic Falsification"**.

## Quick run

```bash
cd code
python run_all.py
pytest -q
```

Outputs are written to `../results/` and `../figures/`.

## Package structure
- `docs/` - main manuscript and supplementary appendix (DOCX and PDF)
- `code/bcq_invariants.py` - invariant primitives and transformation generators
- `code/run_all.py` - complete deterministic synthetic experiment suite
- `code/test_invariants.py` - unit tests for exact invariance properties
- `results/summary.csv` - compact result table
- `results/run_metadata.json` - seed, environment, and experiment metadata
- `results/E2_chi3_collision_vectors.csv` - exhaustive chi=3 collision set
- `results/E6_discriminant_curve.csv` - discriminant/identifiability stress table
- `results/E8_ablation_detection.csv` - ablation detection rates with Wilson 95% CIs
- `figures/` - publication-resolution PNG and vector PDF figures
- `PROTOCOL.md` - locked scope, thresholds, and interpretation rules

## Evidence boundary
The package tests formal invariance and reconstructibility in explicit synthetic testbeds. It does not establish the physical correctness of BCQ, a Standard Model embedding, confinement, a valid compactification, or quark compositeness.
