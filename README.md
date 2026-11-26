# TENSOR-SPEC-CERT-1.0: Constructive Rational and Exterior-Algebra Certificates for Euclidean Geometry in the Asymptotic Spectrum of Tensors

**Author:** Lluis Eriksson (`lluis@airr.science`)  
**Schema:** `TENSOR-SPEC-CERT-1.0`  
**Target Preprint:** OpenAI, *"Euclidean geometry in the asymptotic spectrum of tensors"* (`OAI:Euclidean-geometry-in-the-asymptotic-spectrum-of-tensors-November-24-2026`, unrefereed preprint, November 24, 2026).

## Overview
This repository provides an independently reproducible, zero-floating-point Python verification engine (`src/tensor_spectrum_verifier.py`), a constructive Lean 4.35.0 kernel (`lean/TensorSpectrumCert.lean`, 19 constructive theorems, 0 `sorry`, 0 `axiom`, 0 `Classical.choice`, 0 warnings), 16 canonical JSON certificates (`certificates/C01..C16`), 88 unit and adversarial forgery tests (`tests/test_verifier.py`), 3 publication figures, and a 4-page LaTeX/PDF companion paper (`paper/main.pdf`).

## Quick Reproduction
```bash
# 1. Build Lean 4.35.0 kernel
cd lean && lake build && cd ..

# 2. Verify all 16 canonical certificates
python -m src.tensor_spectrum_verifier certificates/

# 3. Run full pytest & adversarial forgery suite
pytest -q

# 4. Verify SHA-256 manifest
python scripts/build_manifest.py --check
```

## Licenses & Disclosure
- Manuscript & figures: `CC-BY-4.0`
- Code, Lean proofs & certificates: `MIT`
- Produced with Google Antigravity assistance; independently compiled, verified, and reviewed by Lluis Eriksson.
