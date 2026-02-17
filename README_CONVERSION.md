# Julia → Python Conversion Guide

This repository is undergoing conversion from Julia to Python. This document provides quick orientation.

## What This Code Does

**Weak gravitational lensing analysis pipeline** for cosmological simulations (HAGN/Horizon-AGN):
1. Reads deflection fields from N-body simulations
2. Computes lensing observables (convergence κ, shear γ, flexion F/G)
3. Measures tangential correlations around galaxies
4. Fits NFW halo profiles to infer masses
5. Generates publication plots for thesis

## Current State

**Original (Julia)**: 22 files, ~2500 duplicate lines, no tests, scattered organization  
**Target (Python)**: 6 modules, comprehensive tests, production pipeline, HPC-ready

## Documentation

- **Full Plan**: `MODERNIZATION_PLAN.md` (71 pages, all issues & solutions)
- **Status Tracker**: `IMPLEMENTATION_STATUS.md` (checklist & quick reference)
- **This File**: Quick orientation for developers

## Key Decisions (Already Approved)

✅ **Architecture**: Modular Python (6 modules: io, derivatives, observables, correlations, nfw, cosmology)  
✅ **Dependencies**: Pure Python + standard libraries (NumPy, SciPy, Astropy, TreeCorr, HEALPix)  
✅ **Testing**: Hierarchical (unit < debug < validation < production scales)  
✅ **Workflow**: Snakemake pipeline for reproducibility  
✅ **Spatial**: HEALPix decomposition for testing small patches  

## Quick Start (for Dev Agent)

### 1. Review the Plan
```bash
less MODERNIZATION_PLAN.md
# Or open in IDE
```

Key sections:
- **Approved Issues**: All 12 issues with solutions (search for "Issue #")
- **Modern Library Research**: What libraries to use & why (search for "TreeCorr")
- **Ordered Implementation Plan**: 6 phases with tasks (search for "Phase 1:")
- **Test Strategy**: What tests must pass (search for "Unit Tests")

### 2. Set Up Environment
```bash
# Create Python environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install core dependencies
pip install numpy scipy astropy matplotlib pytest

# Install weak lensing tools (CRITICAL)
pip install treecorr healpy

# Install workflow orchestration
pip install snakemake
```

### 3. Create Package Structure
```bash
# Main package
mkdir -p cosmo_lensing
touch cosmo_lensing/__init__.py

# Module files (to be implemented)
touch cosmo_lensing/io.py
touch cosmo_lensing/cosmology.py
touch cosmo_lensing/derivatives.py
touch cosmo_lensing/observables.py
touch cosmo_lensing/correlations.py
touch cosmo_lensing/nfw.py

# Testing
mkdir -p tests
touch tests/conftest.py
touch tests/test_io.py
touch tests/test_observables.py
touch tests/test_against_treecorr.py

# Scripts (production use)
mkdir -p scripts
touch scripts/compute_observables.py
touch scripts/compute_correlations.py
touch scripts/fit_nfw_profiles.py
touch scripts/generate_plots.py

# Workflow
mkdir -p workflow
touch workflow/Snakefile
touch workflow/config.yaml
```

### 4. Start Implementing (Phase 1)

**Recommended order** (see MODERNIZATION_PLAN.md Phase 1 for details):

1. **io.py**: Start with `read_deflection_field()`
   - Reference: `raytrace_tk.jl` lines 77-110
   - Add error handling (Issue #6)
   - Write test: `tests/test_io.py`

2. **cosmology.py**: `CosmologyCalculator` class
   - Wraps `astropy.cosmology.FlatLambdaCDM`
   - Replaces `libsoftlens.so` (Issue #2)
   - Validate against Julia code

3. **derivatives.py**: `compute_jacobian()`
   - Reference: `raytrace_tk.jl` lines 150-326
   - Compute ∂α/∂x (first derivatives)
   - Write tests with synthetic data

4. **observables.py**: κ, γ, F, G functions
   - Reference: `raytrace_tk.jl` lines 659-760
   - Add docstrings with citations (Issue #7)
   - Test against analytic solutions

## File Mapping (Julia → Python)

| Julia File | Purpose | Python Equivalent |
|------------|---------|-------------------|
| `raytrace_tk.jl` (54KB) | Monolithic toolkit | → Split into 6 modules |
| `core.jl` | libsoftlens.so wrapper | → Delete (use astropy) |
| `cosmo.jl` | Cosmology setup | → `cosmology.py` |
| `nfw_tk.jl` | NFW utilities | → `nfw.py` |
| `correlation_function.jl` | Correlations | → `correlations.py` (use TreeCorr) |
| `correlation_plots.jl` | Main analysis | → `scripts/compute_correlations.py` |
| `rebuild_*_maps.jl` (10 files!) | Map generation | → `scripts/compute_observables.py` (one script) |
| `correlation_plots_sis_model*.jl` | Dead code | → **DELETE** |

## Critical Validation Points

**Must achieve < 1% agreement** (Issue #10):
- [ ] Tangential shear vs. TreeCorr
- [ ] Convergence profiles vs. TreeCorr
- [ ] Flexion vs. analytic NFW (no reference implementation)

**Must achieve < 10% agreement**:
- [ ] Fitted halo masses vs. Julia code (systematic OK if documented)

**Must achieve 100% reproducibility**:
- [ ] Same inputs → same outputs (within 0.1% numerical precision)

## Common Pitfalls to Avoid

❌ **Don't**: Reimplement TreeCorr's correlation functions  
✅ **Do**: Wrap TreeCorr with thin interface (Issue #10)

❌ **Don't**: Keep libsoftlens.so dependency  
✅ **Do**: Use astropy.cosmology exclusively (Issue #2)

❌ **Don't**: Process full 20k×20k fields for debugging  
✅ **Do**: Use test hierarchy (100×100 → 1k×1k → 5k×5k → 20k×20k) (Issue #9)

❌ **Don't**: Couple visualization with computation  
✅ **Do**: Separate scripts (compute vs. plot) (Issue #4)

❌ **Don't**: Copy-paste code between scripts  
✅ **Do**: Single parameterized function (Issue #5)

❌ **Don't**: Silent failures (return None)  
✅ **Do**: Explicit exceptions with logging (Issue #6)

## Testing Strategy

Run tests at increasing scales:

```bash
# Unit tests (< 10 sec, synthetic data)
pytest tests/ --config unit

# Debug tests (< 10 sec, small real data patch)
pytest tests/ --config debug

# Integration tests (< 2 min, medium patch)
pytest tests/ --config validation --integration

# Full validation (< 30 min, large subset)
snakemake --config scale=validation --dry-run

# Production (10 min/slice, full dataset)
snakemake --config scale=production --profile slurm
```

## Key References

**Papers to cite** (see MODERNIZATION_PLAN.md for full list):
- **Flexion**: Bacon et al. (2006) MNRAS 365, 414
- **NFW Lensing**: Wright & Brainerd (2000) ApJ 534, 34
- **TreeCorr**: Jarvis et al. (2004) MNRAS 352, 338
- **HEALPix**: Górski et al. (2005) ApJ 622, 759

**Software Documentation**:
- TreeCorr: https://rmjarvis.github.io/TreeCorr/
- Astropy: https://docs.astropy.org/
- HEALPix: https://healpy.readthedocs.io/
- Snakemake: https://snakemake.readthedocs.io/

## Questions?

1. **Check the plan first**: `MODERNIZATION_PLAN.md` (search is your friend)
2. **Check implementation status**: `IMPLEMENTATION_STATUS.md`
3. **Check this guide**: `README_CONVERSION.md` (you are here)
4. **Then ask**: Create issue or discussion

## Timeline

**Phase 1 (Foundation)**: 2 weeks - Core modules + basic tests  
**Phase 2 (Validation)**: 1 week - Synthetic data + analytic tests  
**Phase 3 (TreeCorr)**: 1 week - Correlation functions + validation  
**Phase 4 (Pipeline)**: 2 weeks - Snakemake workflow + provenance  
**Phase 5 (Scaling)**: 1 week - HEALPix + optimization  
**Phase 6 (Final)**: 1 week - Production run + documentation  

**Total**: 8 weeks (expert), 3-4 months (realistic)

---

**Status**: ✅ Planning complete, ready for Phase 1 implementation  
**Next**: Start with `cosmo_lensing/io.py` - read_deflection_field()
