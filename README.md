# Cosmo Lensing: Weak Gravitational Lensing Analysis

Python package for weak gravitational lensing analysis on cosmological N-body simulations (HAGN/Horizon-AGN).

## Overview

This package provides a complete pipeline for:
- Reading deflection fields from Fortran binary files
- Computing lensing Jacobian matrices with finite differences
- Extracting lensing observables: convergence (κ), shear (γ), flexion (F, G)
- Measuring tangential correlations (galaxy-galaxy lensing)
- Fitting NFW halo profiles
- Generating publication-quality plots

**Status**: Phases 1-4 complete - production pipeline ready for deployment.

## Installation

```bash
# Clone repository
git clone <repository-url>
cd laurant-thesis

# Create virtual environment
python3.12 -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Install package in development mode
pip install -e .

# Run tests
pytest tests/
```

## Quick Start

### Using Python API

```python
from cosmo_lensing import io, derivatives, observables

# Read deflection field
alpha, redshift = io.read_deflection_field("deflection.bin")

# Compute Jacobian and second derivatives
jac = derivatives.compute_jacobian(
    alpha,
    pixel_scale=1.0,
    compute_second_derivatives=True
)

# Compute observables
obs = observables.compute_all_observables(jac)

# Access results
kappa = obs['kappa']
gamma1 = obs['gamma1']
gamma2 = obs['gamma2']
```

### Using Production Scripts

```bash
# Compute lensing observables
python scripts/compute_observables.py deflection_050.bin \
    --outdir results/ \
    --observables kappa gamma1 gamma2 \
    --overwrite

# Compute tangential shear correlations
python scripts/compute_correlations.py \
    --gamma1 results/gamma1.fits \
    --gamma2 results/gamma2.fits \
    --catalog galaxies.fits \
    --outdir results/ \
    --rmin 0.05 --rmax 5.0 --nbins 20
```

### Using Snakemake Workflow

```bash
# Dry run to show planned jobs
snakemake -n

# Run locally with 4 cores
snakemake --cores 4

# Run on HPC with SLURM
snakemake --profile slurm --jobs 50
```

## Modules

### Core Computation (`cosmo_lensing/`)

- **`io.py`** (268 lines): Input/output for deflection fields and FITS maps
  - Read Fortran binary deflection fields
  - Write FITS maps with WCS headers
  - Galaxy catalog loading

- **`cosmology.py`** (218 lines): Cosmological calculations
  - Replaced C library (libsoftlens.so) with pure Python
  - Uses `astropy.cosmology.FlatLambdaCDM`
  - Angular diameter distance, critical density, comoving distance

- **`derivatives.py`** (235 lines): Finite-difference derivatives
  - 4th-order accurate central differences
  - 3rd-order one-sided differences at boundaries
  - Extended Jacobian with second derivatives for flexion

- **`observables.py`** (311 lines): Lensing observables
  - Convergence κ, shear γ₁/γ₂
  - First flexion F₁/F₂, second flexion G₁/G₂
  - Rotation field (diagnostic)
  - Literature citations (Bartelmann & Schneider 2001, Bacon et al. 2006)

- **`nfw.py`** (280 lines): NFW halo profiles
  - Analytic convergence κ(r) and tangential shear γₜ(r)
  - Auxiliary functions from Wright & Brainerd (2000)
  - Excess surface density ΔΣ(r)

- **`correlations.py`** (218 lines): Tangential correlation functions
  - TreeCorr wrapper for galaxy-galaxy lensing
  - Automatic npatch adjustment for small samples
  - Jackknife/bootstrap error estimation
  - Flexion correlations

### Production Scripts (`scripts/`)

- **`compute_observables.py`** (300 lines): CLI for observable maps
  - Reads deflection binaries
  - Computes Jacobian and extracts observables
  - Writes FITS files with WCS headers

- **`compute_correlations.py`** (350 lines): CLI for correlation functions
  - Reads FITS maps and galaxy catalogs
  - Computes tangential shear profiles
  - Writes NPZ files with correlation results

### Workflow (`workflow/`)

- **`Snakefile`** (250 lines): Snakemake workflow definition
  - Parallelizable rules for batch processing
  - Provenance tracking (git commit, checksums)
  - HPC-ready (SLURM profile support)

- **`config.yaml`** (65 lines): Analysis configuration
  - Field IDs, observables, correlation parameters
  - Performance tuning options

## Testing

**104 tests, 100% passing, < 4s runtime**

- `test_cosmology.py`: 15 tests for cosmological calculations
- `test_derivatives.py`: 16 tests for finite differences
- `test_observables.py`: 21 tests for lensing observables
- `test_synthetic_data.py`: 18 tests for analytic validation
- `test_nfw.py`: 22 tests for NFW profile
- `test_correlations.py`: 8 tests for TreeCorr integration
- `test_integration.py`: 1 end-to-end pipeline test
- `test_pipeline.py`: 3 production script tests

Run tests:
```bash
# All tests
pytest tests/ -v

# Fast tests only (exclude slow pipeline tests)
pytest tests/ -v -m "not slow"

# With coverage report
pytest tests/ --cov=cosmo_lensing --cov-report=html
```

## Performance

- **Jacobian computation**: ~1s for 100×100 field, ~30s for 1000×1000 field
- **Full pipeline**: ~5s for small test field (150×150 with 30 lenses)
- **Memory**: ~16 GB for 20k×20k deflection fields
- **Parallelization**: Embarrassingly parallel over field IDs (use Snakemake)

## Implementation Status

### ✅ Phase 1: Foundation (Complete)
- Core computational modules
- Comprehensive unit tests
- No external C dependencies

### ✅ Phase 2: Synthetic Data (Complete)
- Analytic test data generators (point mass, SIS, NFW, pure shear/convergence)
- Validation against analytic solutions
- Finite-difference accuracy tests

### ✅ Phase 3: Correlation Functions (Complete)
- Full NFW profile implementation
- TreeCorr integration for tangential correlations
- Flexion correlation support
- Jackknife/bootstrap error estimation

### ✅ Phase 4: Production Pipeline (Complete)
- CLI scripts for observable computation and correlations
- Snakemake workflow for batch processing
- End-to-end integration tests
- Provenance tracking

### ⏸️ Phase 5: HEALPix & Scaling (Pending)
- HEALPix pixel operations
- Memory-efficient chunking
- Performance optimization (numba JIT)

### ⏸️ Phase 6: Final Validation (Pending)
- Comparison to original Julia code
- Documentation (Sphinx, tutorial notebooks)
- Package release (Zenodo DOI)

## Code Statistics

- **Total**: ~5000 lines of Python code
- **Package modules**: 1800 lines
- **Test suite**: 1400 lines
- **Production scripts**: 650 lines
- **Workflow**: 315 lines
- **Test coverage**: >90%

## Features

- ✅ No external C dependencies (replaced libsoftlens.so with astropy)
- ✅ Type hints and comprehensive docstrings throughout
- ✅ Explicit error handling with custom exceptions
- ✅ Citations to literature in code comments
- ✅ Logging at multiple levels (DEBUG, INFO, WARNING, ERROR)
- ✅ 4th-order accurate finite differences
- ✅ Full test coverage with fast execution
- ✅ Production-ready CLI tools
- ✅ Reproducible workflow with Snakemake

## References

- Bartelmann & Schneider (2001): "Weak gravitational lensing", *Physics Reports* 340, 291-472
- Bacon et al. (2006): "Measuring gravitational spin", *MNRAS* 365, 414-428
- Wright & Brainerd (2000): "Gravitational lensing by NFW halos", *ApJ* 534, 34-40
- Jarvis et al. (2016): "TreeCorr: correlation functions with tree codes", *MNRAS* 460, 2245-2281

## License

Academic/Research use

## Author

Laurent Magri-Stella  
Python conversion: Comprehensive testing, modular architecture, production pipeline
