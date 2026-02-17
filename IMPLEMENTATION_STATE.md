# Implementation State: Julia → Python Conversion

**Status**: Phase 1 - Foundation (Not Started)  
**Last Updated**: 2026-02-17  
**Estimated Completion**: 8 weeks (expert) / 3-4 months (realistic)

---

## Summary

This codebase performs **weak gravitational lensing analysis** on cosmological N-body simulations (HAGN/Horizon-AGN). The pipeline:

1. **Reads** deflection fields (α₁, α₂) from Fortran binary files (20k×20k pixels per redshift slice)
2. **Computes** Jacobian matrix ∂αᵢ/∂xⱼ and second derivatives via finite differences
3. **Extracts** lensing observables: convergence (κ), shear (γ₁, γ₂), flexion (F₁, F₂, G₁, G₂)
4. **Measures** tangential correlations around galaxies (galaxy-galaxy lensing)
5. **Fits** NFW halo profiles to infer masses and concentrations
6. **Generates** publication plots for thesis

**Current State (Julia)**: 22 files, ~2500 duplicate lines, no tests, monolithic 1734-line toolkit, hardcoded paths, silent failures

**Target State (Python)**: 6 focused modules, comprehensive tests, production pipeline with Snakemake, HPC-ready, full reproducibility

---

## Issues and Proposed Solutions

### Issue #1: Monolithic 1734-line toolkit with no module boundaries
**Problem**: `raytrace_tk.jl` mixes I/O, derivatives, observables, statistics, plotting in 1734 lines  
**Solution**: Split into 6 Python modules with clear interfaces

**Tasks**:
- [ ] Create package structure: `cosmo_lensing/{__init__,io,cosmology,derivatives,observables,correlations,nfw}.py`
- [ ] Port I/O functions → `io.py` (read_deflection_field, write_fits_map, load_catalog)
- [ ] Port cosmology → `cosmology.py` (CosmologyCalculator class wrapping astropy)
- [ ] Port derivatives → `derivatives.py` (compute_jacobian with finite differences)
- [ ] Port observables → `observables.py` (convergence, shear, flexion functions)
- [ ] Port correlations → `correlations.py` (TangentialCorrelation class wrapping TreeCorr)
- [ ] Port NFW utilities → `nfw.py` (NFWProfile class, mass/concentration conversions)
- [ ] Add type hints to all function signatures
- [ ] Write module-level docstrings

---

### Issue #2: External C library dependency (libsoftlens.so) with no provenance
**Problem**: `core.jl` loads `libsoftlens.so` via FFI for cosmology calculations, no source/build instructions  
**Solution**: Replace with `astropy.cosmology.FlatLambdaCDM`

**Tasks**:
- [ ] Implement `CosmologyCalculator` class in `cosmology.py`:
  - [ ] `angular_diameter_distance(z1, z2=0)` → replaces `jsl_dda`
  - [ ] `critical_density(z)` → replaces hardcoded formula
  - [ ] `Omega_m(z)` → matter density parameter evolution
  - [ ] `distance_scaling(z)` → returns (arcsec2kpc, sigma_crit)
- [ ] Write validation test comparing astropy vs. libsoftlens outputs for 100 test redshifts
- [ ] Document expected numerical agreement (< 0.1%)
- [ ] Delete `core.jl` after validation passes

---

### Issue #3: No computational scaling strategy for large datasets
**Problem**: Cannot test on small patches, always processes full 20k×20k fields (12GB memory)  
**Solution**: Implement HEALPix hierarchical pixelization

**Tasks**:
- [ ] Create `spatial.py` module with `HEALPixDeflectionField` class
- [ ] Implement Cartesian → HEALPix conversion
- [ ] Implement patch extraction: `get_patch(center_ra, center_dec, radius_deg)`
- [ ] Define test hierarchy in `config.py`:
  - [ ] `unit`: NSIDE=16 (3k pixels), synthetic data, < 1 sec
  - [ ] `debug`: NSIDE=64 (49k pixels), real patch, < 10 sec
  - [ ] `validation`: NSIDE=512 (3M pixels), large subset, < 2 min
  - [ ] `production`: NSIDE=2048 (50M pixels), full data, ~10 min/slice
- [ ] Add `--config {unit,debug,validation,production}` CLI flag to all scripts

---

### Issue #4: Data processing coupled with visualization
**Problem**: `correlation_plots.jl` interleaves compute/stats/fitting/plotting, cannot reuse components  
**Solution**: Separate into compute/analyze/visualize layers

**Tasks**:
- [ ] Create `scripts/` directory for production scripts
- [ ] `compute_observables.py`: Read deflection → compute κ/γ/F/G → save FITS (no plotting)
- [ ] `compute_correlations.py`: Read observables + catalog → compute correlations → save .npz
- [ ] `fit_models.py`: Read correlations → fit NFW → save parameters CSV
- [ ] `generate_plots.py`: Read .npz + CSV → generate figures (separate from compute)
- [ ] Add argparse CLI to all scripts
- [ ] Enable headless execution (no X11 required)

---

### Issue #5: Massive DRY violation - 10 duplicate map generation scripts
**Problem**: `rebuild_{kappa,gamma1,gamma2,F1,F2,G1,G2,shear,flexion,rot}_maps.jl` are 95% identical (2500 duplicate lines)  
**Solution**: Single parameterized script with observable selection

**Tasks**:
- [ ] Implement `scripts/generate_maps.py` with CLI:
  - [ ] `--ids 50 100 150` (deflection IDs to process)
  - [ ] `--observables kappa gamma1 gamma2` (which to compute)
  - [ ] `--output-dir results/maps/` (where to save)
  - [ ] `--rebin-factor 5` (downsampling factor)
- [ ] Define `OBSERVABLE_MAP = {'kappa': observables.convergence, 'gamma1': observables.shear_1, ...}`
- [ ] Delete all 10 `rebuild_*_maps.jl` files after testing
- [ ] Add `--observables all` shortcut

---

### Issue #6: Zero error handling in I/O operations
**Problem**: `read_bin()` silently returns `(nothing, nothing)` on error, crashes propagate downstream  
**Solution**: Explicit exceptions with validation

**Tasks**:
- [ ] Define custom exceptions in `io.py`:
  - [ ] `DeflectionFieldError(Exception)` for invalid/corrupted files
- [ ] Implement `read_deflection_field()` with checks:
  - [ ] Raise `FileNotFoundError` if file missing
  - [ ] Raise `DeflectionFieldError` if size < 1000 bytes (corrupted)
  - [ ] Raise `DeflectionFieldError` if dimensions invalid (≤ 0)
  - [ ] Raise `DeflectionFieldError` if NaN/Inf in data
  - [ ] Add logging (info/warning/error levels)
- [ ] Add try/except at call sites with explicit handling:
  - [ ] Log warning and skip on FileNotFoundError (missing optional data)
  - [ ] Log error and raise on DeflectionFieldError (fail fast)
- [ ] Write tests for error conditions: `test_read_missing_file()`, `test_read_corrupted_file()`

---

### Issue #7: Undocumented mathematical formulas with no validation
**Problem**: Flexion formulas (F₁, F₂, G₁, G₂) lack citations, no tests vs. analytic solutions  
**Solution**: Add docstrings with paper references, implement synthetic data generators, write unit tests

**Tasks**:
- [ ] Add citations to all physics functions in `observables.py`:
  - [ ] `convergence()` → Bartelmann & Schneider (2001) Eq. 3.47
  - [ ] `shear_1()`, `shear_2()` → Bartelmann & Schneider (2001) Eq. 3.48
  - [ ] `flexion_F1()`, `flexion_F2()` → Bacon et al. (2006) MNRAS 365, 414, Eq. 8
  - [ ] `flexion_G1()`, `flexion_G2()` → Bacon et al. (2006), Eq. 9
- [ ] Create `synthetic_data.py` module:
  - [ ] `generate_point_mass_deflection(theta_E, grid_size)` → analytic solution
  - [ ] `generate_sis_deflection(sigma_v, grid_size)` → Singular Isothermal Sphere
  - [ ] `generate_nfw_deflection(mass, c, z, grid_size)` → NFW from Wright & Brainerd (2000)
- [ ] Write unit tests in `tests/test_observables.py`:
  - [ ] `test_point_mass_convergence()` → compare to κ(θ) = θ_E²/(2θ²)
  - [ ] `test_nfw_convergence_normalization()` → ∫ κ(r) 2πr dr = enclosed mass
  - [ ] `test_flexion_symmetry()` → F_radial ≈ 0 for spherical lens
  - [ ] `test_derivatives_finite_difference()` → compare to analytic derivatives
- [ ] Target: 10+ analytic validation tests, all passing

---

### Issue #8: 50KB dead code from abandoned SIS model attempts
**Problem**: `correlation_plots_sis_model.jl` and `_2.jl` are unused (2 abandoned experiments)  
**Solution**: Delete dead code, refactor correlation function into reusable class

**Tasks**:
- [ ] DELETE `correlation_plots_sis_model.jl` (25KB)
- [ ] DELETE `correlation_plots_sis_model_2.jl` (25KB)
- [ ] Implement `TangentialCorrelation` class in `correlations.py`:
  - [ ] Wrap TreeCorr's `NGCorrelation` for standardization
  - [ ] `__init__(rmin, rmax, nbins, sep_units, var_method='jackknife')`
  - [ ] `compute(ra_lens, dec_lens, ra_src, dec_src, g1, g2, weights=None)` → returns dict
  - [ ] Return: `{'r': radii, 'xi': correlation, 'xi_err': errors, 'npairs': counts}`
- [ ] Implement `NFWProfile` class in `nfw.py`:
  - [ ] `__init__(ks, rs)` → convergence scale, scale radius
  - [ ] `convergence(r)` → κ(r) from Wright & Brainerd (2000) Eq. 11
  - [ ] `shear_tangential(r)` → γ_t(r) from Wright & Brainerd (2000) Eq. 12
  - [ ] Static method `_F_function(x)` → auxiliary function Eq. 13
- [ ] Test against current `comp_gs_corr()` output (should match within 1%)

---

### Issue #9: No test mode for small patches (development velocity)
**Problem**: Debugging requires full 10-minute production runs, no fast iteration  
**Solution**: Hierarchical test pyramid with synthetic + real data at multiple scales

**Tasks**:
- [ ] Define `TEST_CONFIGS` in `config.py`:
```python
{
  'unit': {'npix': 100, 'ngal': 5, 'nslices': 1, 'data': 'synthetic', 'runtime': '<1s'},
  'debug': {'npix': 1000, 'ngal': 50, 'nslices': 2, 'data': 'real_patch', 'runtime': '<10s'},
  'validation': {'npix': 5000, 'ngal': 500, 'nslices': 5, 'data': 'real', 'runtime': '<2min'},
  'production': {'npix': 20000, 'ngal': None, 'nslices': 10, 'data': 'real', 'runtime': '~10min/slice'}
}
```
- [ ] Implement synthetic data generators (see Issue #7)
- [ ] Implement patch extraction for real data:
  - [ ] `extract_patch(deflection_full, center_ra, center_dec, patch_size)` in `spatial.py`
- [ ] Add `--config {unit,debug,validation,production}` to all scripts
- [ ] Write pytest fixtures using `@pytest.fixture(params=['unit', 'debug'])`
- [ ] Test full pipeline on unit scale (should complete in < 1 sec)

---

### Issue #10: No validation against established weak lensing software
**Problem**: Custom correlation implementation never compared to TreeCorr (gold standard)  
**Solution**: Cross-validation tests requiring < 1% agreement

**Tasks**:
- [ ] Install TreeCorr: `pip install treecorr` (version ≥ 4.3)
- [ ] Implement `tests/test_against_treecorr.py`:
  - [ ] `test_tangential_shear_vs_treecorr()` → load synthetic NFW, compare outputs
  - [ ] `test_convergence_azimuthal_vs_treecorr()` → test κ(r) averaging
  - [ ] Use `np.testing.assert_allclose(our_result, treecorr_result, rtol=0.01)`
- [ ] Generate test data:
  - [ ] Synthetic NFW deflection (M=10¹⁴ Msun, c=5, z=0.5)
  - [ ] 100 random galaxy positions
  - [ ] Compute γ_t(r) both ways
- [ ] If discrepancy > 1%:
  - [ ] Debug: check tangential projection formula
  - [ ] Debug: check binning scheme (log vs linear)
  - [ ] Debug: check coordinate systems (RA/Dec vs Cartesian)
- [ ] CRITICAL: Must pass before production use

---

### Issue #11: Unclear goal (research exploration vs. production pipeline)
**Problem**: Code exhibits identity crisis - half interactive, half automated  
**Solution**: Commit to production pipeline with Snakemake for reproducibility

**Tasks**:
- [ ] Create `workflow/Snakefile` with rules:
  - [ ] `rule compute_observables`: deflection.bin → kappa.fits, gamma.fits, ...
  - [ ] `rule compute_correlations`: observables + catalog → tangential_shear.npz
  - [ ] `rule fit_models`: tangential_shear.npz → nfw_params.csv
  - [ ] `rule generate_plots`: nfw_params.csv → correlation_plots.png
  - [ ] `rule all`: aggregate all outputs
- [ ] Create `workflow/config.yaml`:
```yaml
deflection_ids: [50, 100, 150, 200, 250, 300, 350, 400, 417, 450]
observables: ['kappa', 'gamma1', 'gamma2', 'F1', 'F2', 'G1', 'G2']
galaxy_catalog: "Data/Galaxies_0-6_lensed.v2.0_cut_i27.fits"
correlation: {rmin: 0.05, rmax: 5.0, nbins: 50, var_method: 'jackknife'}
output_dir: "results/run_{timestamp}/"
```
- [ ] Add provenance tracking:
  - [ ] Save git commit hash to `{output_dir}/git_commit.txt`
  - [ ] Save config snapshot to `{output_dir}/config.yaml`
  - [ ] Compute SHA256 checksums of outputs → `{output_dir}/checksums.txt`
- [ ] Test workflow:
  - [ ] `snakemake --configfile workflow/config.yaml --dry-run` (validate DAG)
  - [ ] `snakemake --config scale=debug --cores 4` (local test)
  - [ ] `snakemake --config scale=production --profile slurm --jobs 50` (HPC)

---

### Issue #12: Inefficient jackknife implementation (80% of runtime)
**Problem**: Current "jackknife" is actually bootstrap with 5000 samples, very slow  
**Solution**: Use TreeCorr's spatial jackknife (50 patches, 100× faster)

**Tasks**:
- [ ] Modify `TangentialCorrelation` class to enable jackknife:
  - [ ] Pass `var_method='jackknife'` to `treecorr.NGCorrelation(...)`
  - [ ] TreeCorr handles spatial partitioning automatically
  - [ ] Errors include spatial correlations (proper treatment)
- [ ] Implement covariance matrix export (for proper χ² fitting):
```python
def compute_covariance_matrix(ng_correlation, method='jackknife'):
    """Extract full covariance matrix from TreeCorr NGCorrelation object."""
    # TreeCorr provides varxi (diagonal), can also bootstrap for full cov
    if method == 'jackknife':
        # Use jackknife samples to compute cov(xi_i, xi_j)
        ...
    return cov_matrix  # (nbins, nbins)
```
- [ ] Update fitting code to use covariance:
```python
def fit_nfw_with_covariance(r, xi, cov, initial_params):
    """Fit NFW with proper χ² = (data - model)ᵀ Cov⁻¹ (data - model)"""
    cov_inv = np.linalg.inv(cov)
    def chi2(params):
        model = nfw_profile(r, *params)
        residual = xi - model
        return residual @ cov_inv @ residual
    result = scipy.optimize.minimize(chi2, initial_params, bounds=[(0.001, 10), (10, 1000)])
    return result
```
- [ ] Remove old jackknife loop from Julia equivalent (lines 296-346 in correlation_plots.jl)
- [ ] Validate: Compare error bars (old bootstrap vs. new spatial jackknife)
- [ ] Expected: ~100× speedup (5000 samples → 50 patches)

---

## Test Strategy

Tests must be **atomic** (single feature), **fast** (suite < 30 sec), and **deterministic** (no flaky tests).

### Unit Tests (< 10 seconds total)
Run on every commit via CI/CD. Use synthetic data (100×100 pixels).

**I/O Module** (`tests/test_io.py`):
- [ ] `test_read_deflection_field_success()` → reads valid file, returns array + redshift
- [ ] `test_read_deflection_field_missing()` → raises FileNotFoundError
- [ ] `test_read_deflection_field_corrupted()` → raises DeflectionFieldError
- [ ] `test_write_fits_map()` → writes FITS, header has WCS keywords + redshift

**Cosmology Module** (`tests/test_cosmology.py`):
- [ ] `test_angular_diameter_distance()` → compare to astropy directly (< 0.1% diff)
- [ ] `test_critical_density()` → check units, compare to literature values
- [ ] `test_distance_scaling()` → returns (arcsec2kpc, sigma_crit) with correct units

**Derivatives Module** (`tests/test_derivatives.py`):
- [ ] `test_jacobian_point_mass()` → analytic derivatives for point mass
- [ ] `test_jacobian_determinant_positive()` → no multiple images (det > 0 everywhere)
- [ ] `test_finite_difference_accuracy()` → compare to analytic derivatives (2nd order accurate)

**Observables Module** (`tests/test_observables.py`):
- [ ] `test_convergence_point_mass()` → κ(θ) = θ_E²/(2θ²) for point mass
- [ ] `test_convergence_sis()` → κ(θ) = θ_E/(2θ) for SIS
- [ ] `test_nfw_convergence_normalization()` → ∫ κ(r) 2πr dr = M_enclosed
- [ ] `test_shear_rotation_invariance()` → γ(θ+π) = -γ(θ) for point mass
- [ ] `test_flexion_F_symmetry()` → F_radial ≈ 0 for spherical lens
- [ ] `test_flexion_G_symmetry()` → G(φ) = G(-φ) for spherical lens

**NFW Module** (`tests/test_nfw.py`):
- [ ] `test_nfw_F_function_limits()` → F(x→1) = 1, F(x→0) ∝ 1/x, F(x→∞) → 0
- [ ] `test_mass_concentration_roundtrip()` → (M,c) → (ks,rs) → (M,c) exact
- [ ] `test_nfw_convergence_vs_wright2000()` → compare to paper formulas

**Pass Criteria**: All tests pass, > 80% code coverage (`pytest --cov=cosmo_lensing`)

---

### Integration Tests (< 2 minutes)
Run before merging to main. Use 1000×1000 real data patch OR large synthetic field.

**Full Pipeline** (`tests/test_integration.py`):
- [ ] `test_deflection_to_observables()` → read deflection → compute all observables → save FITS
- [ ] `test_observables_to_correlation()` → load observables + catalog → compute γ_t(r)
- [ ] `test_correlation_to_fit()` → load correlation → fit NFW → return (M, c, χ²)
- [ ] `test_end_to_end_pipeline()` → deflection → observables → correlation → fit → plot

**TreeCorr Validation** (`tests/test_against_treecorr.py`):
- [ ] `test_tangential_shear_vs_treecorr()` → must agree within 1% (CRITICAL)
- [ ] `test_convergence_azimuthal_vs_treecorr()` → must agree within 1%
- [ ] If > 1% discrepancy: FAIL and investigate before proceeding

**Error Propagation** (`tests/test_errors.py`):
- [ ] `test_jackknife_error_scaling()` → σ(xi) ∝ 1/√N_patches (expected behavior)
- [ ] `test_covariance_matrix_positive_definite()` → all eigenvalues > 0

**Pass Criteria**: All integration tests pass, correlations smooth (no NaN/artifacts)

---

### Validation Tests (< 30 minutes)
Run before production, weekly during development. Use 5000×5000 pixels, 500 galaxies, 5 redshift slices.

**Regression Tests** (`tests/test_regression.py`):
- [ ] `test_correlation_vs_saved_reference()` → compare to saved .npz files (rtol=1e-6)
- [ ] `test_fitted_masses_vs_reference()` → compare to saved CSV (rtol=0.01)
- [ ] Store reference outputs in `tests/reference_outputs/` (checked into git)

**Comparison to Julia Code**:
- [ ] Run both pipelines on identical input (5 redshift slices)
- [ ] Compare correlation functions: must match < 5%
- [ ] Compare fitted masses: must match < 10% (systematic OK if documented)
- [ ] If > 10% discrepancy: investigate before finalizing

**Visual Inspection**:
- [ ] Generate plots for all observables (κ, γ, F, G)
- [ ] Check for: NaN values, edge artifacts, smoothness
- [ ] Compare to Julia plots side-by-side

**Pass Criteria**:
- All validation tests pass
- Correlation functions match Julia < 5%
- Fitted masses match Julia < 10%
- No visual artifacts in plots

---

### Production Tests (One-time before thesis submission)
Full dataset: 20k×20k pixels, all galaxies, 10 redshift slices.

**Completeness**:
- [ ] All expected outputs generated (10 redshift × 7 observables = 70 FITS files)
- [ ] No missing files, no corrupted FITS headers
- [ ] All correlation .npz files present
- [ ] All NFW fit CSVs present
- [ ] All plots present

**Reproducibility**:
- [ ] Run pipeline twice with identical config → identical outputs (SHA256 match)
- [ ] Re-run after code change → only affected outputs regenerated (Snakemake dependency tracking)

**Scientific Validation**:
- [ ] Fitted masses reasonable (10¹³ - 10¹⁵ Msun for galaxy halos)
- [ ] Reduced χ² ≈ 1 for NFW fits (good fit quality)
- [ ] Error bars reasonable (not too small, not too large)
- [ ] Trends with redshift/mass make physical sense

**Pass Criteria**:
- Pipeline completes without errors
- All outputs present and valid
- Results reproducible within 0.1%
- Thesis figures generated from outputs
- External user can reproduce from scratch (tested by advisor/colleague)

---

## Key Dependencies

**Core Scientific Stack**:
- `numpy >= 1.24` - Numerical arrays
- `scipy >= 1.10` - Curve fitting, interpolation, optimization
- `astropy >= 5.0` - Cosmology, FITS I/O, units, coordinates
- `matplotlib >= 3.6` - Publication plots

**Weak Lensing (CRITICAL)**:
- `treecorr >= 4.3` - Gold standard correlation functions (DES/HSC/LSST use this)
- `healpy >= 1.16` - HEALPix spatial pixelization

**Workflow & Testing**:
- `snakemake >= 7.0` - Pipeline orchestration
- `pytest >= 7.0` - Testing framework
- `hypothesis >= 6.0` - Property-based testing (optional)

**Installation**:
```bash
pip install numpy scipy astropy matplotlib treecorr healpy snakemake pytest
```

---

## Implementation Phases

### Phase 1: Foundation (Weeks 1-2) - **CURRENT PHASE**
**Goal**: Core modules functional, basic tests passing

**Deliverables**:
- [ ] Package structure created (`cosmo_lensing/` with 7 modules)
- [ ] I/O module working (read/write with error handling)
- [ ] Cosmology module working (validated vs astropy)
- [ ] Derivatives module working (Jacobian computation)
- [ ] Observables module working (κ, γ, F, G with citations)
- [ ] Pytest infrastructure set up (fixtures, config)
- [ ] 20+ unit tests passing

**Test**: Process one deflection field, compute all observables, save FITS maps

---

### Phase 2: Synthetic Data & Validation (Week 3)
**Goal**: Test suite comprehensive, validated against analytic solutions

**Deliverables**:
- [ ] Synthetic data generators (point mass, SIS, NFW)
- [ ] 10+ analytic validation tests (all passing)
- [ ] Test hierarchy implemented (unit/debug/validation/production)
- [ ] All tests pass at unit and debug scales

**Test**: `pytest tests/ --config unit` completes in < 10 seconds, all passing

---

### Phase 3: Correlation Functions (Week 4)
**Goal**: TreeCorr integration, validation < 1%

**Deliverables**:
- [ ] `TangentialCorrelation` class wrapping TreeCorr
- [ ] `NFWProfile` class with lensing functions
- [ ] TreeCorr validation tests (< 1% agreement achieved)
- [ ] Fitting code with covariance matrix support

**Test**: Synthetic NFW correlation matches analytic γ_t(r) from Wright & Brainerd (2000)

---

### Phase 4: Production Pipeline (Weeks 5-6)
**Goal**: Snakemake workflow, batch processing, reproducibility

**Deliverables**:
- [ ] Snakefile with rules (compute_observables, compute_correlations, fit_models, generate_plots)
- [ ] config.yaml template
- [ ] Production scripts with CLI (compute_observables.py, etc.)
- [ ] Provenance tracking (git commit, config snapshot, checksums)
- [ ] Workflow tested at debug and validation scales

**Test**: `snakemake --config scale=validation` completes successfully, outputs reproducible

---

### Phase 5: Scaling & Optimization (Week 7)
**Goal**: HEALPix support, handle production scale efficiently

**Deliverables**:
- [ ] `HEALPixDeflectionField` class
- [ ] Cartesian ↔ HEALPix conversion
- [ ] Patch extraction working
- [ ] Parallel processing (multiprocessing for embarrassingly parallel tasks)
- [ ] Memory optimization (chunked processing, cleanup)
- [ ] Can process NSIDE=2048 (50M pixels)

**Test**: Production scale run (one redshift slice) completes in < 1 hour

---

### Phase 6: Final Validation (Week 8)
**Goal**: Production run, comparison to Julia, documentation

**Deliverables**:
- [ ] Full production run (10 redshift slices) completed
- [ ] Comparison to Julia code (< 5% correlation, < 10% masses)
- [ ] Comprehensive documentation (README, API docs, tutorial)
- [ ] Code cleanup (delete Julia files, linting with black/isort)
- [ ] Release v1.0 (git tag, Zenodo DOI)

**Test**: External user (advisor/colleague) reproduces thesis results from scratch

---

## Current Status: Phase 1 - Foundation

**Next Actions**:
1. Create package structure: `mkdir -p cosmo_lensing tests scripts workflow`
2. Set up Python environment: `python -m venv venv && source venv/bin/activate`
3. Install dependencies: `pip install numpy scipy astropy matplotlib treecorr healpy pytest snakemake`
4. Start implementing `cosmo_lensing/io.py` → `read_deflection_field()` function
5. Reference Julia code: `raytrace_tk.jl` lines 77-110

**Blocking Issues**: None  
**Dependencies Met**: All (Python 3.10+, standard packages available)

---

## Success Criteria

Must achieve **ALL** of the following:

**Technical**:
- [ ] All 12 issues resolved (no exceptions)
- [ ] 100% of Julia functionality ported
- [ ] Test suite > 80% coverage, all passing
- [ ] Pipeline 10× faster than Julia (TreeCorr optimization)
- [ ] Memory efficient (32GB laptop for debug scale)
- [ ] Fully reproducible (same inputs → same outputs within 0.1%)

**Scientific**:
- [ ] TreeCorr validation < 1% difference (CRITICAL)
- [ ] Fitted masses match Julia < 10% (systematic documented if present)
- [ ] Error bars reasonable (spatial jackknife accounts for correlations)
- [ ] Thesis results reproducible by external user
- [ ] Code suitable for publication (Zenodo DOI)

**Engineering**:
- [ ] Clean modular architecture (6 modules, single responsibility)
- [ ] No external C dependencies (pure Python + standard libs)
- [ ] Explicit error handling (no silent failures)
- [ ] Comprehensive documentation (README + API docs + docstrings)
- [ ] Production pipeline (Snakemake, provenance, HPC-ready)
- [ ] Maintainable (future student can extend)

---

**Last Updated**: 2026-02-17  
**Status**: Phase 1 - Not Started  
**Estimated Completion**: 2026-04-17 (8 weeks) or 2026-06-17 (4 months)
