# Implementation Progress Checklist

Track your progress through the modernization roadmap. Check off items as completed.

---

## Phase 0: Setup (Days 1-2)

**Goal**: Initialize project structure, testing framework, CI/CD

### Repository & Dependencies
- [ ] Create `pyproject.toml` (Poetry or setuptools)
  - [ ] Package name: `lensing-pipeline`
  - [ ] Python >=3.9
  - [ ] Core deps: numpy, scipy, astropy, treecorr, h5py, numba
  - [ ] Dev deps: pytest, black, mypy, sphinx (optional)
- [ ] Create `.gitignore` (data, venv, __pycache__, *.pyc)
- [ ] Create `requirements.txt` or use Poetry lock file
- [ ] Initial commit to git

### Project Structure
- [ ] Create `/lensing` package directory with `__init__.py`
- [ ] Create `/tests` with `conftest.py`
- [ ] Create `/scripts` directory
- [ ] Create `/data/{config,raw,processed}` directories
- [ ] Create `/docs` directory

### Testing & Quality
- [ ] `pytest.ini`: pytest configuration
- [ ] `.pre-commit-config.yaml`: black, isort, mypy, flake8
- [ ] `setup.cfg` or `pyproject.toml`: flake8 + mypy config
- [ ] GitHub Actions: `.github/workflows/tests.yml`
  - [ ] Run on: push, pull_request
  - [ ] Matrix: Python 3.9, 3.10, 3.11
  - [ ] Jobs: pytest, mypy, coverage

### Documentation
- [ ] `/docs/conf.py`: Sphinx config (if using)
- [ ] `/docs/index.rst`: TOC
- [ ] `README.md`: Project overview + installation
- [ ] `LICENSE` file (MIT or similar)

### Verify Phase 0
```bash
pytest tests/ -v
# → PASSED (0 tests, that's OK for now)

mypy lensing/
# → Success: no issues found

pip install -e .
# → Successfully installed lensing-pipeline
```

---

## Phase 1: Core Physics (Weeks 1-2)

### 1a: Cosmology & Configuration (Days 3-5)

**Owner**: Person A  
**Goal**: Centralized config system + astropy cosmology wrapper

#### Cosmology Module
- [ ] Create `lensing/cosmology.py`
  - [ ] `CosmologyCalculator` class
  - [ ] `angular_diameter_distance(z1, z2)` → backing astropy
  - [ ] `comoving_distance(z)`
  - [ ] `critical_surface_density(zl, zs)` or `sigma_crit()`
  - [ ] `hubble_parameter(z)` → E(z)
  - [ ] Full LaTeX docstrings for each function
  - [ ] Input validation (z ranges, etc.)

#### Configuration Module
- [ ] Create `lensing/config.py`
  - [ ] `CosmologyConfig` dataclass (h, omega_m, omega_lambda, z_lens)
  - [ ] `PipelineConfig` dataclass (z_ids, data_paths, etc.)
  - [ ] YAML loader: `load_config(filepath: str) -> dict`
  - [ ] Create sample config: `data/config/default.yaml`

#### Tests
- [ ] Create `tests/unit/test_cosmology.py`
  - [ ] `test_angular_diameter_distance()` → known values (PLANCK 2015)
  - [ ] `test_aad_symmetry()`: D_A(z1,z) + D_A(z,z2) relation
  - [ ] `test_comoving_distance_zeros()`: D_C(0) = 0
  - [ ] `test_critical_density_limit()`: Σ_crit → ∞ as z_s → z_l
  - [ ] `test_hubble_parameter_bounds()`: E(z) ≥ 1 for all z
  - [ ] 10-15 tests total

- [ ] Create `tests/unit/test_config.py`
  - [ ] `test_config_load_yaml()`: can parse YAML
  - [ ] `test_config_validation()`: omega_m + omega_lambda ≤ 1
  - [ ] `test_dataclass_fields()`: all expected fields present

#### Deliverables
- [ ] `lensing/cosmology.py` ✓
- [ ] `lensing/config.py` ✓
- [ ] `data/config/default.yaml` ✓
- [ ] `tests/unit/test_cosmology.py` ✓ (15+ tests)
- [ ] `tests/unit/test_config.py` ✓ (5+ tests)
- [ ] Full docstrings with paper citations

#### Checkpoint
```bash
pytest tests/unit/test_cosmology.py tests/unit/test_config.py -v
# → All pass (20+ tests)
```

---

### 1b: NFW Profile Consolidation (Days 5-9)

**Owner**: Person A  
**Goal**: Single source of truth for all NFW lensing functions

#### NFW Module
- [ ] Create `lensing/lensing_models/nfw_profile.py`
- [ ] `NFWProfile` class with static methods:
  - [ ] `delta_c(c: float) -> float`
    - Formula + citation (Navarro+ 1997, Eq. 1)
    - Range check: c > 0
  - [ ] `ksrs2mc(ks, rs, q, cosmo, z) -> (M200, c)`
    - Full formula + units consistent
    - Input validation (ks, rs > 0), etc.
  - [ ] `mc2ksrs(m, c, q, cosmo, z) -> (ks, rs)`
    - Inverse of above; round-trip test
  - [ ] `c_m_K16(m, z) -> c` (concentration from K16 formula)
    - Interpolation setup
  - [ ] Flexion functions: `Flens(r, rs)`, `Fprime()`, `Hlens()`, `Glens()`
    - All with LaTeX docstrings + references (Bartelmann & Schneider 2001)
  - [ ] `kappafunc(r, ks, rs)` → convergence profile
  - [ ] All functions vectorized (accept arrays)

#### Tests
- [ ] Create `tests/unit/test_nfw_profile.py` (20-30 tests)
  - [ ] `test_delta_c_value(c=4)` → precomputed from paper
  - [ ] `test_mc2ksrs_roundtrip()`: mc2ksrs → ksrs2mc → recover input
  - [ ] `test_ksrs2mc_roundtrip()`: reverse
  - [ ] `test_delta_c_limits()`: δ_c(c→0) behavior, δ_c(c→∞)
  - [ ] `test_c_m_interpolation()`: c_m_K16 matches table values
  - [ ] `test_kappafunc_positive()`: κ > 0 for valid inputs
  - [ ] `test_kappafunc_units_consistent()`: κ [dimensionless] from ks [shear], rs [physical]
  - [ ] Closed-form validation: δ_c(1.0) = 1/3 * some_value
  - [ ] Array input tests: all functions handle np.ndarray

#### Validation Against Julia
- [ ] Compare output with Julia `nfw_tk.jl` functions on test inputs
- [ ] Document any formula differences (if any)

#### Deliverables
- [ ] `lensing/lensing_models/nfw_profile.py` ✓
- [ ] `lensing/lensing_models/__init__.py` ✓ (exports)
- [ ] `tests/unit/test_nfw_profile.py` ✓ (25+ tests)
- [ ] Every function has LaTeX docstring + min 2 paper citations

#### Checkpoint
```bash
pytest tests/unit/test_nfw_profile.py -v
# → All pass (25+ tests)
# → No Julia value mismatches
```

---

### 1c: Synthetic Test Data & Julia Reference (Days 9-12)

**Owner**: Person B  
**Goal**: Create test fixtures; generate Julia baseline for regression testing

#### Synthetic Map Generators
- [ ] Create `tests/fixtures/synth_maps.py`
  - [ ] `@pytest.fixture simple_constant_map()` → uniform α = [x, y]
    - 64×64, simple geometry, Δ ≈ 0 (nearly identity Jacobian)
  - [ ] `@pytest.fixture synthetic_nfw_map(seed=42, **kwargs)` → single NFW halo at center
    - Returns: (alpha, redshift)
    - Parameters: size, rs (scale radius), ks (shear norm), seed
    - Validatable: exact known κ, γ from formula
  - [ ] `@pytest.fixture synthetic_gaussian_map()` → smooth random field
    - For robustness testing
  - [ ] `@pytest.fixture real_corner_patch()` → 64×64 extract from real FITS
    - Load `Data/Galaxies_0-6_lensed.v2.0_cut_i27.fits`
    - Extract corner [0:64, 0:64, :]
    - One fixture per z-snapshot (id 0050, 0100, etc.)

#### Julia Reference Outputs
- [ ] Create Julia script: `scripts/generate_julia_reference.jl`
  - [ ] Load `synthetic_nfw_map(seed=42)` (Python → save to FITS, read in Julia)
  - [ ] Or: hard-code synthetic_nfw_map in Julia
  - [ ] Compute:
    - [ ] jac = alpha2jac(alpha)
    - [ ] kappa = jac2kappa(jac)
    - [ ] gamma1 = jac2gamma1(jac)
    - [ ] gamma2 = jac2gamma2(jac)
    - [ ] F1 = jac2F1(jac), F2 = jac2F2(jac)
    - [ ] G1 = jac2G1(jac), G2 = jac2G2(jac)
    - [ ] rot = jac2rot(jac)
  - [ ] Save to HDF5: `tests/fixtures/reference_outputs.h5`
    - [ ] Groups: kappa, gamma1, gamma2, F1, F2, G1, G2, rot, jacobian
    - [ ] Attrs: seed, alpha_params

#### Reference Validation
- [ ] Create `tests/fixtures/reference_values.py`
  - [ ] Hard-coded known values for ultra-simple cases (4×4 grids, constant gradient)
  - [ ] For closed-form tests (e.g., uniform field → uniform κ)

#### Tests
- [ ] Create `tests/unit/test_fixtures.py`
  - [ ] `test_synthetic_nfw_map_shape()`: correct dimensions
  - [ ] `test_synthetic_nfw_map_reproducible()`: same seed → identical output
  - [ ] `test_reference_outputs_exist()`: HDF5 file readable, has expected keys
  - [ ] `test_real_corner_patch_loads()`: real FITS loader works

#### Deliverables
- [ ] `tests/fixtures/synth_maps.py` ✓ (3-5 fixtures)
- [ ] `tests/fixtures/reference_outputs.h5` ✓ (Julia baseline)
- [ ] `scripts/generate_julia_reference.jl` ✓
- [ ] `tests/unit/test_fixtures.py` ✓ (5+ tests)

#### Checkpoint
```bash
python -c "from tests.fixtures.synth_maps import simple_constant_map; m = simple_constant_map(); print(m.shape)"
# → (64, 64, 2)

python -c "import h5py; f = h5py.File('tests/fixtures/reference_outputs.h5'); print(list(f.keys()))"
# → ['F1', 'F2', 'G1', 'G2', 'gamma1', 'gamma2', 'jacobian', 'kappa', 'rot']
```

---

### 1d: Jacobian Computation with Numba (Days 12-17)

**Owner**: Person A (heavily performance-critical)  
**Goal**: Core Jacobian computation; must match Julia speed + accuracy

#### Jacobian Module
- [ ] Create `lensing/lensing_models/jacobian.py`

- [ ] Understand Julia's `deriv!()` function
  - [ ] Is it 3-point central difference? 5-point? One-sided boundaries?
  - [ ] Exact stencil: (f[i+1] - f[i-1]) / 2 or similar?
  - [ ] Document findings in comments

- [ ] Implement `finite_difference_1d()` (Numba JIT)
  - [ ] Signature: `(f: array, out: array) -> None` (in-place update)
  - [ ] Stencil matching Julia exactly
  - [ ] Boundary handling (natural or one-sided)
  - [ ] Numba `@jit(nopython=True)` decorator

- [ ] Implement `jacobian_from_deflection(alpha, scale=1.0)`
  - [ ] Input: α field, shape (ny, nx, 2)
  - [ ] Output: jac, shape (ny, nx, 12)
  - [ ] Components [0:4]: first derivatives, [4:8]: second, [8:12]: mixed
  - [ ] Hand-optimized loops (Numba-compiled)
  - [ ] Detailed docstring explaining Jacobian structure
  - [ ] **LaTeX formula**: ∂α_i/∂x_j, ∂²α_i/∂x_j², etc.

- [ ] Implement `_jacobian_kernel()` (Numba-compiled inner loop)
  - [ ] All critical loops here (will be compiled to native code)
  - [ ] Allocate temp arrays `d`, `dd`, `ddd` as in Julia
  - [ ] Iterate: extract 1D slice, compute derivatives, store back

#### Tests
- [ ] Create `tests/unit/test_jacobian.py` (15-20 tests)
  - [ ] `test_jacobian_constant_field()`: α=[c_x, c_y] → jac ≈ 0
  - [ ] `test_jacobian_linear_field()`: α=[ax, cy] → jac[..., 0] ≈ a, etc.
  - [ ] `test_jacobian_shape()`: output shape (64, 64, 12)
  - [ ] `test_jacobian_dtype()`: float32
  - [ ] `test_jacobian_matches_julia_reference(simple_constant_map)`
    - Load reference HDF5, compare all 12 components, rtol=1e-5
  - [ ] `test_jacobian_matches_julia_reference(synthetic_nfw_map)`
  - [ ] `test_jacobian_symmetry()`: If purely radial α, certain derivatives should vanish
  - [ ] `test_jacobian_finite_difference_order()`: ∂²α ≈ (∂α[i+1] - ∂α[i-1]) / 2

#### Performance
- [ ] Create `tests/performance/test_jacobian_perf.py`
  - [ ] `test_jacobian_64x64_under_100ms()`
    - Input: 64×64 map
    - Time: jacobian_from_deflection(map)
    - Assert: elapsed < 0.100 s
  - [ ] Optional: profile with cProfile
    ```python
    import cProfile
    cProfile.run("jacobian_from_deflection(alpha)", sort='cumulative')
    ```
  - [ ] Compare: Julia version timing (if available)

#### Validation vs Julia
- [ ] Run both Julia + Python on `synthetic_nfw_map`
- [ ] Compute relative error: `|jac_py - jac_julia| / |jac_julia|`
- [ ] Document in test; assert <1e-4 (typical floating-point tolerance)

#### Deliverables
- [ ] `lensing/lensing_models/jacobian.py` ✓
- [ ] Numba JIT implementation (speed ≈ Julia)
- [ ] `tests/unit/test_jacobian.py` ✓ (15+ tests)
- [ ] `tests/performance/test_jacobian_perf.py` ✓ (SLA verified)
- [ ] Performance documented: "64×64: XXXms (vs Julia: YYYms)"

#### Checkpoint
```bash
pytest tests/unit/test_jacobian.py -v -k "matches_julia"
# → All pass, error <1e-5

pytest tests/performance/test_jacobian_perf.py -v
# → test_jacobian_64x64_under_100ms PASSED
```

---

### 1e: Observable Extraction (Days 17-20)

**Owner**: Person A  
**Goal**: Extract κ, γ, F, G from Jacobian; add input validation

#### Observables Module
- [ ] Create `lensing/lensing_models/observables.py`

- [ ] Input validation helper: `_validate_jacobian(jac, ...)`
  - [ ] Check dtype == float32
  - [ ] Check ndim == 3
  - [ ] Check jac.shape[2] >= 4 (or exact = expected_size)
  - [ ] Check isfinite (no NaN/Inf)
  - [ ] Helpful error messages

- [ ] Implement observable extraction functions:
  - [ ] `jac2kappa(jac)` → (ny, nx)
    - Formula: κ = (jac[..., 0] + jac[..., 3]) / 2
    - LaTeX + ref: Schneider et al. 1998
    - Input validation
  - [ ] `jac2gamma1(jac)` → (ny, nx)
    - Formula: γ1 = (jac[..., 0] - jac[..., 3]) / 2
    - Orthogonal to κ
  - [ ] `jac2gamma2(jac)` → (ny, nx)
    - Formula: γ2 = jac[..., 2] (or [jac[..., 1] + jac[..., 2]] / 2)?
    - Check Julia code for exact formula
  - [ ] `jac2F1(jac)`, `jac2F2(jac)` → flexion components
    - Formula: F_i = (higher derivatives)
    - See Julia code for exact formula
  - [ ] `jac2G1(jac)`, `jac2G2(jac)` → second flexion
  - [ ] `jac2rot(jac)` → rotation component (if applicable)
  - [ ] `jac2muinv(jac)` → μ_inv (inverse magnification)
    - μ_inv = det(A) = (1-κ)² - (γ1² + γ2²)

#### Tests
- [ ] Create `tests/unit/test_observables.py` (20-30 tests)
  - [ ] Physics bounds:
    - [ ] `test_kappa_range()`: |κ| ≤ 1 (physical bound)
    - [ ] `test_gamma_range()`: |γ1|, |γ2| ≤ 1
    - [ ] `test_mu_inv_positive()`: μ_inv > 0 (no stuck particles)
  - [ ] Regression to Julia:
    - [ ] `test_kappa_matches_julia_reference()` → rtol=1e-5
    - [ ] `test_gamma1_matches_julia_reference()`
    - [ ] `test_gamma2_matches_julia_reference()`
    - [ ] etc. for F1, F2, G1, G2
  - [ ] Input validation:
    - [ ] `test_invalid_dtype_raises()`
    - [ ] `test_invalid_shape_raises()`
    - [ ] `test_nan_input_raises()`
  - [ ] Edge cases:
    - [ ] `test_zero_jacobian_gives_zero_observables()`
    - [ ] `test_identity_jacobian()`: jac = [1,0,0,1,0,0,...] → κ=1, γ=0
  - [ ] Array broadcasting:
    - [ ] Input shape (1, 64, 64, 12) → output (1, 64, 64)
    - [ ] Input shape (10, 10, 12) → output (10, 10)

#### Validation vs Julia
- [ ] For synthetic_nfw_map:
  - Python: jac_py → κ_py → compare to κ_julia from reference
  - Assert all observables match <1e-5

#### Deliverables
- [ ] `lensing/lensing_models/observables.py` ✓
- [ ] All 8+ observable extraction functions
- [ ] `tests/unit/test_observables.py` ✓ (25+ tests)
- [ ] Every function has LaTeX formula + references in docstring
- [ ] Input validation on all functions

#### Checkpoint
```bash
pytest tests/unit/test_observables.py -v
# → All pass (25+ tests)

pytest tests/unit/test_observables.py -v -k "matches_julia"
# → All regression tests pass
```

---

### **Phase 1 Complete** ✓

**Summary Checkpoint**:
```bash
pytest tests/unit/ -v
# → PASSED (60+ tests)
# Breakdown:
#   - test_cosmology.py: 15 tests
#   - test_config.py: 5 tests
#   - test_nfw_profile.py: 25 tests
#   - test_jacobian.py: 15 tests
#   - test_observables.py: 25 tests

pytest tests/performance/ -v
# → PASSED (SLAs met)

coverage run -m pytest tests/unit
coverage report
# → Coverage >85%

cd tests/fixtures && ls -la reference_outputs.h5 synth_maps.py
# → Both exist, HDF5 readable
```

**Definition of "Phase 1 Done"**:
- All 16 core functions implemented (cosmology, NFW, Jacobian, 8 observables)
- 60+ unit tests passing
- Julia reference outputs stored
- Synthetic fixtures ready for integration testing
- Coverage >85%

**Time elapsed**: 2-3 weeks (Days 1-20)

---

## Phase 2: Data I/O & Pipeline Integration (Weeks 3-4)

*[Detailed checklist for Phase 2]*

### 2a: I/O Layer (Days 22-27)

- [ ] Create `lensing/io/fortran_binary.py`
  - [ ] `read_fortran_deflection_map(filepath)`
  - [ ] Tests: correct shape, dtype, values sample

- [ ] Create `lensing/io/hdf5_backend.py`
  - [ ] `HDF5DataStore` class
  - [ ] `write_deflection_map()`, `read_deflection_map()`
  - [ ] Compression validation (>80% size reduction)
  - [ ] Tests: roundtrip I/O, metadata preservation

- [ ] Create `scripts/convert_fortran_to_hdf5.py`
  - [ ] One-time migration: `Data/lensing_maps/*.bin` → `data/processed/deflection_maps.h5`
  - [ ] Verify file sizes before/after
  - [ ] Test roundtrip read: original ≈ converted

- [ ] Tests: `tests/unit/test_hdf5_backend.py` (10+ tests)

### 2b: Data Loader (Days 27-30)

- [ ] Create `lensing/pipeline/data_loader.py`
  - [ ] `ObservationData` dataclass
  - [ ] `DataPipeline` class
  - [ ] `load_observation(z_id)`
  - [ ] `_load_galaxy_catalog()`
  - [ ] Filter by z, type, mass

- [ ] Tests: `tests/integration/test_data_pipeline.py` (5+ tests)

### 2c: Correlator (Days 30-36)

- [ ] Create `lensing/pipeline/correlator.py`
  - [ ] `Correlator` class
  - [ ] `compute_kk_correlation()` via treecorr
  - [ ] `_extract_at_positions()` (interpolation)

- [ ] Tests: `tests/unit/test_correlator.py` (5+ tests)

### 2d: Fitter (Days 36-40)

- [ ] Create `lensing/pipeline/fitter.py`
  - [ ] `NFWFitter` class
  - [ ] `fit_nfw()` using scipy.optimize.curve_fit

- [ ] Tests: `tests/unit/test_fitter.py` (5+ tests)

---

## Phase 3: Orchestration (Week 5)

### 3a: Main Entry Point

- [ ] Create `scripts/run_pipeline.py`
  - [ ] CLI: config + z-ids + output dir
  - [ ] Full pipeline orchestration

- [ ] Create `lensing/pipeline/__init__.py` (Pipeline class)

- [ ] Create `data/config/thesis.yaml`

### 3b: Integration Tests

- [ ] Create `tests/integration/test_pipeline_e2e.py`
  - [ ] Single z-slice run
  - [ ] Multi z-slice run
  - [ ] Real data validation

**Checkpoint**:
```bash
python scripts/run_pipeline.py data/config/thesis.yaml --z-ids 50
# → Success, output/ created with jacobians/, correlations/, fits/
```

---

## Phase 4: Plotting (Week 6)

- [ ] `lensing/plotting/correlation_plots.py`
- [ ] `lensing/plotting/map_plots.py`
- [ ] `scripts/plot_correlations.py`
- [ ] Example figures generated

---

## Phase 5: Documentation (Week 7)

- [ ] `docs/physics.md`: formulas + citations
- [ ] `docs/api.md`: function reference
- [ ] `docs/architecture.md`: design overview
- [ ] `docs/quickstart.md`: how to run
- [ ] `README.md`: top-level overview
- [ ] All docstrings finalized

---

## Phase 6: Validation & Polish (Week 8)

- [ ] Julia regression tests (`tests/integration/test_real_data.py`)
- [ ] Performance profiling
- [ ] Final code review
- [ ] CI/CD green across all Python versions (3.9, 3.10, 3.11)

**Final Checkpoint**:
```bash
pytest tests/ -v --tb=short
# → ALL PASS (100+ tests)

coverage report
# → >85% coverage

mypy lensing/
# → Success

python scripts/run_pipeline.py data/config/thesis.yaml --z-ids 50 100 150
# → <30 seconds, produces full results

python scripts/validate_against_julia.py tests/fixtures/reference_outputs.h5
# → All observables match <0.1% error
```

---

## Final Checklist

### Before Thesis Submission
- [ ] All code committed to git
- [ ] All tests passing (`pytest tests/ -v`)
- [ ] mypy clean (`mypy lensing/`)
- [ ] Coverage >85%
- [ ] README + quickstart accessible
- [ ] Figures reproducible from scripts
- [ ] All formulas cited in docstrings

### For Open-Source Release (Optional)
- [ ] LICENSE file (MIT/GPL)
- [ ] CONTRIBUTING.md (dev guidelines)
- [ ] PyPI setup (optional: `pip install lensing-pipeline`)
- [ ] Zenodo DOI (if desired)

---

**Last Updated**: February 17, 2026  
**Status**: Ready for Phase 0 kickoff

