# Weak-Lensing Pipeline Modernization Roadmap

**Project**: Convert & refactor Laurant's Julia weak-lensing ray-tracing code to modern, testable Python  
**Timeline**: 8 weeks (can be accelerated with parallelization)  
**Start date**: February 2026  
**Status**: Planning phase  

---

## Executive Summary

This roadmap modernizes a scattered Julia codebase (raytrace + correlations + fitting) into a clean, testable Python pipeline. Key improvements:

| Dimension | Current | Target |
|-----------|---------|--------|
| **Testability** | None | Comprehensive pytest suite + synthetic fixtures |
| **Code reuse** | ~60% duplication | ~100% DRY (single NFW/cosmology module) |
| **Reproducibility** | Hardcoded paths/params | Config-driven + documented formulas |
| **Dependencies** | Julia + custom C libs | Standard Python stack (numpy, scipy, astropy, treecorr) |
| **Architecture** | Monolithic scripts | Modular pipeline (I/O → Physics → Plotting) |
| **Iteration speed** | ~5 min/test (full data) | ~100ms/test (synthetic data) |
| **Documentation** | Scattered comments | Full LaTeX docstrings + paper citations |

---

## Directory Structure (Target)

```
lensing-pipeline/
├── pyproject.toml                 # Poetry/setuptools config
├── README.md                      # High-level overview
├── MODERNIZATION_ROADMAP.md       # This file
│
├── lensing/                       # Main package
│   ├── __init__.py
│   ├── config.py                  # Config management (YAML loader)
│   ├── cosmology.py               # Cosmology module (astropy wrapper)
│   ├── lensing_models/            # Lensing theory
│   │   ├── __init__.py
│   │   ├── nfw_profile.py         # NFW halo (all functions consolidated)
│   │   ├── jacobian.py            # α → Jacobian computation (Numba JIT)
│   │   └── observables.py         # Jacobian → κ, γ, F, G (extraction)
│   ├── pipeline/                  # Orchestration
│   │   ├── __init__.py
│   │   ├── data_loader.py         # Read FITS/HDF5/Fortran binary
│   │   ├── correlator.py          # Compute 2-point correlations
│   │   └── fitter.py              # NFW parameter fitting
│   ├── io/                        # I/O utilities
│   │   ├── __init__.py
│   │   ├── fortran_binary.py      # Read Fortran deflection maps
│   │   ├── hdf5_backend.py        # HDF5 read/write
│   │   └── converters.py          # Format converters (bin→HDF5, etc.)
│   └── plotting/                  # Visualization layer (decoupled)
│       ├── __init__.py
│       ├── correlation_plots.py
│       └── map_plots.py
│
├── tests/                         # Comprehensive test suite
│   ├── conftest.py                # Pytest fixtures
│   ├── fixtures/                  # Test data
│   │   ├── synth_maps.py          # Synthetic map generators
│   │   ├── reference_outputs.h5   # Julia baseline (for regression)
│   │   └── real_deflection_z0050_corner.fits
│   ├── unit/                      # Unit tests (no I/O)
│   │   ├── test_cosmology.py      # Cosmology functions
│   │   ├── test_nfw_profile.py    # NFW functions (closed-form tests)
│   │   ├── test_jacobian.py       # Jacobian computation
│   │   └── test_observables.py    # Observable extraction
│   ├── integration/               # End-to-end tests
│   │   ├── test_pipeline_e2e.py   # Full pipeline on synthetic data
│   │   └── test_real_data.py      # Validation against Julia outputs
│   └── performance/               # Timing benchmarks
│       └── test_timings.py        # SLA checks
│
├── scripts/                       # Entry points for users
│   ├── run_pipeline.py            # Main: convert config → results
│   ├── convert_fortran_to_hdf5.py # One-time Fortran→HDF5 migration
│   ├── validate_against_julia.py  # Regression test (compare to Julia)
│   └── plot_correlations.py       # Plotting convenience
│
├── data/                          # Input data + outputs
│   ├── config/
│   │   ├── default.yaml           # Default cosmology/parameters
│   │   └── thesis.yaml            # Thesis-specific config
│   ├── raw/                       # Original Fortran binaries, FITS
│   │   ├── deflection_maps/
│   │   ├── galaxy_catalogs/
│   │   └── correlation_data/
│   └── processed/                 # Converted to HDF5, outputs
│       ├── deflection_maps.h5
│       ├── correlation_functions.h5
│       └── fit_results.csv
│
├── docs/                          # Documentation
│   ├── physics.md                 # Physics formulas & citations
│   ├── architecture.md            # System design
│   ├── api.md                     # Function reference
│   └── quickstart.md              # Getting started
│
└── requirements.txt               # Dependencies (pip)
```

---

## Phase Breakdown (8 Weeks)

### **Phase 0: Preparation (Days 1-2)**

**Goal**: Set up repository, environment, testing infrastructure.

**Tasks**:
- [ ] Create Python package structure (pyproject.toml, setup.py)
- [ ] Initialize pytest + conftest.py
- [ ] Set up pre-commit hooks (black, isort, mypy)
- [ ] Create GitHub Actions CI/CD pipeline
- [ ] Set up documentation skeleton (Sphinx or MkDocs)

**Deliverables**:
- Empty package structure with passing CI
- All developers can `pip install -e .` and run `pytest`

**Time**: 1-2 days  
**Team**: 1 person (setup)

---

### **Phase 1: Core Physics + Testing Foundations (Weeks 1-2)**

**Goal**: Build testable core modules; establish regression baseline.

#### **1a: Cosmology & Configuration (Days 3-5)**

**Tasks**:
1. Create `lensing/config.py`:
   - YAML loader for configuration
   - Centralized cosmology parameters
   - Type-validated dataclasses

2. Create `lensing/cosmology.py`:
   - Wrap astropy.cosmology
   - Implement: D_A(z), D_C(z), Σ_crit, E(z)
   - Add unit tests for known values (e.g., PLANCK 2018)

**Code outline**:
```python
# config.py
@dataclass
class CosmologyConfig:
    h: float = 0.7
    omega_m: float = 0.3
    omega_lambda: float = 0.7
    
    @property
    def astropy_cosmology(self):
        return FlatLambdaCDM(H0=self.h*100*u.km/u.s/u.Mpc, Om0=self.omega_m)

# cosmology.py
class CosmologyCalculator:
    def __init__(self, config: CosmologyConfig):
        self.cosmo = config.astropy_cosmology
    
    def angular_diameter_distance(self, z1: float, z2: float) -> u.Quantity:
        """D_A(z1, z2) with validation."""
        if not (0 <= z1 < z2 < 100):
            raise ValueError(f"Invalid redshifts: z1={z1}, z2={z2}")
        return self.cosmo.angular_diameter_distance_z1z2(z1*u.dimensionless_unscaled, 
                                                          z2*u.dimensionless_unscaled)
    # ...
```

**Tests** (5-10):
- `test_dda_matches_astropy()`: Verify against astropy baseline
- `test_dda_boundary_cases()`: z1=z2, z1=0, z1>z2 errors
- `test_comoving_distance()`: Known values from Hogg 1999

**Deliverables**:
- `lensing/config.py` ✓
- `lensing/cosmology.py` ✓
- 10-15 passing tests
- `tests/unit/test_cosmology.py` ✓

**Time**: 2-3 days  
**Owner**: 1 person

---

#### **1b: NFW Profile (Consolidated, Days 5-9)**

**Goal**: Single source of truth for all NFW functions (consolidating from nfw_tk.jl, correlation_function.jl, correlation_plots.jl).

**Tasks**:
1. Create `lensing/lensing_models/nfw_profile.py`:
   - `NFWProfile` class with static methods
   - Functions: `delta_c(c)`, `ksrs2mc()`, `mc2ksrs()`, `c_m_K16()`, `Flens()`, `Fprime()`, `Hlens()`, `Glens()`
   - Add full LaTeX docstrings + paper citations

2. Add formula validation docstrings:
```python
@staticmethod
def delta_c(c: float) -> float:
    r"""
    Critical density contrast for NFW profile.
    
    .. math::
        \delta_c(c) = \frac{c^3}{\ln(1+c) - c/(1+c)} / 3
    
    References:
        Navarro, Frenk & White (1997), Eq. 1
        Bartelmann & Schneider (2001), Eq. 42
    """
    return c**3 / (np.log(1+c) - c/(1+c)) / 3
```

3. Create `tests/unit/test_nfw_profile.py`:
   - Closed-form tests (e.g., delta_c(c=4) = exact value from paper)
   - Round-trip tests: `mc2ksrs() → ksrs2mc()` should recover input
   - Edge case validation: mass/c ranges

**Deliverables**:
- `lensing/lensing_models/nfw_profile.py` ✓
- All functions from nfw_tk.jl + correlation_*.jl consolidated
- 15-20 unit tests
- LaTeX docstrings for every function

**Time**: 4-5 days  
**Owner**: 1 person

---

#### **1c: Synthetic Test Data Generators (Days 9-12)**

**Goal**: Enable fast iteration without full-size files.

**Tasks**:
1. Create `tests/fixtures/synth_maps.py`:

```python
@pytest.fixture
def simple_constant_map() -> np.ndarray:
    """Trivial: uniform deflection field α = [x, y]."""
    size = 64
    x = np.linspace(0, 1, size)
    y = np.linspace(0, 1, size)
    xx, yy = np.meshgrid(x, y)
    alpha = np.stack([xx, yy], axis=2).astype(np.float32)
    return alpha

@pytest.fixture  
def synthetic_nfw_map(seed=42) -> Tuple[np.ndarray, float]:
    """64×64 deflection map from single NFW halo at center.
    
    Returns:
        (alpha, z_source): deflection field + source redshift
    """
    rng = np.random.default_rng(seed)
    size = 64
    center = size / 2
    
    # Create coordinate grid in arcseconds (1 pixel = 1 arcsec)
    x = np.arange(size) - center
    y = np.arange(size) - center
    xx, yy = np.meshgrid(x, y)
    r = np.sqrt(xx**2 + yy**2) + 0.1  # Avoid singularity
    
    # NFW deflection: α_r = 4π(rs × ks/r) × F(r/rs)
    rs = 10.0  # scale radius in arcsec
    ks = 0.1   # shear normalization
    
    # Deflection magnitude from NFW profile
    alpha_r = compute_nfw_deflection(r, rs, ks)
    
    # Decompose to components (radial symmetry):
    # α_x = (α_r / r) * x, α_y = (α_r / r) * y
    alpha_x = (alpha_r / r) * xx
    alpha_y = (alpha_r / r) * yy
    
    alpha = np.stack([alpha_x, alpha_y, alpha_x, alpha_y], axis=2).astype(np.float32)
    return alpha, 1.0
```

2. Create `tests/fixtures/reference_outputs.h5`:
   - Run Laurant's Julia code on `synthetic_nfw_map()`
   - Save outputs: jac, kappa, gamma1, gamma2, F1, F2, G1, G2
   - This becomes regression baseline

**Script**: `scripts/generate_julia_reference.jl`
```julia
using HDF5
include("raytrace_tk.jl")

# Load synthetic map
alpha = load_synthetic_nfw_map()

# Compute all observables
jac = alpha2jac(alpha)
kappa = jac2kappa(jac)
gamma1 = jac2gamma1(jac)
# ... etc

# Save as HDF5
h5write("tests/fixtures/reference_outputs.h5", "kappa", kappa)
# ... etc
```

**Deliverables**:
- `tests/fixtures/synth_maps.py` ✓ (3-5 fixtures)
- `tests/fixtures/reference_outputs.h5` ✓ (Julia baseline)
- `scripts/generate_julia_reference.jl` ✓

**Time**: 3-4 days  
**Owner**: 1 person (requires Julia environment)

---

#### **1d: Jacobian Computation with Numba (Days 12-17)**

**Goal**: Core performance-critical function; must match Julia speed + accuracy.

**Tasks**:
1. Examine Julia's `deriv!()` function → understand finite-difference stencil

2. Create `lensing/lensing_models/jacobian.py`:

```python
from numba import jit
import numpy as np

@jit(nopython=True)
def finite_difference_1d(f: np.ndarray, out: np.ndarray) -> None:
    """In-place 1st derivative via 3-point stencil.
    
    f'[i] ≈ (f[i+1] - f[i-1]) / 2
    Boundary handled by one-sided stencils.
    """
    n = len(f)
    # Central differences
    for i in range(1, n-1):
        out[i] = (f[i+1] - f[i-1]) / 2.0
    # Boundaries
    out[0] = f[1] - f[0]
    out[n-1] = f[n-1] - f[n-2]

def jacobian_from_deflection(alpha: np.ndarray, 
                              scale: float = 1.0) -> np.ndarray:
    """Convert deflection angles → full Jacobian matrix.
    
    Args:
        alpha: Deflection field, shape (ny, nx, 2) with components [α_x, α_y]
        scale: Pixel scale (grid spacing)
    
    Returns:
        jac: Jacobian array, shape (ny, nx, 12) with:
             [0:4] = 1st derivatives (∂α_i/∂x_j)
             [4:8] = 2nd derivatives (∂²α_i/∂x_j²)
             [8:12] = mixed derivatives (∂²α_i/∂x_j∂x_k)
    
    References:
        Thesis Section 3.1: Jacobian computation via finite differences
    """
    ny, nx = alpha.shape[:2]
    jac = np.zeros((ny, nx, 12), dtype=np.float32)
    
    # Implementation: hand-optimized loops (Numba JIT will compile to machine code)
    _jacobian_kernel(alpha, scale, jac)
    
    return jac

@jit(nopython=True)
def _jacobian_kernel(alpha, scale, jac):
    """Actual computation loop (Numba-compiled)."""
    ny, nx = alpha.shape[:2]
    conv_x = nx / scale
    conv_y = ny / scale
    
    # Allocate temp arrays
    d = np.zeros(nx, dtype=np.float32)
    dd = np.zeros(nx, dtype=np.float32)
    ddd = np.zeros(nx, dtype=np.float32)
    
    # 1st derivatives: ∂α_x/∂x (index 0), ∂α_y/∂x (index 1)
    for j in range(ny):
        for i in range(nx):
            d[i] = alpha[i, j, 0] * conv_x
        finite_difference_1d(d, dd)  # 1st deriv
        finite_difference_1d(dd, ddd)  # 2nd deriv
        for i in range(nx):
            jac[i, j, 0] = dd[i]  # ∂α_x/∂x
            jac[i, j, 4] = ddd[i]  # ∂²α_x/∂x²
    # ... (similar for other components)
```

3. Validation against Julia:

```python
# tests/unit/test_jacobian.py
def test_jacobian_matches_julia_reference(synthetic_nfw_map):
    """Regression: ensure Jacobian matches Julia baseline."""
    alpha, _ = synthetic_nfw_map
    
    jac_py = jacobian_from_deflection(alpha, scale=1.0)
    
    with h5py.File('tests/fixtures/reference_outputs.h5', 'r') as f:
        jac_ref = f['jacobian'][:]
    
    # Typical tolerance: 1e-5 (floating point precision)
    assert np.allclose(jac_py, jac_ref, rtol=1e-5, atol=1e-10)

def test_jacobian_constant_field():
    """Edge case: zero deflection → zero Jacobian."""
    alpha = np.zeros((64, 64, 2), dtype=np.float32)
    jac = jacobian_from_deflection(alpha)
    assert np.allclose(jac, 0)
```

4. Benchmark:
```python
# tests/performance/test_timings.py
@pytest.mark.benchmark
def test_jacobian_perf(synthetic_nfw_map):
    """64×64 Jacobian computation should be <100ms."""
    alpha, _ = synthetic_nfw_map
    
    start = timer()
    jac = jacobian_from_deflection(alpha, scale=1.0)
    elapsed = timer() - start
    
    assert elapsed < 0.100, f"Too slow: {elapsed:.3f}s (target: <100ms)"
```

**Deliverables**:
- `lensing/lensing_models/jacobian.py` ✓
- Numba JIT implementation (speed ~ Julia)
- `tests/unit/test_jacobian.py` ✓ (10+ tests)
- Performance baseline documented

**Time**: 5-6 days  
**Owner**: 1 person (requires performance profiling)

**Notes**:
- Profile against Julia version to confirm speed
- If Numba disappoints, fall back to scipy.ndimage.convolve1d (still >100× faster than naive Python)

---

#### **1e: Observable Extraction (Days 17-20)**

**Goal**: Extract κ, γ, F, G from Jacobian.

**Tasks**:
1. Create `lensing/lensing_models/observables.py`:

```python
def jac2kappa(jac: np.ndarray) -> np.ndarray:
    r"""
    Convergence from Jacobian trace.
    
    .. math::
        \kappa = \frac{1}{2}(\partial_1 \alpha_1 + \partial_2 \alpha_2)
             = \frac{1}{2} \text{Tr}(A)
    
    where A is the magnification matrix (Jacobian).
    
    References:
        Schneider et al. (1998), Eq. 2.7
    """
    if jac.shape[2] < 4:
        raise ValueError(f"Jacobian must have ≥4 components, got {jac.shape[2]}")
    
    # κ = (jac[..., 0] + jac[..., 3]) / 2
    kappa = (jac[..., 0] + jac[..., 3]) / 2.0
    
    # Validation
    if not np.isfinite(kappa).all():
        n_bad = (~np.isfinite(kappa)).sum()
        raise ValueError(f"Output contains {n_bad} non-finite values")
    
    return kappa

# Similar for jac2gamma1, jac2gamma2, jac2F1, jac2F2, jac2G1, jac2G2
```

2. Add input validation contracts (Issue 2-A from Code Quality):

```python
def _validate_jacobian(jac: np.ndarray, 
                       expected_shape: Tuple[int, int, int] = None,
                       allow_nan: bool = False) -> None:
    """Helper: validate Jacobian structure."""
    if jac.dtype != np.float32:
        raise TypeError(f"Expected float32, got {jac.dtype}")
    if len(jac.shape) != 3:
        raise ValueError(f"Expected 3D array, got {len(jac.shape)}D")
    if jac.shape[2] < 4:
        raise ValueError(f"jac[..., :] must have ≥4 components, got {jac.shape[2]}")
    if not allow_nan and not np.isfinite(jac).all():
        raise ValueError(f"Input contains {np.isnan(jac).sum()} NaNs, {np.isinf(jac).sum()} Infs")
```

3. Tests:

```python
# tests/unit/test_observables.py
def test_kappa_range(synthetic_nfw_map):
    """Physical κ must satisfy |κ| ≤ 1."""
    alpha, _ = synthetic_nfw_map
    jac = jacobian_from_deflection(alpha)
    kappa = jac2kappa(jac)
    assert np.all(np.abs(kappa) <= 1.0)

def test_kappa_matches_julia(synthetic_nfw_map):
    """Regression: κ matches Julia baseline."""
    alpha, _ = synthetic_nfw_map
    jac = jacobian_from_deflection(alpha)
    kappa = jac2kappa(jac)
    
    with h5py.File('tests/fixtures/reference_outputs.h5', 'r') as f:
        kappa_ref = f['kappa'][:]
    
    assert np.allclose(kappa, kappa_ref, rtol=1e-5)

def test_kappa_invalid_input():
    """Invalid input raises helpful error."""
    bad_jac = np.ones((64, 64, 2))  # Too few components
    with pytest.raises(ValueError, match="must have ≥4 components"):
        jac2kappa(bad_jac)
```

**Deliverables**:
- `lensing/lensing_models/observables.py` ✓
- 20-30 tests (one per observable, plus edge cases)
- Input validation on all functions
- All functions matched against Julia reference

**Time**: 3-4 days  
**Owner**: 1 person

---

#### **Phase 1 Summary**

| Task | Days | Owner | Status |
|------|------|-------|--------|
| 1a. Cosmology | 2-3 | Person A | ▢ |
| 1b. NFW Profile | 4-5 | Person A | ▢ |
| 1c. Synthetic test data | 3-4 | Person B | ▢ |
| 1d. Jacobian (Numba) | 5-6 | Person A | ▢ |
| 1e. Observables | 3-4 | Person A | ▢ |
| **Total** | **17-21 days** | | |

**Milestones**:
- ✓ All core physics functions implemented + regression-tested
- ✓ 60+ unit tests passing
- ✓ Synthetic fixtures ready for integration testing
- ✓ Performance baseline documented (Jacobian <100ms for 64×64)

**Definition of done for Phase 1**:
```bash
pytest tests/unit/ -v  # All pass
pytest tests/performance/ -v  # All SLAs met
coverage report --min-coverage=85  # >85% code coverage
```

---

### **Phase 2: Data I/O & Pipeline Integration (Weeks 3-4)**

**Goal**: Build data loading, HDF5 conversion, correlation computation.

#### **2a: I/O Layer (Days 22-27)**

**Tasks**:
1. Create `lensing/io/fortran_binary.py`:
   - Read Fortran binary deflection maps
   - Validate output shape/type

2. Create `lensing/io/hdf5_backend.py`:
   - High-level read/write interface
   - Hierarchical structure (z-slices, observables)

3. Create `scripts/convert_fortran_to_hdf5.py`:
   - One-time migration script
   - Convert all Fortran bins → single HDF5 file
   - Compress (reduce 1.7GB → ~200MB)

**Code outline**:
```python
# lensing/io/hdf5_backend.py
class HDF5DataStore:
    def __init__(self, filepath: str):
        self.path = filepath
    
    def write_deflection_map(self, z_id: int, alpha: np.ndarray, redshift: float):
        """Store deflection field for snapshot z_id."""
        with h5py.File(self.path, 'a') as f:
            grp = f.create_group(f'deflection_maps/{z_id:04d}')
            grp.create_dataset('alpha', data=alpha, compression='gzip', compression_opts=4)
            grp.attrs['redshift'] = redshift
    
    def read_deflection_map(self, z_id: int) -> Tuple[np.ndarray, float]:
        """Load deflection field."""
        with h5py.File(self.path, 'r') as f:
            grp = f[f'deflection_maps/{z_id:04d}']
            alpha = grp['alpha'][:]
            z = grp.attrs['redshift']
        return alpha, z
```

4. Tests:
```python
# tests/unit/test_hdf5_backend.py
def test_write_read_roundtrip(tmp_path):
    """Write → read should recover data exactly."""
    store = HDF5DataStore(tmp_path / 'test.h5')
    alpha = np.random.randn(64, 64, 2).astype(np.float32)
    
    store.write_deflection_map(50, alpha, 0.21)
    alpha_read, z = store.read_deflection_map(50)
    
    assert np.allclose(alpha, alpha_read)
    assert np.isclose(z, 0.21)

def test_compression_ratio(tmp_path):
    """Compression should reduce size by ~80%."""
    # Create large array
    alpha = np.random.randn(1024, 1024, 2).astype(np.float32)
    
    store = HDF5DataStore(tmp_path / 'test.h5')
    store.write_deflection_map(50, alpha, 0.21)
    
    file_size = (tmp_path / 'test.h5').stat().st_size
    uncompressed_size = alpha.nbytes
    
    ratio = file_size / uncompressed_size
    assert ratio < 0.2, f"Poor compression: {ratio:.1%} (target: <20%)"
```

**Deliverables**:
- `lensing/io/fortran_binary.py` ✓
- `lensing/io/hdf5_backend.py` ✓
- `scripts/convert_fortran_to_hdf5.py` ✓
- Full data migration (Fortran → HDF5) ✓
- 10-15 I/O tests

**Time**: 4-6 days  
**Owner**: 1 person

---

#### **2b: Data Loader (Days 27-30)**

**Goal**: High-level interface to load catalogs + maps.

**Tasks**:
1. Create `lensing/pipeline/data_loader.py`:

```python
@dataclass
class ObservationData:
    """Container for single snapshot."""
    z_id: int
    redshift: float
    deflection_map: np.ndarray  # (ny, nx, 2)
    galaxy_positions: np.ndarray  # (n_gal, 2) in arcsec
    galaxy_masses: np.ndarray  # (n_gal,) in Msun
    galaxy_types: np.ndarray  # (n_gal,) 0=central, 1=satellite

class DataPipeline:
    def __init__(self, config: CosmologyConfig, data_store: HDF5DataStore):
        self.config = config
        self.store = data_store
    
    def load_observation(self, z_id: int) -> ObservationData:
        """Load all data for snapshot z_id."""
        alpha, z = self.store.read_deflection_map(z_id)
        galaxies = self._load_galaxy_catalog(z_id, z)
        return ObservationData(z_id, z, alpha, **galaxies)
    
    def _load_galaxy_catalog(self, z_id: int, z: float) -> Dict:
        """Load + filter galaxy catalog for this snapshot."""
        # Read from FITS (use astropy)
        # Filter by redshift, type, mass
        return {...}
```

2. Tests:
```python
# tests/integration/test_data_pipeline.py
def test_load_observation_complete():
    """Load observation has all required fields."""
    pipeline = DataPipeline(cosmology_config, store)
    obs = pipeline.load_observation(50)
    
    assert obs.z_id == 50
    assert 0.20 < obs.redshift < 0.23  # z0050 ~ 0.21
    assert obs.deflection_map.shape[:2] == (1024, 1024)  # or whatever size
    assert len(obs.galaxy_positions) > 100
```

**Deliverables**:
- `lensing/pipeline/data_loader.py` ✓
- 5-10 integration tests

**Time**: 2-3 days  
**Owner**: 1 person

---

#### **2c: Correlator (Days 30-36)**

**Goal**: Compute 2-point correlations using treecorr.

**Tasks**:
1. Create `lensing/pipeline/correlator.py`:

```python
import treecorr

class Correlator:
    """Compute 2-point correlation functions."""
    
    def __init__(self, config: CosmologyConfig):
        self.config = config
        self.rmin = 0.1  # arcsec
        self.rmax = 100  # arcsec
        self.nbins = 20
    
    def compute_kk_correlation(self, obs: ObservationData) -> CorrelationResult:
        """κ × κ correlation at galaxy positions."""
        # Extract convergence at galaxy locations
        # (Requires map interpolation or nearest-neighbor lookup)
        kappa_at_gal = self._extract_at_positions(obs.convergence_map, 
                                                   obs.galaxy_positions)
        
        # Build treecorr catalog
        cat = treecorr.Catalog(x=obs.galaxy_positions[:, 0],
                               y=obs.galaxy_positions[:, 1],
                               k=kappa_at_gal)
        
        # Compute correlation
        kk = treecorr.KKCorrelation(self.rmin, self.rmax, self.nbins)
        kk.process(cat)
        
        return CorrelationResult(
            r_bins=kk.meanr,
            xi=kk.xi,
            error=kk.sigma,
            estimator='natural'
        )
    
    def _extract_at_positions(self, field_map: np.ndarray, 
                               positions: np.ndarray) -> np.ndarray:
        """Extract field values at galaxy positions."""
        # Bilinear interpolation or nearest-neighbor
        from scipy.interpolate import RegularGridInterpolator
        # ...
```

2. Tests:
```python
# tests/unit/test_correlator.py
def test_correlation_symmetric(synthetic_nfw_map):
    """κ×κ should be symmetric: ξ(r) ≈ ξ(-r)."""
    # ...

def test_correlation_positive_mass_positive():
    """Massive object → positive κ×κ correlation."""
    # ...
```

**Deliverables**:
- `lensing/pipeline/correlator.py` ✓
- 5-10 tests
- Documentation of Landy-Szalay estimator

**Time**: 5-6 days  
**Owner**: 1 person (may overlap with 2a/2b)

---

#### **2d: Fitter (Days 36-40)**

**Goal**: Fit NFW parameters to observed correlations.

**Tasks**:
1. Create `lensing/pipeline/fitter.py`:

```python
class NFWFitter:
    """Fit NFW model to observed correlation functions."""
    
    def __init__(self, cosmology: CosmologyCalculator):
        self.cosmo = cosmology
    
    def fit_nfw(self, 
                obs: CorrelationResult,
                initial_guess: Tuple[float, float] = (1e14, 4.0)) -> FitResult:
        """
        Fit M200 and concentration c to observed κ×κ correlation.
        
        Uses scipy.optimize.minimize with Levenberg-Marquardt.
        """
        def model(r, m200, c):
            """Predict κ×κ from NFW parameters."""
            ks, rs = self.nfw_model.mc2ksrs(m200, c, 
                                             q=1.0, z=obs.redshift)
            kappa = self.nfw_model.kappafunc(r, ks, rs)
            return kappa
        
        popt, pcov = curve_fit(model, obs.r_bins, obs.xi, 
                                p0=initial_guess,
                                sigma=obs.error)
        
        return FitResult(m200=popt[0], c=popt[1], cov=pcov)
```

2. Tests:
```python
# tests/unit/test_fitter.py
def test_fitter_recovers_injected_params():
    """Fit synthetic correlation generated from known NFW → recover params."""
    m_true, c_true = 1e14, 4.0
    
    # Generate synthetic correlation from NFW
    r = np.logspace(-1, 2, 20)
    kk_synthetic = nfw_model.kappafunc(r, *mc2ksrs(m_true, c_true, ...))
    
    # Fit
    result = fitter.fit_nfw(CorrelationResult(r, kk_synthetic + noise, ...))
    
    # Should recover (within fitting error)
    assert np.isclose(result.m200, m_true, rtol=0.1)
    assert np.isclose(result.c, c_true, rtol=0.1)
```

**Deliverables**:
- `lensing/pipeline/fitter.py` ✓
- 5-10 tests
- Uncertainty quantification (covariance propagation)

**Time**: 3-4 days  
**Owner**: 1 person

---

#### **Phase 2 Summary**

| Task | Days | Owner | Status |
|------|------|-------|--------|
| 2a. I/O layer | 4-6 | Person B | ▢ |
| 2b. Data loader | 2-3 | Person B | ▢ |
| 2c. Correlator | 5-6 | Person C | ▢ |
| 2d. Fitter | 3-4 | Person A | ▢ |
| **Total** | **14-19 days** | | |

**Milestones**:
- ✓ All data loads from HDF5 correctly
- ✓ Correlations computed (validated vs treecorr examples)
- ✓ NFW fitting works end-to-end
- ✓ Integration tests passing

---

### **Phase 3: Full Pipeline Orchestration (Week 5)**

**Goal**: Wire everything together; run full pipeline end-to-end on small data.

#### **3a: Main Entry Point (Days 41-44)**

**Tasks**:
1. Create `scripts/run_pipeline.py`:

```python
#!/usr/bin/env python3
"""Main entry point: config → results."""

import argparse
from lensing.config import load_config
from lensing.pipeline import Pipeline

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('config', help='Path to config.yaml')
    parser.add_argument('--output', help='Output directory', default='./output')
    parser.add_argument('--z-ids', nargs='+', type=int, 
                        help='Redshift snapshots to process (default: all)')
    args = parser.parse_args()
    
    config = load_config(args.config)
    pipeline = Pipeline(config)
    
    results = pipeline.run(z_ids=args.z_ids, output_dir=args.output)
    
    print(f"✓ Pipeline complete. Results in {args.output}/")
    print(f"  - Jacobians: {len(results.jacobians)} z-slices")
    print(f"  - Correlations: {results.n_correlations}")
    print(f"  - Fits: {results.n_fits}")

if __name__ == '__main__':
    main()
```

2. Create `lensing/pipeline/__init__.py` (Pipeline orchestrator):

```python
class Pipeline:
    def __init__(self, config: CosmologyConfig):
        self.config = config
        self.data_loader = DataPipeline(config)
        self.correlator = Correlator(config)
        self.fitter = NFWFitter(CosmologyCalculator(config))
    
    def run(self, z_ids: List[int], output_dir: str) -> PipelineResult:
        results = PipelineResult()
        
        for z_id in z_ids:
            print(f"Processing z_id={z_id}...")
            
            # Load
            obs = self.data_loader.load_observation(z_id)
            
            # Compute Jacobian + observables
            jac = jacobian_from_deflection(obs.deflection_map)
            obs.jacobian = jac
            obs.convergence = jac2kappa(jac)
            obs.gamma1 = jac2gamma1(jac)
            obs.gamma2 = jac2gamma2(jac)
            
            results.jacobians[z_id] = jac
            
            # Correlate
            corr = self.correlator.compute_kk_correlation(obs)
            results.correlations[z_id] = corr
            
            # Fit
            fit = self.fitter.fit_nfw(corr)
            results.fits[z_id] = fit
            
            # Save intermediate
            self.data_loader.store.write_jacobian(z_id, jac)
            self.data_loader.store.write_correlation(z_id, corr)
        
        return results
```

3. Config file (YAML):

```yaml
# data/config/thesis.yaml
cosmology:
  h: 0.7
  omega_m: 0.3
  omega_lambda: 0.7
  z_lens: 0.3

pipeline:
  z_ids: [50, 100, 150, 200, 250, 300, 350, 400, 417, 450]
  data_path: data/processed/deflection_maps.h5
  
correlator:
  rmin: 0.1  # arcsec
  rmax: 100
  nbins: 20
  
fitter:
  mass_prior: [1e12, 1e15]  # Msun
  c_prior: [1.0, 10.0]
```

**Deliverables**:
- `scripts/run_pipeline.py` ✓
- `lensing/pipeline/__init__.py` (orchestrator) ✓
- `data/config/thesis.yaml` ✓
- Example run on synthetic z_ids [50]

**Time**: 3-4 days  
**Owner**: 1 person (integration)

---

#### **3b: End-to-End Integration Tests (Days 44-48)**

**Tasks**:
```python
# tests/integration/test_pipeline_e2e.py

def test_full_pipeline_synthetic_single_z():
    """Run complete pipeline on 1 synthetic snapshot."""
    config = load_config('data/config/test.yaml')
    config.pipeline.z_ids = [50]  # Single z
    
    pipeline = Pipeline(config)
    result = pipeline.run(z_ids=[50], output_dir='/tmp/test_output')
    
    # Assertions
    assert 50 in result.jacobians
    assert 50 in result.correlations
    assert 50 in result.fits
    assert result.fits[50].m200 > 1e13  # Physics check
    assert result.fits[50].c > 1 and result.fits[50].c < 10  # Reasonable c

def test_pipeline_multiple_z():
    """Run on multiple z-slices."""
    config = load_config('data/config/test.yaml')
    config.pipeline.z_ids = [50, 100, 150]
    
    pipeline = Pipeline(config)
    result = pipeline.run(z_ids=[50, 100, 150], output_dir='/tmp/test_output')
    
    assert len(result.fits) == 3
    # Fit results should be roughly consistent across z
    # (light check: don't constrain too tightly)

def test_pipeline_matches_julia_on_real_data():
    """Validate against Julia on real (small) data."""
    # Use corner patch of real FITS file
    config = load_config('data/config/test.yaml')
    
    # Load real small patch
    alpha_real, z = load_real_corner_patch()
    
    # Compute Python pipeline
    pipeline = Pipeline(config)
    # ... inject real data ...
    
    # Load Julia reference
    with h5py.File('tests/fixtures/reference_outputs.h5', 'r') as f:
        jac_julia = f['jacobian'][:]
        kappa_julia = f['kappa'][:]
    
    # Compare
    assert np.allclose(jac_python, jac_julia, rtol=1e-5)
    assert np.allclose(kappa_python, kappa_julia, rtol=1e-5)
```

**Deliverables**:
- `tests/integration/test_pipeline_e2e.py` ✓ (3-5 full-pipeline tests)
- All integration tests passing
- Output files validated

**Time**: 4 days  
**Owner**: 1 person

---

#### **Phase 3 Summary**

| Task | Days | Owner | Status |
|------|------|-------|--------|
| 3a. Pipeline entry point | 3-4 | Person A | ▢ |
| 3b. Integration tests | 4 | Person B | ▢ |
| **Total** | **7-8 days** | | |

**Definition of done**:
```bash
python scripts/run_pipeline.py data/config/thesis.yaml --z-ids 50
# Should complete in <10 seconds, produce results/

pytest tests/integration/ -v --tb=short
# All ~5 integration tests pass
```

---

### **Phase 4: Plotting & Visualization (Week 6)**

**Goal**: Clean plotting layer (decoupled from physics).

**Tasks**:
1. Create `lensing/plotting/correlation_plots.py`:

```python
class CorrelationPlotter:
    """Matplotlib-based correlation visualization."""
    
    @staticmethod
    def plot_correlation_with_fit(corr: CorrelationResult,
                                   fit: FitResult,
                                   z: float,
                                   ax=None) -> plt.Axes:
        """Plot ξ(r) + NFW fit overlay."""
        if ax is None:
            fig, ax = plt.subplots()
        
        # Data
        ax.errorbar(corr.r_bins, corr.xi, yerr=corr.error,
                    fmt='o', label='Measured', capsize=5)
        
        # Fit
        r_smooth = np.logspace(np.log10(corr.r_bins.min()),
                               np.log10(corr.r_bins.max()), 100)
        xi_fit = fit.predict(r_smooth)
        ax.plot(r_smooth, xi_fit, 'r--', label=f'NFW (M={fit.m200:.2e}, c={fit.c:.1f})')
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('r [arcsec]')
        ax.set_ylabel('ξ(r)')
        ax.legend()
        ax.set_title(f'κ×κ Correlation at z={z:.2f}')
        return ax
```

2. Create `lensing/plotting/map_plots.py` (convergence maps, etc.)

3. Create `scripts/plot_correlations.py`:

```python
#!/usr/bin/env python3
"""Convenience: plot all correlation functions from results."""

if __name__ == '__main__':
    results = load_results('output/pipeline_results.h5')
    plotter = CorrelationPlotter()
    
    for z_id, (corr, fit) in results.items():
        fig, ax = plt.subplots()
        plotter.plot_correlation_with_fit(corr, fit, results.redshifts[z_id], ax=ax)
        plt.savefig(f'output/correlation_z{z_id:04d}.png', dpi=150, bbox_inches='tight')
        plt.close()
```

**Deliverables**:
- `lensing/plotting/` module ✓
- `scripts/plot_correlations.py` ✓
- Example plots generated

**Time**: 2-3 days  
**Owner**: 1 person

---

### **Phase 5: Documentation & Thesis Integration (Week 7)**

**Goal**: Complete documentation package; prepare for thesis.

**Tasks**:
1. Create `docs/physics.md`:
   - All formulas with LaTeX
   - Paper citations for each
   - Thesis references

2. Create `docs/api.md`:
   - Class/function reference
   - Auto-generated from docstrings

3. Create `docs/quickstart.md`:
   - How to run pipeline
   - Example: process single z-slice
   - Interpret output

4. Create `README.md`:
   - Project overview
   - Installation: `pip install -e .`
   - Quick example

**Example structure**:
```markdown
# Weak-Lensing Pipeline

## Installation
```bash
git clone <repo>
cd lensing-pipeline
pip install -e .
```

## Quick Start
```bash
python scripts/run_pipeline.py data/config/thesis.yaml --z-ids 50
```

## References
- Bartelmann & Schneider (2001) - Weak lensing theory
- Navarro, Frenk & White (1997) - NFW profile
- See `docs/physics.md` for complete formula compendium.

---
```

**Deliverables**:
- `docs/` complete ✓
- `README.md` ✓
- All docstrings finalized ✓
- Sphinx/MkDocs setup (if desired) ✓

**Time**: 3-4 days  
**Owner**: 1 person

---

### **Phase 6: Validation & Performance Tuning (Week 8)**

**Goal**: Regression test against Julia; optimize hotspots; prepare for production.

**Tasks**:
1. Create `scripts/validate_against_julia.jl`:
   - Run Julia pipeline on test data
   - Save reference outputs
   - Python script compares

2. Run full regression:
```bash
python scripts/validate_against_julia.py --compare-to tests/fixtures/reference_outputs.h5
# Output: PASS/FAIL for each component
```

3. Profiling:
```bash
python -m cProfile -s cumulative scripts/run_pipeline.py ... | head -20
# Identify remaining bottlenecks
```

4. Optional optimizations:
   - GPU acceleration (JAX) for Jacobian (if needed)
   - C++ extensions for hot loops (if profiling shows >50% time spent)

**Deliverables**:
- Regression test suite ✓
- Performance report ✓
- Optimization recommendations ✓

**Time**: 3-4 days  
**Owner**: 1 person

---

## Implementation Timeline Summary

```
Week 1-2 (Phase 1: Core Physics)
├─ Mon: Setup, Cosmology config
├─ Tue-Wed: NFW consolidation
├─ Thu-Fri: Synthetic test data + Julia baseline
├─ Week 2 Mon-Wed: Jacobian (Numba)
└─ Thu-Fri: Observable extraction + tests
   Status: ✓ 60+ tests, regression baseline ready

Week 3-4 (Phase 2: I/O & Integration)
├─ Mon-Wed: I/O layer + HDF5 conversion
├─ Thu-Fri: Data loader
├─ Week 4 Mon-Tue: Correlator (treecorr)
├─ Wed-Thu: Fitter
└─ Fri: Integration glue
   Status: ✓ All data pipelines working

Week 5 (Phase 3: Orchestration)
├─ Mon-Tue: Main entry point + config system
├─ Wed-Thu: Full integration tests
└─ Fri: End-to-end validation on synthetic z=[50]
   Status: ✓ Full pipeline runs end-to-end

Week 6 (Phase 4: Plotting)
├─ Mon-Tue: Correlation plotting
├─ Wed: Map visualization
└─ Thu-Fri: Example plots + scripts
   Status: ✓ Publication-quality figures

Week 7 (Phase 5: Documentation)
├─ Mon-Tue: Physics formulas + citations
├─ Wed: API documentation
├─ Thu: Quickstart + README
└─ Fri: Polish + final review
   Status: ✓ Thesis-ready documentation

Week 8 (Phase 6: Validation)
├─ Mon-Tue: Julia regression tests
├─ Wed-Thu: Performance profiling
├─ Fri: Final optimization pass
   Status: ✓ Production-ready
```

---

## Critical Path & Dependencies

**Critical chain** (blocks everything else):
1. Phase 0 (setup) → Phase 1a (config) → Phase 1b (NFW) → Phase 1d (Jacobian) → Phase 2a (I/O) → Phase 3a (orchestration)

**Can parallelize**:
- Phase 1c (synthetic data) alongside 1a/1b
- Phase 2a/2b/2c alongside Phase 1d/1e
- Phase 4 anytime after Phase 3a
- Phase 5 anytime after Phase 1b

**Estimated team composition**:
- **Person A** (physics): Cosmology, NFW, Jacobian, Fitter (~25 days)
- **Person B** (data/I/O): Synthetic data, loader, HDF5 (~15 days)
- **Person C** (integration): Correlator, plots, orchestration (~20 days)
- **Total**: ~4-5 weeks elapsed time (8 weeks wall time with parallelization)

---

## Definition of Done (Cumulative Checkpoints)

### End of Phase 1
- [ ] `pytest tests/unit/ -v` → All pass (60+ tests)
- [ ] Physics validated against Julia reference (H5 comparison)
- [ ] Coverage >85% (`pytest --cov=lensing`)
- [ ] Docstrings complete with citations

### End of Phase 2
- [ ] Data loads correctly from HDF5
- [ ] Correlation computation matches treecorr examples
- [ ] NFW fitting recovers injected parameters
- [ ] Integration tests passing (10-15)

### End of Phase 3
- [ ] `python scripts/run_pipeline.py data/config/thesis.yaml --z-ids 50 --output /tmp/test`
- [ ] Output contains: `jacobians/`, `correlations/`, `fits/`
- [ ] Full pipeline runs in <2 min (64×64 map)
- [ ] Results match Julia to ~1e-5 relative error

### End of Phase 4
- [ ] Plots are publication-quality (high DPI, good labels, legends)
- [ ] Can reproduce all figures from thesis with one script

### End of Phase 5
- [ ] `pip install -e .` works cleanly
- [ ] `pytest tests/ -v` → All pass (100+ tests)
- [ ] All docstrings have LaTeX math + references
- [ ] README + quickstart complete

### End of Phase 6
- [ ] Regression tests passing (Julia vs Python <0.1% difference)
- [ ] Performance baseline documented
- [ ] Code review clearance
- [ ] Ready for public release

---

## Risks & Mitigation

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|-----------|
| Numba not matching Julia speed | Medium | Medium | Fallback to scipy.ndimage.convolve1d |
| HDF5 migration loses data | Low | High | Validate all conversions; keep Fortran backup |
| treecorr performance not good | Low | Low | Fall back to scipy.spatial.cKDTree |
| Julia environment unavailable | Low | Medium | Pre-generate reference outputs once early |
| Scope creep on plotting | High | Low | Prioritize minimal plots (thesis-only) |
| Underestimated Jacobian complexity | Medium | Medium | Allocate extra time in Phase 1d |

---

## Success Metrics

**Code quality**:
- [ ] 100+ tests, >85% coverage
- [ ] Zero warnings from mypy (type checking)
- [ ] All formulas cited

**Performance**:
- [ ] 64×64 Jacobian: <100ms
- [ ] 64×64 full pipeline: <10s
- [ ] 1024×1024 Jacobian: <5s

**Reproducibility**:
- [ ] Python results match Julia to ~1e-5 (relative error)
- [ ] Anyone can run `pip install -e . && python scripts/run_pipeline.py`
- [ ] Thesis figures regeneratable from scripts

**Documentation**:
- [ ] Every formula has a source citation
- [ ] API complete + searchable
- [ ] Quickstart takes <30 min to complete

---

## Post-Project Ideas (Future)

1. **GPU acceleration**: JAX for Jacobian computation (10-100× speedup)
2. **Batch processing**: Run multiple z-slices in parallel (embarrassingly parallel)
3. **Web app**: Interactive visualization of correlations + fitting
4. **Publication**: Release as PyPI package + Zenodo DOI
5. **Extensibility**: Plug in alternative lensing profiles (Einasto, truncated NFW, etc.)

---

This roadmap is **living**: adjust phases/timelines as you discover implementation details. Mark completed items as you go!

