# Laurant Thesis Codebase Modernization Plan
## Python Conversion & Production Pipeline Implementation

**Date**: February 17, 2026  
**Status**: Planning Phase - All 12 Issues Approved  
**Goal**: Convert scattered Julia thesis code to maintainable Python production pipeline

---

## EXECUTIVE SUMMARY

### What This Codebase Does

This is a **weak gravitational lensing analysis pipeline** for cosmological simulations. It:
1. Processes deflection fields from HAGN (Horizon-AGN) N-body simulations
2. Computes lensing observables: convergence (κ), shear (γ), flexion (F, G)
3. Measures tangential correlations around galaxies (galaxy-galaxy lensing)
4. Fits NFW (Navarro-Frenk-White) halo profiles to infer masses
5. Generates publication-quality plots for thesis

### Current State Assessment

**Strengths:**
- Core physics implementation appears sound (produces reasonable results)
- Handles large datasets (20k×20k pixel deflection fields)
- Parallel processing capability (Distributed.jl)

**Critical Issues (12 identified, all approved for fix):**
- **Architecture**: Monolithic 1734-line toolkit, external C dependency, no spatial decomposition, coupled viz/compute
- **Code Quality**: Massive DRY violations (2500 duplicate lines), zero error handling, undocumented formulas, dead code
- **Performance**: No test mode, no validation vs. TreeCorr, unclear goal (research vs. production), inefficient jackknife

**Technical Debt**: ~50KB dead code, no unit tests, no citations, hardcoded paths, silent failures

---

## MODERN PYTHON LIBRARY RESEARCH

### Core Scientific Stack

**1. NumPy (numpy.org) - Numerical arrays**
```python
import numpy as np
# Replaces: Julia base arrays
# Use for: All numerical operations, FFTs, linear algebra
# Version: ≥1.24 (has better typing, performance)
```

**2. SciPy (scipy.org) - Scientific algorithms**
```python
from scipy import ndimage, interpolate, optimize, special
# Replaces: Julia Interpolations.jl, Roots.jl, Polynomials.jl
# Use for: Image processing (ndimage), curve fitting (optimize.curve_fit)
# Key functions: ndimage.map_coordinates() for lensing distortion
```

**3. Astropy (astropy.org) - Astronomy utilities**
```python
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy.io import fits
# Replaces: Cosmology.jl, FITSIO.jl, libsoftlens.so
# Use for: Distance calculations, FITS I/O, coordinate transforms
# Example: cosmo = FlatLambdaCDM(H0=70, Om0=0.3); cosmo.angular_diameter_distance(z)
```

### Weak Lensing Specialized Tools

**4. TreeCorr (github.com/rmjarvis/TreeCorr) - CRITICAL**
```python
import treecorr
# Gold standard for correlation functions (DES, HSC, KiDS, LSST all use this)
# Replaces: comp_gs_corr() custom implementation
# Features: Optimized spatial tree, automatic jackknife, GPU support
# References: Jarvis, Bernstein, Jain (2004) MNRAS 352, 338
cat = treecorr.Catalog(ra=ra, dec=dec, g1=gamma1, g2=gamma2, ra_units='deg')
ng = treecorr.NGCorrelation(min_sep=0.05, max_sep=5.0, nbins=50, 
                             sep_units='arcmin', var_method='jackknife')
ng.process(cat_lens, cat_source)
gamma_t = ng.xi  # Tangential shear with proper errors
```

**5. HEALPix (healpy.readthedocs.io) - Spatial decomposition**
```python
import healpy as hp
# Industry standard for spherical pixelization (Planck, CMB analyses)
# Replaces: Custom tiling, enables hierarchical testing
# Use for: NSIDE=16 (3k pixels) → NSIDE=2048 (50M pixels) scaling
# Key functions: hp.ang2pix(), hp.get_all_neighbours(), hp.ud_grade()
npix = hp.nside2npix(NSIDE)  # Equal-area pixels
ipix = hp.ang2pix(NSIDE, ra, dec, lonlat=True)
```

**6. LensTools (github.com/apetri/LensTools) - Lensing toolkit**
```python
from lenstools import ConvergenceMap, ShearMap
# Provides: κ/γ map manipulation, power spectra, ray-tracing
# Use for: Validation, complementary analysis
# May use for: Rebinning, smoothing operations
```

### Workflow Orchestration (Production Pipeline)

**7. Snakemake (snakemake.readthedocs.io) - RECOMMENDED**
```python
# Workflow: Snakefile defines dependency graph
rule compute_observables:
    input: "deflections/{id}.bin"
    output: "observables/kappa_{id}.fits"
    shell: "python compute_observables.py {input} {output}"

rule correlations:
    input: expand("observables/kappa_{id}.fits", id=IDS)
    output: "correlations/tangential_shear.npz"
    script: "scripts/compute_correlations.py"
```
- **Why Snakemake**: Pythonic syntax, automatic parallelization, cluster integration
- **Alternatives**: Nextflow (more DSL-like), Luigi (older, more boilerplate)
- **Features**: Automatic re-run on parameter changes, built-in provenance, conda integration

### Testing & Validation

**8. pytest (pytest.org) - Unit testing**
```python
import pytest
import numpy.testing as npt

def test_nfw_convergence_normalization():
    """NFW convergence should integrate to enclosed mass"""
    kappa = nfw_convergence(r, ks=1.0, rs=100)
    enclosed_mass = 2 * np.pi * np.trapz(r * kappa, r)
    npt.assert_allclose(enclosed_mass, expected_mass, rtol=0.01)
```

**9. Hypothesis (hypothesis.readthedocs.io) - Property testing**
```python
from hypothesis import given, strategies as st

@given(st.floats(min_value=0.1, max_value=10.0))
def test_jacobian_determinant_positive(redshift):
    """Jacobian determinant should always be positive (no multiple images)"""
    alpha = generate_test_deflection(z=redshift)
    jac = compute_jacobian(alpha)
    det = (1 - jac[..., 0]) * (1 - jac[..., 3]) - jac[..., 1] * jac[..., 2]
    assert np.all(det > 0), "Negative determinant indicates multiple imaging"
```

### Visualization

**10. Matplotlib + Seaborn**
```python
import matplotlib.pyplot as plt
import seaborn as sns
# Standard scientific plotting
# Replaces: Plots.jl

from matplotlib import rc
rc('text', usetex=True)  # LaTeX rendering for publications
```

---

## APPROVED ISSUES & SOLUTIONS

### ARCHITECTURE (Issues #1-4)

#### Issue #1: Monolithic 1734-line toolkit → Modular Python structure
**Approved Solution**: Split into 6 focused modules

**Implementation:**
```
cosmo_lensing/
├── __init__.py
├── io.py              # read_deflection_field(), write_fits_map(), load_catalog()
├── derivatives.py     # compute_jacobian(), finite_difference_2d()
├── observables.py     # convergence(), shear(), flexion_F(), flexion_G()
├── correlations.py    # TangentialCorrelation class, wraps TreeCorr
├── nfw.py            # NFWProfile class, mass_concentration_conversions()
└── cosmology.py      # CosmologyCalculator (wraps astropy)
```

**Tasks:**
- [ ] Create package structure with proper `__init__.py`
- [ ] Port I/O functions (read_bin → io.read_deflection_field)
- [ ] Port derivatives (alpha2jac → derivatives.compute_jacobian)
- [ ] Port observables (jac2kappa → observables.convergence)
- [ ] Write docstrings with math formulas + citations
- [ ] Add type hints (Python 3.10+ style)

---

#### Issue #2: libsoftlens.so dependency → Pure Python with astropy
**Approved Solution**: Eliminate C library, use astropy.cosmology

**Implementation:**
```python
# cosmology.py
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u

class CosmologyCalculator:
    """Cosmological distance calculations using Astropy."""
    
    def __init__(self, H0=70, Om0=0.3, Tcmb0=2.725):
        self.cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Tcmb0=Tcmb0)
    
    def angular_diameter_distance(self, z1, z2=0):
        """Angular diameter distance between redshifts [Mpc]"""
        if z2 == 0:
            return self.cosmo.angular_diameter_distance(z1).to(u.Mpc).value
        # For z1 < z2, use angular_diameter_distance_z1z2
        return self.cosmo.angular_diameter_distance_z1z2(z1, z2).to(u.Mpc).value
    
    def critical_density(self, z):
        """Critical density at redshift z [Msun/Mpc³]"""
        return self.cosmo.critical_density(z).to(u.Msun / u.Mpc**3).value
    
    def Omega_m(self, z):
        """Matter density parameter at redshift z"""
        return self.cosmo.Om(z)
```

**Validation task:**
- [ ] Compare astropy vs. libsoftlens for 100 test redshifts
- [ ] Assert agreement within 0.1% (numerical precision)
- [ ] Document any systematic differences

---

#### Issue #3: No spatial decomposition → HEALPix hierarchy
**Approved Solution**: Implement HEALPix-based processing

**Implementation:**
```python
# spatial.py
import healpy as hp
import numpy as np

class HEALPixDeflectionField:
    """Deflection field stored in HEALPix format."""
    
    def __init__(self, nside, coordsys='C'):
        self.nside = nside
        self.npix = hp.nside2npix(nside)
        self.alpha = np.zeros((self.npix, 2))  # (α₁, α₂) per pixel
    
    @classmethod
    def from_cartesian(cls, alpha_cart, wcs, nside):
        """Convert Cartesian 2D array to HEALPix."""
        # Project WCS grid onto sphere
        # Interpolate onto HEALPix pixels
        ...
    
    def get_patch(self, center_ra, center_dec, radius_deg):
        """Extract small patch around position."""
        vec = hp.ang2vec(center_ra, center_dec, lonlat=True)
        ipix_disk = hp.query_disc(self.nside, vec, np.radians(radius_deg))
        return ipix_disk, self.alpha[ipix_disk]
    
    def compute_derivatives_patch(self, ipix_patch):
        """Compute Jacobian only for pixels in patch."""
        # Use HEALPix neighbors for finite differences
        neighbors = [hp.get_all_neighbours(self.nside, ip) for ip in ipix_patch]
        # Finite difference using neighbor pixels
        ...

# Test hierarchy
TEST_SCALES = {
    'unit': {'nside': 16, 'npix': 3072},      # 3k pixels, < 1 sec
    'debug': {'nside': 64, 'npix': 49152},    # 49k pixels, ~10 sec
    'validation': {'nside': 512, 'npix': 3M}, # 3M pixels, ~2 min
    'production': {'nside': 2048, 'npix': 50M} # 50M pixels, ~10 min
}
```

**Tasks:**
- [ ] Implement HEALPixDeflectionField class
- [ ] Write conversion from Cartesian → HEALPix
- [ ] Implement patch extraction
- [ ] Test on NSIDE=16 synthetic data
- [ ] Validate against full Cartesian computation

---

#### Issue #4: Coupled data/viz → Separate layers
**Approved Solution**: Clean separation into compute/analysis/visualization

**Implementation:**
```
scripts/
├── compute_observables.py    # Pure computation, no plots
├── compute_correlations.py   # Correlation functions, save .npz
├── fit_models.py             # Model fitting, save parameters
└── generate_plots.py         # Read results, make figures

# Example: compute_observables.py
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--deflection', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--observables', nargs='+', 
                        choices=['kappa', 'gamma1', 'gamma2', 'F1', 'F2', 'G1', 'G2'])
    args = parser.parse_args()
    
    # COMPUTE ONLY - no plotting
    alpha, z = io.read_deflection_field(args.deflection)
    jac = derivatives.compute_jacobian(alpha)
    
    for obs in args.observables:
        data = OBSERVABLE_MAP[obs](jac)
        io.write_fits(data, f"{args.output_dir}/{obs}_z{z:.3f}.fits")
    
    logging.info(f"Computed {len(args.observables)} observables, saved to {args.output_dir}")
```

**Tasks:**
- [ ] Separate compute scripts (no matplotlib imports)
- [ ] Separate visualization scripts (read from .npz/.fits)
- [ ] Add CLI interfaces (argparse)
- [ ] Enable headless computation (for HPC)

---

### CODE QUALITY (Issues #5-8)

#### Issue #5: DRY violation - 10 duplicate map scripts → Single parameterized script
**Approved Solution**: One script with observable selection

**Implementation:**
```python
# scripts/generate_maps.py
OBSERVABLES = {
    'kappa': lambda jac: observables.convergence(jac),
    'gamma1': lambda jac: observables.shear_1(jac),
    'gamma2': lambda jac: observables.shear_2(jac),
    'shear_magnitude': lambda jac: np.hypot(observables.shear_1(jac), observables.shear_2(jac)),
    'F1': lambda jac: observables.flexion_F1(jac),
    'F2': lambda jac: observables.flexion_F2(jac),
    'G1': lambda jac: observables.flexion_G1(jac),
    'G2': lambda jac: observables.flexion_G2(jac),
    'rotation': lambda jac: observables.rotation(jac),
}

def generate_maps(deflection_ids, observable_types, output_dir, rebin_factor=5):
    for did in deflection_ids:
        alpha, z = io.read_deflection_field(deflection_file(did))
        jac = derivatives.compute_jacobian(alpha, scale=1.0)
        
        for obs_type in observable_types:
            data = OBSERVABLES[obs_type](jac)
            if rebin_factor > 1:
                data = rebin(data, factor=rebin_factor)
            
            fits_file = f"{output_dir}/{obs_type}_id{did:04d}_z{z:.3f}.fits"
            io.write_fits_map(data, fits_file, redshift=z)
            logging.info(f"Wrote {obs_type} map: {fits_file}")

# CLI: python generate_maps.py --ids 50 100 150 --observables kappa gamma1 gamma2
```

**Tasks:**
- [ ] Delete 9 duplicate rebuild_*_maps.jl files
- [ ] Implement generate_maps.py with argparse CLI
- [ ] Add `--observables all` shortcut
- [ ] Test on single ID first, then batch

---

#### Issue #6: Zero error handling → Explicit validation
**Approved Solution**: Add error handling with custom exceptions

**Implementation:**
```python
# io.py
class DeflectionFieldError(Exception):
    """Raised when deflection field data is invalid."""
    pass

def read_deflection_field(filepath):
    """
    Read Fortran binary deflection field.
    
    Returns:
        alpha: (nx, ny, 2) array of (α₁, α₂) deflections [degrees]
        redshift: Source plane redshift
    
    Raises:
        FileNotFoundError: If file doesn't exist
        DeflectionFieldError: If file format invalid or data corrupted
        PermissionError: If file not readable
    """
    from pathlib import Path
    import logging
    
    filepath = Path(filepath)
    if not filepath.exists():
        logging.error(f"Deflection field not found: {filepath}")
        raise FileNotFoundError(f"Missing: {filepath}")
    
    if filepath.stat().st_size < 1000:  # Minimum valid size
        raise DeflectionFieldError(f"File too small, likely corrupted: {filepath}")
    
    try:
        with FortranFile(filepath, 'r') as f:
            size = f.read_ints(dtype=np.int32)
            if len(size) != 2 or np.any(size <= 0):
                raise DeflectionFieldError(f"Invalid dimensions: {size}")
            
            dummy = f.read_reals(dtype=np.float32)  # Header padding
            redshift = f.read_reals(dtype=np.float64)[0]
            
            alpha = np.zeros((size[0], size[1], 2), dtype=np.float32)
            # Read chunked data...
            
    except (struct.error, EOFError) as e:
        raise DeflectionFieldError(f"Corrupted file {filepath}: {e}")
    
    # Validate data
    if np.any(np.isnan(alpha)) or np.any(np.isinf(alpha)):
        raise DeflectionFieldError(f"NaN/Inf in {filepath}")
    
    if np.abs(alpha).max() > 10.0:  # Sanity check: deflections > 10° unlikely
        logging.warning(f"Large deflections (max={np.abs(alpha).max():.2f}°) in {filepath}")
    
    logging.info(f"Loaded deflection: {size} pixels, z={redshift:.3f}, "
                 f"α_rms={np.std(alpha):.4f}°")
    return alpha, redshift
```

**Tasks:**
- [ ] Define custom exception classes
- [ ] Add error handling to all I/O functions
- [ ] Add data validation (NaN/Inf checks, range checks)
- [ ] Add logging statements (info, warning, error levels)
- [ ] Write tests for error conditions

---

#### Issue #7: Undocumented formulas → Citations + unit tests
**Approved Solution**: Add paper references and analytic test cases

**Implementation:**
```python
# observables.py
def flexion_F1(jacobian):
    """
    Compute first flexion component F₁ (arcmin⁻¹).
    
    Flexion describes third-order lensing distortions beyond shear.
    
    Formula (Bacon et al. 2006, Eq. 8):
        F₁ = ½(∂²α₁/∂θ₁² - ∂²α₂/∂θ₁∂θ₂) + ½(∂²α₁/∂θ₂² + ∂²α₂/∂θ₁∂θ₂)
    
    Args:
        jacobian: (nx, ny, 12) array from compute_jacobian()
                  [4:8] contains second derivatives
    
    Returns:
        F1: (nx, ny) array of first flexion component
    
    References:
        Bacon, D. J., Goldberg, D. M., Rowe, B. T. P., & Taylor, A. N. (2006).
        "Weak gravitational flexion", MNRAS, 365, 414.
        https://doi.org/10.1111/j.1365-2966.2005.09624.x
        
        Schneider, P., & Er, X. (2008).
        "Weak lensing goes bananas: what flexion really measures", 
        A&A, 485, 363. https://doi.org/10.1051/0004-6361:20078631
    """
    d11a1 = jacobian[:, :, 4]   # ∂²α₁/∂θ₁²
    d22a1 = jacobian[:, :, 6]   # ∂²α₁/∂θ₂²
    d12a2 = jacobian[:, :, 11]  # ∂²α₂/∂θ₁∂θ₂
    
    return 0.5 * (d11a1 - d12a2) + 0.5 * (d22a1 + d12a2)

# tests/test_observables.py
def test_point_mass_convergence():
    """
    Convergence for point mass lens (Einstein radius θ_E).
    
    Analytic solution (Schneider, Ehlers, Falco 1992, Eq. 4.15):
        κ(θ) = θ_E²/(2θ²) * H(θ_E - θ)
    where H is Heaviside function (1 inside Einstein radius, 0 outside).
    """
    theta_E = 1.0  # arcsec
    nx, ny = 200, 200
    theta = np.linspace(0.01, 5.0, nx)
    
    # Generate point mass deflection field
    alpha = generate_point_mass_deflection(theta_E, grid_size=nx)
    jac = derivatives.compute_jacobian(alpha)
    kappa_numeric = observables.convergence(jac)
    
    # Analytic solution
    r = np.hypot(*np.meshgrid(theta - theta.mean(), theta - theta.mean()))
    kappa_analytic = np.where(r < theta_E, theta_E**2 / (2 * r**2), 0)
    
    # Should agree to 1% (limited by finite differencing)
    np.testing.assert_allclose(kappa_numeric, kappa_analytic, rtol=0.01)

def test_nfw_flexion_symmetry():
    """
    For circularly symmetric NFW lens, flexion should be tangential.
    
    Test: F_radial ≈ 0, G obeys φ → -φ symmetry.
    """
    # Generate NFW deflection (ks=1.0, rs=100 kpc)
    alpha_nfw = generate_nfw_deflection(ks=1.0, rs=100, grid_size=500)
    jac = derivatives.compute_jacobian(alpha_nfw)
    F1 = observables.flexion_F1(jac)
    G1 = observables.flexion_G1(jac)
    
    # Convert to polar: F_r, F_φ
    x, y = np.meshgrid(np.arange(500) - 250, np.arange(500) - 250)
    phi = np.arctan2(y, x)
    F_radial = F1 * np.cos(phi)  # Simplified
    
    # Radial flexion should be small for spherical lens
    assert np.abs(F_radial).max() < 0.01 * np.abs(F1).max()
    
    # G should be symmetric: G(φ) = G(-φ)
    np.testing.assert_allclose(G1[:, :250], G1[:, 250::-1], rtol=0.05)
```

**Tasks:**
- [ ] Add docstrings with formulas + citations to all functions
- [ ] Implement synthetic data generators (point mass, SIS, NFW)
- [ ] Write 10+ analytic test cases
- [ ] Run tests in CI (GitHub Actions)

---

#### Issue #8: Dead SIS code → Delete + refactor correlations
**Approved Solution**: Remove 50KB dead code, proper correlation class

**Implementation:**
```python
# correlations.py
import treecorr
import logging

class TangentialCorrelation:
    """
    Tangential correlation function wrapper around TreeCorr.
    
    Provides unified interface for γ_t, κ, F_t, G_t correlations.
    """
    
    def __init__(self, rmin=0.05, rmax=5.0, nbins=50, sep_units='arcmin',
                 bin_scheme='log', var_method='jackknife', n_patches=50):
        """
        Args:
            rmin, rmax: Separation range
            nbins: Number of radial bins
            sep_units: 'arcmin', 'arcminutes', 'degrees'
            bin_scheme: 'log' or 'linear'
            var_method: 'jackknife', 'bootstrap', or None
            n_patches: Number of jackknife/bootstrap patches
        """
        self.config = {
            'min_sep': rmin,
            'max_sep': rmax,
            'nbins': nbins,
            'sep_units': sep_units,
            'bin_slop': 0.1,  # Trade accuracy for speed
        }
        if bin_scheme == 'log':
            self.config['bin_type'] = 'Log'
        
        self.var_method = var_method
        self.n_patches = n_patches
    
    def compute(self, ra_lens, dec_lens, ra_source, dec_source,
                g1, g2, weights=None):
        """
        Compute tangential shear/convergence around lenses.
        
        Args:
            ra_lens, dec_lens: Lens positions [degrees]
            ra_source, dec_source: Source positions [degrees]
            g1, g2: Shear/convergence components at source positions
            weights: Optional weights per source
        
        Returns:
            result: dict with keys 'r', 'xi', 'xi_err', 'npairs'
        """
        # Lens catalog (points)
        cat_lens = treecorr.Catalog(
            ra=ra_lens, dec=dec_lens,
            ra_units='deg', dec_units='deg'
        )
        
        # Source catalog (with shear/convergence field)
        cat_source = treecorr.Catalog(
            ra=ra_source, dec=dec_source,
            g1=g1, g2=g2, w=weights,
            ra_units='deg', dec_units='deg',
            patch_centers=self._get_patch_centers() if self.var_method else None
        )
        
        # Compute NG correlation (lens-source)
        ng = treecorr.NGCorrelation(**self.config, var_method=self.var_method)
        ng.process(cat_lens, cat_source)
        
        result = {
            'r': np.exp(ng.meanlogr),  # Mean separation per bin
            'xi': ng.xi,               # Tangential component
            'xi_err': np.sqrt(ng.varxi) if self.var_method else None,
            'npairs': ng.npairs,
            'weight': ng.weight
        }
        
        logging.info(f"Computed correlation: {len(ra_lens)} lenses, "
                     f"{len(ra_source)} sources, {ng.npairs.sum()} pairs")
        return result
    
    def _get_patch_centers(self):
        """Generate spatial patches for jackknife/bootstrap."""
        # Use k-means clustering for irregular patches
        from sklearn.cluster import KMeans
        # ... implementation
        return patch_centers

# models.py
class NFWProfile:
    """
    Navarro-Frenk-White halo density profile.
    
    References:
        Navarro, J. F., Frenk, C. S., & White, S. D. M. (1997).
        "A Universal Density Profile from Hierarchical Clustering",
        ApJ, 490, 493. https://doi.org/10.1086/304888
        
        Wright, C. O., & Brainerd, T. G. (2000).
        "Gravitational Lensing by NFW Halos", ApJ, 534, 34.
        https://doi.org/10.1086/308744
    """
    
    def __init__(self, ks, rs):
        """
        Args:
            ks: Convergence scale parameter [dimensionless]
            rs: Scale radius [kpc or arcsec, specify in runit]
        """
        self.ks = ks
        self.rs = rs
    
    def convergence(self, r):
        """Convergence κ(r) for NFW profile (Wright & Brainerd 2000, Eq. 11)."""
        x = r / self.rs
        f = self._F_function(x)
        return 2 * self.ks * self.rs**2 * (1 - f) / (r**2 - self.rs**2)
    
    def shear_tangential(self, r):
        """Tangential shear γ_t(r) (Wright & Brainerd 2000, Eq. 12)."""
        # ... implementation
    
    @staticmethod
    def _F_function(x):
        """Auxiliary function for NFW lensing (Wright & Brainerd 2000, Eq. 13)."""
        if x > 1:
            s = np.sqrt(x**2 - 1)
            return np.arctan(s) / s
        elif x < 1:
            s = np.sqrt(1 - x**2)
            return np.arctanh(s) / s
        else:
            return 1.0
```

**Tasks:**
- [ ] DELETE correlation_plots_sis_model.jl (25KB)
- [ ] DELETE correlation_plots_sis_model_2.jl (25KB)
- [ ] Implement TangentialCorrelation class
- [ ] Implement NFWProfile model class
- [ ] Test against current comp_gs_corr() output
- [ ] Validate TreeCorr vs. custom implementation (Issue #10)

---

### PERFORMANCE (Issues #9-12)

#### Issue #9: No test mode → Hierarchical test pyramid
**Approved Solution**: Unit/debug/validation/production scales

**Implementation:**
```python
# config.py
TEST_CONFIGS = {
    'unit': {
        'description': 'Unit tests with synthetic data',
        'deflection_size': 100,      # 100×100 pixels
        'n_galaxies': 5,
        'n_redshift_slices': 1,
        'runtime_target': '< 1 second',
        'data_type': 'synthetic',    # Generated NFW
        'use_case': 'pytest, CI/CD'
    },
    'debug': {
        'description': 'Development iteration',
        'deflection_size': 1000,     # 1000×1000 = 1M pixels
        'n_galaxies': 50,
        'n_redshift_slices': 2,
        'runtime_target': '< 10 seconds',
        'data_type': 'real_patch',   # Extract from full field
        'use_case': 'Local development'
    },
    'validation': {
        'description': 'Pre-production validation',
        'deflection_size': 5000,     # 25M pixels
        'n_galaxies': 500,
        'n_redshift_slices': 5,
        'runtime_target': '< 2 minutes',
        'data_type': 'real',
        'use_case': 'Pre-HPC validation'
    },
    'production': {
        'description': 'Full science run',
        'deflection_size': 20000,    # 400M pixels
        'n_galaxies': None,          # All available
        'n_redshift_slices': 10,
        'runtime_target': '~10 minutes per slice',
        'data_type': 'real',
        'use_case': 'Science results'
    }
}

# synthetic_data.py
def generate_nfw_deflection(mass, concentration, redshift, 
                             grid_size=100, fov_deg=1.0):
    """
    Generate synthetic NFW deflection field with known parameters.
    
    Use for unit tests - can compare numerical derivatives to analytic.
    """
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    
    # Convert mass/concentration to lensing parameters
    Dls = cosmo.angular_diameter_distance_z1z2(redshift, redshift*2).value
    # ... NFW deflection angle formula (Wright & Brainerd 2000)
    
    # Create grid
    x = np.linspace(-fov_deg/2, fov_deg/2, grid_size)
    xx, yy = np.meshgrid(x, x)
    r = np.hypot(xx, yy)
    
    # NFW deflection (radial)
    alpha_r = nfw_deflection_angle(r, ks, rs)
    alpha_x = alpha_r * xx / r
    alpha_y = alpha_r * yy / r
    
    alpha = np.stack([alpha_x, alpha_y], axis=-1)
    return alpha, redshift

# CLI usage:
# pytest tests/ --config unit  # Run all tests with synthetic data
# python analyze.py --config debug --patch 0,0,1000  # Extract 1000×1000 patch
# python analyze.py --config validation  # Full validation run
# snakemake --config scale=production  # Production pipeline
```

**Tasks:**
- [ ] Implement synthetic NFW/SIS/point-mass generators
- [ ] Add `--config` flag to all scripts
- [ ] Implement patch extraction for real data
- [ ] Test full pipeline on each scale
- [ ] Document runtime/memory for each scale

---

#### Issue #10: No validation vs. established tools → TreeCorr cross-check
**Approved Solution**: Automated tests against TreeCorr

**Implementation:**
```python
# tests/test_against_treecorr.py
import pytest
import numpy as np
import treecorr
from cosmo_lensing import correlations, io, derivatives, observables

@pytest.fixture
def test_deflection():
    """Load or generate test deflection field."""
    # Use synthetic NFW with known parameters
    from cosmo_lensing.synthetic_data import generate_nfw_deflection
    return generate_nfw_deflection(mass=1e14, concentration=5, 
                                    redshift=0.5, grid_size=500)

@pytest.fixture
def test_galaxies():
    """Generate test galaxy positions."""
    np.random.seed(42)
    n_gal = 100
    ra = np.random.uniform(-0.5, 0.5, n_gal)
    dec = np.random.uniform(-0.5, 0.5, n_gal)
    return ra, dec

def test_tangential_shear_vs_treecorr(test_deflection, test_galaxies):
    """
    Our tangential shear should match TreeCorr's within 1%.
    
    This is a CRITICAL validation test.
    """
    alpha, z = test_deflection
    ra_lens, dec_lens = test_galaxies
    
    # OUR implementation
    jac = derivatives.compute_jacobian(alpha)
    gamma1 = observables.shear_1(jac)
    gamma2 = observables.shear_2(jac)
    
    # Create source grid (same as deflection grid)
    nx, ny = alpha.shape[:2]
    x = np.linspace(-0.5, 0.5, nx)
    ra_grid, dec_grid = np.meshgrid(x, x)
    
    our_corr = correlations.TangentialCorrelation(
        rmin=0.05, rmax=1.0, nbins=20, sep_units='deg'
    )
    our_result = our_corr.compute(
        ra_lens, dec_lens,
        ra_grid.flatten(), dec_grid.flatten(),
        gamma1.flatten(), gamma2.flatten()
    )
    
    # TREECORR implementation (direct)
    cat_lens = treecorr.Catalog(ra=ra_lens, dec=dec_lens, ra_units='deg', dec_units='deg')
    cat_source = treecorr.Catalog(
        ra=ra_grid.flatten(), dec=dec_grid.flatten(),
        g1=gamma1.flatten(), g2=gamma2.flatten(),
        ra_units='deg', dec_units='deg'
    )
    
    ng = treecorr.NGCorrelation(min_sep=0.05, max_sep=1.0, nbins=20, sep_units='deg')
    ng.process(cat_lens, cat_source)
    
    treecorr_result = ng.xi
    
    # Should agree within 1% (allows for minor implementation differences)
    np.testing.assert_allclose(our_result['xi'], treecorr_result, rtol=0.01,
                                err_msg="Tangential shear differs from TreeCorr")
    
    print(f"✓ Tangential shear matches TreeCorr (max diff: "
          f"{np.abs(our_result['xi'] - treecorr_result).max():.4f})")

def test_convergence_azimuthal_average_vs_treecorr(test_deflection, test_galaxies):
    """Azimuthally averaged convergence should also match."""
    # Similar test but for κ(r) instead of γ_t(r)
    ...

def test_against_published_des_profiles():
    """
    Load published DES Y3 tangential shear profiles (Amon et al. 2022).
    
    Our pipeline should produce similar shapes/amplitudes for comparable
    lens samples (not exact match since different simulations, but ballpark).
    """
    # Load DES Y3 data from paper
    # Run our pipeline on comparable mass bin
    # Check shape is similar (power-law at small r, flattening at large r)
    ...
```

**Tasks:**
- [ ] Implement TreeCorr comparison tests
- [ ] Run tests on 3 different synthetic profiles (point mass, SIS, NFW)
- [ ] Document any systematic differences (if > 1%)
- [ ] Add tests to CI pipeline
- [ ] Optionally: compare to published DES/HSC/KiDS profiles

---

#### Issue #11: Research vs. production → Production pipeline with Snakemake
**Approved Solution**: Workflow orchestration for reproducibility

**Implementation:**
```python
# workflow/Snakefile
configfile: "config.yaml"

# Extract config
DEFLECTION_IDS = config['deflection_ids']  # [50, 100, 150, ...]
OBSERVABLES = config['observables']        # ['kappa', 'gamma1', 'gamma2']
OUTPUT_DIR = config['output_dir']

# Rule: compute all observables for all deflection fields
rule all:
    input:
        # Observable maps
        expand("{outdir}/observables/{obs}_id{id:04d}.fits",
               outdir=OUTPUT_DIR, obs=OBSERVABLES, id=DEFLECTION_IDS),
        # Correlation functions
        f"{OUTPUT_DIR}/correlations/tangential_shear.npz",
        # Model fits
        f"{OUTPUT_DIR}/fits/nfw_parameters.csv",
        # Plots
        expand("{outdir}/plots/{obs}_correlation.png",
               outdir=OUTPUT_DIR, obs=['kappa', 'gamma'])

# Rule: read deflection field and compute observables
rule compute_observables:
    input:
        deflection = lambda wildcards: get_deflection_path(wildcards.id)
    output:
        expand("{{outdir}}/observables/{obs}_id{{id}}.fits", obs=OBSERVABLES)
    params:
        observables = OBSERVABLES
    threads: 4  # Parallel processing within job
    resources:
        mem_mb = 16000  # 16GB RAM
    log:
        "logs/observables_{id}.log"
    script:
        "scripts/compute_observables.py"

# Rule: compute correlation functions from observables + galaxy catalog
rule compute_correlations:
    input:
        observables = expand("{outdir}/observables/kappa_id{id:04d}.fits",
                              outdir=OUTPUT_DIR, id=DEFLECTION_IDS),
        catalog = config['galaxy_catalog']
    output:
        "{outdir}/correlations/tangential_shear.npz"
    params:
        rmin = config['correlation']['rmin'],
        rmax = config['correlation']['rmax'],
        nbins = config['correlation']['nbins']
    threads: 8
    resources:
        mem_mb = 32000
    log:
        "logs/correlations.log"
    script:
        "scripts/compute_correlations.py"

# Rule: fit NFW profiles to correlations
rule fit_models:
    input:
        "{outdir}/correlations/tangential_shear.npz"
    output:
        "{outdir}/fits/nfw_parameters.csv",
        "{outdir}/fits/fit_quality.png"
    log:
        "logs/fitting.log"
    script:
        "scripts/fit_nfw_profiles.py"

# Rule: generate publication plots
rule generate_plots:
    input:
        correlations = "{outdir}/correlations/tangential_shear.npz",
        fits = "{outdir}/fits/nfw_parameters.csv"
    output:
        expand("{{outdir}}/plots/{obs}_correlation.png", obs=['kappa', 'gamma'])
    script:
        "scripts/generate_plots.py"

# Helper function
def get_deflection_path(deflection_id):
    return f"/data/deflections/test2_deflection_ipl_{deflection_id:04d}_propage.bin"
```

```yaml
# config.yaml
# Deflection field IDs to process
deflection_ids: [50, 100, 150, 200, 250, 300, 350, 400, 417, 450]

# Which observables to compute
observables: ['kappa', 'gamma1', 'gamma2', 'F1', 'F2', 'G1', 'G2']

# Galaxy catalog
galaxy_catalog: "Data/Galaxies_0-6_lensed.v2.0_cut_i27.fits"

# Correlation function parameters
correlation:
  rmin: 0.05     # arcmin
  rmax: 5.0      # arcmin
  nbins: 50
  bin_scheme: 'log'
  var_method: 'jackknife'
  n_patches: 50

# Output directory (timestamped for reproducibility)
output_dir: "results/run_20260217_183000"

# Computational resources
threads_per_job: 4
max_jobs: 10  # Parallel jobs on cluster
```

**Usage:**
```bash
# Dry run (see what would be executed)
snakemake --configfile config.yaml -n

# Local execution (use all cores)
snakemake --configfile config.yaml --cores 16

# Cluster execution (SLURM)
snakemake --configfile config.yaml --profile slurm --jobs 50

# Only recompute what changed (if config.yaml modified)
snakemake --configfile config.yaml --cores 16

# Generate workflow diagram
snakemake --configfile config.yaml --dag | dot -Tpng > workflow.png
```

**Tasks:**
- [ ] Install Snakemake (`pip install snakemake`)
- [ ] Create Snakefile with dependency graph
- [ ] Create config.yaml template
- [ ] Test workflow on `--config scale=debug`
- [ ] Add SLURM cluster profile for HPC
- [ ] Document workflow execution
- [ ] Add provenance tracking (git commit, config snapshot)

---

#### Issue #12: Inefficient jackknife → Spatial jackknife with TreeCorr
**Approved Solution**: Proper spatial jackknife, 100× faster

**Implementation:**
```python
# Already implemented in Issue #8 TangentialCorrelation class
# TreeCorr handles jackknife automatically with var_method='jackknife'

# Additional: Covariance matrix export
def compute_covariance_matrix(correlation_results, n_bootstrap=1000):
    """
    Compute full covariance matrix for correlation function bins.
    
    Needed for proper χ² fitting with off-diagonal correlations.
    """
    n_bins = len(correlation_results['r'])
    
    # TreeCorr provides varxi (diagonal errors), but we can also
    # use bootstrap for full covariance
    if 'bootstrap_samples' in correlation_results:
        samples = correlation_results['bootstrap_samples']  # (n_bootstrap, n_bins)
        cov = np.cov(samples.T)
    else:
        # Use diagonal approximation (assumes independent bins)
        cov = np.diag(correlation_results['xi_err']**2)
        logging.warning("Using diagonal covariance (no bootstrap samples)")
    
    return cov

# Fitting with covariance
def fit_nfw_with_covariance(r, xi, cov, initial_params):
    """
    Fit NFW profile with proper χ² using covariance matrix.
    """
    cov_inv = np.linalg.inv(cov)
    
    def chi2(params):
        ks, rs = params
        model = nfw_model(r, ks, rs)
        residual = xi - model
        return residual @ cov_inv @ residual
    
    result = scipy.optimize.minimize(
        chi2, initial_params,
        bounds=[(0.001, 10), (10, 1000)]  # Reasonable ks, rs ranges
    )
    
    # Parameter uncertainties from inverse Hessian
    hessian = scipy.optimize.approx_fprime(result.x, lambda p: chi2(p))
    param_cov = np.linalg.inv(hessian / 2)  # Factor of 2 from χ²
    param_errors = np.sqrt(np.diag(param_cov))
    
    return {
        'params': result.x,
        'errors': param_errors,
        'chi2': result.fun,
        'reduced_chi2': result.fun / (len(r) - 2),
        'success': result.success
    }
```

**Comparison:**
```
OLD (current Julia):
- Method: Bootstrap (misnamed "jackknife")
- N_samples: 5000
- Runtime: ~80% of total (dominant cost)
- Output: Marginal errors only (diagonal)
- Spatial correlations: Ignored

NEW (Python):
- Method: Spatial jackknife (TreeCorr)
- N_patches: 50
- Runtime: < 5% of total (TreeCorr optimized)
- Output: Full covariance matrix
- Spatial correlations: Properly accounted for
- Speed: 100× faster
```

**Tasks:**
- [ ] Remove jackknife loop from correlation_plots.jl equivalent
- [ ] Use TreeCorr's `var_method='jackknife'`
- [ ] Implement covariance matrix computation
- [ ] Update fitting to use covariance (not just marginal errors)
- [ ] Validate: compare old vs. new error estimates
- [ ] Document: spatial jackknife methodology

---

## ORDERED IMPLEMENTATION PLAN

### Phase 1: Foundation (Week 1-2)
**Goal**: Set up Python infrastructure, port core computations

**Tasks:**
1. [ ] Create package structure (`cosmo_lensing/`)
2. [ ] Set up development environment (conda/venv, requirements.txt)
3. [ ] Port I/O functions (`io.py`)
   - [ ] read_deflection_field() with error handling
   - [ ] write_fits_map() with WCS headers
   - [ ] load_galaxy_catalog()
4. [ ] Port cosmology (`cosmology.py`)
   - [ ] CosmologyCalculator class (astropy)
   - [ ] Validation test vs. libsoftlens
5. [ ] Port derivatives (`derivatives.py`)
   - [ ] compute_jacobian() (first derivatives)
   - [ ] compute_higher_order_derivatives() (second derivatives)
6. [ ] Port observables (`observables.py`)
   - [ ] convergence(), shear_1(), shear_2()
   - [ ] flexion_F1(), flexion_F2(), flexion_G1(), flexion_G2()
   - [ ] Add docstrings with citations
7. [ ] Set up pytest infrastructure
   - [ ] tests/ directory
   - [ ] conftest.py with fixtures
   - [ ] pytest.ini configuration

**Deliverable**: Core library that can read deflection, compute observables
**Test**: Process one deflection field, save FITS map

---

### Phase 2: Synthetic Data & Validation (Week 3)
**Goal**: Generate test datasets, validate against analytic solutions

**Tasks:**
8. [ ] Implement synthetic data generators (`synthetic_data.py`)
   - [ ] generate_point_mass_deflection()
   - [ ] generate_sis_deflection()
   - [ ] generate_nfw_deflection()
9. [ ] Write analytic validation tests
   - [ ] test_point_mass_convergence()
   - [ ] test_nfw_convergence_normalization()
   - [ ] test_flexion_symmetry()
   - [ ] test_derivatives_vs_finite_difference()
10. [ ] Implement hierarchical test configs
    - [ ] TEST_CONFIGS dict (unit/debug/validation/production)
    - [ ] CLI flag `--config`
11. [ ] Test full pipeline on synthetic data
    - [ ] Unit scale (100×100, < 1 sec)
    - [ ] Debug scale (1000×1000, < 10 sec)

**Deliverable**: Test suite with 10+ passing tests
**Test**: `pytest tests/ --config unit` completes in < 10 seconds

---

### Phase 3: Correlation Functions & TreeCorr (Week 4)
**Goal**: Port correlation computation, integrate TreeCorr

**Tasks:**
12. [ ] Implement correlation module (`correlations.py`)
    - [ ] TangentialCorrelation class
    - [ ] Wrap TreeCorr for standardization
13. [ ] Implement NFW model (`nfw.py`)
    - [ ] NFWProfile class
    - [ ] mass_concentration_conversions()
    - [ ] Lensing functions (κ, γ_t from Wright & Brainerd 2000)
14. [ ] Write TreeCorr validation tests (`test_against_treecorr.py`)
    - [ ] test_tangential_shear_vs_treecorr()
    - [ ] test_convergence_azimuthal_average()
15. [ ] Port fitting code
    - [ ] fit_nfw_with_covariance()
    - [ ] Use scipy.optimize

**Deliverable**: Correlation functions match TreeCorr within 1%
**Test**: Run on synthetic NFW, compare to analytic γ_t(r)

---

### Phase 4: Production Pipeline (Week 5-6)
**Goal**: Snakemake workflow, batch processing

**Tasks:**
16. [ ] Create Snakefile workflow
    - [ ] Rule: compute_observables
    - [ ] Rule: compute_correlations
    - [ ] Rule: fit_models
    - [ ] Rule: generate_plots
17. [ ] Create config.yaml template
18. [ ] Implement production scripts
    - [ ] scripts/compute_observables.py (CLI)
    - [ ] scripts/compute_correlations.py
    - [ ] scripts/fit_nfw_profiles.py
    - [ ] scripts/generate_plots.py (separate from compute)
19. [ ] Add provenance tracking
    - [ ] Save git commit hash
    - [ ] Save config.yaml snapshot
    - [ ] Compute output checksums
20. [ ] Test workflow
    - [ ] Dry run (`snakemake -n`)
    - [ ] Local execution (debug scale)
    - [ ] Cluster execution (validation scale)

**Deliverable**: Fully automated workflow, reproducible results
**Test**: `snakemake --config scale=validation` completes successfully

---

### Phase 5: Scaling & Optimization (Week 7)
**Goal**: HEALPix decomposition, parallel processing

**Tasks:**
21. [ ] Implement HEALPix support (`spatial.py`)
    - [ ] HEALPixDeflectionField class
    - [ ] Cartesian → HEALPix conversion
    - [ ] Patch extraction
22. [ ] Add parallel processing
    - [ ] Multiprocessing for embarrassingly parallel tasks
    - [ ] Dask for large arrays (optional)
23. [ ] Optimize performance
    - [ ] Profile code (cProfile)
    - [ ] Optimize hotspots (numba JIT if needed)
24. [ ] Memory optimization
    - [ ] Chunked processing for large files
    - [ ] Intermediate cleanup (gc.collect())

**Deliverable**: Can process NSIDE=2048 HEALPix (50M pixels)
**Test**: Production scale run completes in < 1 hour

---

### Phase 6: Final Validation & Documentation (Week 8)
**Goal**: Publication-ready, fully documented

**Tasks:**
25. [ ] Run full production pipeline on all 10 redshift slices
26. [ ] Compare results to original Julia code
    - [ ] Check correlation functions match (< 5% difference)
    - [ ] Check fitted masses match (< 10% difference)
27. [ ] Write comprehensive documentation
    - [ ] README.md with quickstart
    - [ ] API documentation (Sphinx)
    - [ ] Tutorial notebooks (optional)
28. [ ] Code cleanup
    - [ ] Delete all Julia files (backup to `archive/`)
    - [ ] Remove dead code
    - [ ] Final linting (black, isort, mypy)
29. [ ] Create release
    - [ ] Tag version v1.0
    - [ ] Archive on Zenodo (DOI)
    - [ ] Submit code with thesis

**Deliverable**: Publication-ready code + documentation
**Test**: External user can reproduce thesis results from scratch

---

## TEST STRATEGY

### Unit Tests (Fast, < 10 seconds total)
Run on every commit (CI/CD).

**Coverage:**
- [ ] I/O functions (read/write, error handling)
- [ ] Cosmology calculations (distances, critical density)
- [ ] Derivatives (Jacobian, finite differences)
- [ ] Observables (κ, γ, F, G with synthetic data)
- [ ] NFW model functions (analytic limits)

**Test data**: Synthetic 100×100 pixel fields

**Command**: `pytest tests/ --config unit --cov=cosmo_lensing`

**Pass criteria**: All tests pass, > 80% code coverage

---

### Integration Tests (Medium, < 2 minutes)
Run before merging to main branch.

**Coverage:**
- [ ] Full pipeline (deflection → observables → correlation → fit)
- [ ] TreeCorr validation (match within 1%)
- [ ] Error propagation (jackknife errors reasonable)
- [ ] File I/O (FITS headers correct, data not corrupted)

**Test data**: 1000×1000 real data patch OR large synthetic field

**Command**: `pytest tests/ --config debug --integration`

**Pass criteria**: All integration tests pass, correlation functions smooth

---

### Validation Tests (Slow, < 30 minutes)
Run before production runs, weekly during development.

**Coverage:**
- [ ] Full 10 redshift slices (validation scale, not production)
- [ ] Compare to original Julia output (< 5% difference)
- [ ] Check fit quality (reduced χ² ≈ 1)
- [ ] Visual inspection of plots (no artifacts)

**Test data**: Subset of full dataset (5000×5000 pixels, 500 galaxies)

**Command**: `snakemake --config scale=validation --dry-run; snakemake --config scale=validation`

**Pass criteria**:
- Pipeline completes without errors
- Fitted masses within 10% of Julia code
- Plots look reasonable (no NaNs, smooth profiles)

---

### Regression Tests (Continuous)
Run on every change to catch regressions.

**Coverage:**
- [ ] Saved reference outputs (correlation functions, fitted parameters)
- [ ] Check new code produces identical results (within floating point precision)

**Test data**: Fixed test cases (checked into git)

**Command**: `pytest tests/test_regression.py --reference-dir tests/reference_outputs/`

**Pass criteria**: All outputs match reference (np.testing.assert_allclose with rtol=1e-6)

---

### Production Validation (One-time)
Run before finalizing thesis.

**Coverage:**
- [ ] All 10 redshift slices at full resolution (20k×20k)
- [ ] All galaxy samples (centrals, satellites, mass bins)
- [ ] Compare final results to Julia code extensively

**Test data**: Full dataset

**Command**: `snakemake --config scale=production --profile slurm --jobs 50`

**Pass criteria**:
- All expected outputs generated
- Masses match Julia within 10% (systematic difference OK if documented)
- No NaN/Inf in any output file
- Thesis figures reproducible from outputs

---

## SUCCESS CRITERIA

### Technical Success
✅ All 12 issues resolved (no exceptions, no compromises)
✅ 100% of Julia functionality ported to Python
✅ Test suite with > 80% code coverage, all passing
✅ Pipeline runs 10× faster than Julia code (TreeCorr optimization)
✅ Memory efficient (can run on 32GB laptop for debug scale)
✅ Reproducible (same inputs → same outputs within 0.1%)

### Scientific Success
✅ Correlation functions validated against TreeCorr (< 1% difference)
✅ Fitted halo masses match Julia code (< 10% difference, systematic documented)
✅ Error bars reasonable (jackknife properly accounts for spatial correlations)
✅ Thesis results reproducible by external user from scratch
✅ Code suitable for publication (clean, documented, tested)

### Engineering Success  
✅ Modular architecture (6 clean modules, single responsibility)
✅ No external C dependencies (pure Python + standard libraries)
✅ Explicit error handling (no silent failures)
✅ Comprehensive documentation (README + API docs + docstrings)
✅ Production pipeline (Snakemake, provenance tracking, HPC ready)
✅ Maintainable (future student can extend without pain)

---

## MODERN LIBRARY SUMMARY

| Purpose | Library | Version | Why This One? |
|---------|---------|---------|---------------|
| Arrays | NumPy | ≥1.24 | Standard, optimized, type hints |
| Scientific | SciPy | ≥1.10 | Curve fitting, interpolation, optimization |
| Astronomy | Astropy | ≥5.0 | Cosmology, FITS I/O, coordinates, units |
| Weak Lensing | **TreeCorr** | ≥4.3 | **CRITICAL**: Gold standard correlations, used by all surveys |
| Spatial | HEALPix | ≥1.16 | Hierarchical spherical pixelization, CMB standard |
| Workflow | Snakemake | ≥7.0 | Python-based, automatic parallelization, HPC integration |
| Testing | pytest | ≥7.0 | Standard, plugins, good CI/CD support |
| Plotting | Matplotlib | ≥3.6 | Publication quality, LaTeX support |
| Optional | LensTools | ≥1.1 | Complementary lensing utilities |

**Key insight**: Use TreeCorr for correlations (don't reinvent), use HEALPix for spatial decomposition (enables testing), use Snakemake for workflows (reproducibility).

---

## NOTES & CONSIDERATIONS

### What Gets Deleted
- [ ] All 22 Julia (.jl) files → archive to `backup/julia_original/`
- [ ] correlation_plots_sis_model.jl (25KB dead code)
- [ ] correlation_plots_sis_model_2.jl (25KB dead code)
- [ ] core.jl (libsoftlens.so dependency)
- [ ] 9 duplicate rebuild_*_maps.jl scripts

**Retained for reference:**
- `papers/` directory (PDF references)
- `Data/` directory (input data)
- Thesis PDF (`papers/laurant_flexion_HAGN_2024.pdf`)

### Git History
Archive original state before starting:
```bash
git checkout -b julia-original
git push origin julia-original  # Preserve for posterity

git checkout main
git checkout -b python-conversion
# ... implement plan
```

### Backward Compatibility
**Not a goal.** Clean break from Julia.

If comparison needed:
- Keep Julia environment in separate directory
- Run both pipelines in parallel for validation
- Document any systematic differences (> 1%)

### Dependencies Management
```bash
# requirements.txt
numpy>=1.24
scipy>=1.10
astropy>=5.0
treecorr>=4.3
healpy>=1.16
matplotlib>=3.6
pytest>=7.0
snakemake>=7.0

# Optional
lenstools>=1.1
hypothesis>=6.0
```

### Estimated Timeline
- **Minimum (expert Python dev)**: 6 weeks
- **Realistic (grad student, part-time)**: 3 months
- **Conservative (first time, debugging)**: 4 months

**Critical path**: TreeCorr validation (Issue #10) - must match within 1% or investigate discrepancy.

---

**END OF PLAN**

This plan is ready for implementation. All 12 issues approved with Option A (recommended solutions).

**Next step**: Begin Phase 1 (Foundation) - create package structure and port I/O functions.

**Questions before starting?** Review plan, clarify uncertainties, then proceed.

