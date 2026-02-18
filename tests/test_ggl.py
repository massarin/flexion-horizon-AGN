"""
Tests for the galaxy-galaxy lensing (GGL) module.
"""

import numpy as np
import pytest

from cosmo_lensing import ggl
from cosmo_lensing.profiles import NFWProfile, SISProfile


class TestComputeGGLCorrelation:
    """Tests for compute_ggl_correlation function."""
    
    def test_basic_shear_correlation(self):
        """Test basic shear correlation computation."""
        # Create simple test data
        n_lens = 10
        n_src = 100
        
        ra_lens = np.random.uniform(0, 10, n_lens)
        dec_lens = np.random.uniform(0, 10, n_lens)
        ra_src = np.random.uniform(0, 10, n_src)
        dec_src = np.random.uniform(0, 10, n_src)
        
        gamma1 = np.random.randn(n_src) * 0.01
        gamma2 = np.random.randn(n_src) * 0.01
        
        result = ggl.compute_ggl_correlation(
            ra_lens, dec_lens, ra_src, dec_src,
            gamma1, gamma2,
            min_sep=0.1, max_sep=10, nbins=5, sep_units="arcmin"
        )
        
        # Check output structure
        assert "r" in result
        assert "xi" in result
        assert "xim" in result
        assert "xi_err" in result
        assert "npairs" in result
        
        # Check shapes
        assert result["r"].shape == (5,)
        assert result["xi"].shape == (5,)
    
    def test_with_weights(self):
        """Test correlation with source weights."""
        n_lens = 5
        n_src = 50
        
        ra_lens = np.random.uniform(0, 5, n_lens)
        dec_lens = np.random.uniform(0, 5, n_lens)
        ra_src = np.random.uniform(0, 5, n_src)
        dec_src = np.random.uniform(0, 5, n_src)
        
        gamma1 = np.random.randn(n_src) * 0.01
        gamma2 = np.random.randn(n_src) * 0.01
        weights = np.random.uniform(0.5, 1.5, n_src)
        
        result = ggl.compute_ggl_correlation(
            ra_lens, dec_lens, ra_src, dec_src,
            gamma1, gamma2, weights=weights,
            nbins=3
        )
        
        assert "weight" in result
        assert np.all(result["weight"] >= 0)
    
    def test_convergence_correlation(self):
        """Test correlation with scalar convergence field."""
        n_lens = 5
        n_src = 50
        
        ra_lens = np.random.uniform(0, 5, n_lens)
        dec_lens = np.random.uniform(0, 5, n_lens)
        ra_src = np.random.uniform(0, 5, n_src)
        dec_src = np.random.uniform(0, 5, n_src)
        
        kappa = np.random.randn(n_src) * 0.01
        zeros = np.zeros_like(kappa)
        
        result = ggl.compute_ggl_correlation(
            ra_lens, dec_lens, ra_src, dec_src,
            kappa, zeros,
            nbins=3
        )
        
        assert result["r"].shape == (3,)


class TestComputeAllGGLCorrelations:
    """Tests for compute_all_ggl_correlations function."""
    
    def test_all_observables(self):
        """Test computing correlations for all observables."""
        n_lens = 5
        n_src = 50
        
        ra_lens = np.random.uniform(0, 5, n_lens)
        dec_lens = np.random.uniform(0, 5, n_lens)
        ra_src = np.random.uniform(0, 5, n_src)
        dec_src = np.random.uniform(0, 5, n_src)
        
        gamma1 = np.random.randn(n_src) * 0.01
        gamma2 = np.random.randn(n_src) * 0.01
        F1 = np.random.randn(n_src) * 0.001
        F2 = np.random.randn(n_src) * 0.001
        G1 = np.random.randn(n_src) * 0.001
        G2 = np.random.randn(n_src) * 0.001
        kappa = np.random.randn(n_src) * 0.01
        
        results = ggl.compute_all_ggl_correlations(
            ra_lens, dec_lens, ra_src, dec_src,
            gamma1, gamma2, F1, F2, G1, G2,
            kappa=kappa,
            nbins=3
        )
        
        # Check all observables computed
        assert "shear" in results
        assert "F" in results
        assert "G" in results
        assert "kappa" in results
        
        # Check structure of each result
        for name, result in results.items():
            assert "r" in result
            assert "xi" in result
            assert result["r"].shape == (3,)
    
    def test_without_convergence(self):
        """Test computing correlations without convergence."""
        n_lens = 5
        n_src = 50
        
        ra_lens = np.random.uniform(0, 5, n_lens)
        dec_lens = np.random.uniform(0, 5, n_lens)
        ra_src = np.random.uniform(0, 5, n_src)
        dec_src = np.random.uniform(0, 5, n_src)
        
        gamma1 = np.random.randn(n_src) * 0.01
        gamma2 = np.random.randn(n_src) * 0.01
        F1 = np.random.randn(n_src) * 0.001
        F2 = np.random.randn(n_src) * 0.001
        G1 = np.random.randn(n_src) * 0.001
        G2 = np.random.randn(n_src) * 0.001
        
        results = ggl.compute_all_ggl_correlations(
            ra_lens, dec_lens, ra_src, dec_src,
            gamma1, gamma2, F1, F2, G1, G2,
            kappa=None,
            nbins=3
        )
        
        # Check convergence not computed
        assert "shear" in results
        assert "F" in results
        assert "G" in results
        assert "kappa" not in results


class TestFitGGLProfile:
    """Tests for fit_ggl_profile function."""
    
    def test_nfw_shear_fit(self):
        """Test NFW fitting to synthetic shear profile."""
        # Create synthetic NFW profile
        ks = 0.1
        rs = 5.0  # arcmin
        
        profile = NFWProfile(ks, rs)
        
        r_arcmin = np.logspace(0, 1.5, 10)  # 1 to ~30 arcmin
        gamma_t_true = profile.shear_tangential(r_arcmin)
        
        # Add small noise
        noise = np.random.randn(len(r_arcmin)) * 0.001
        gamma_t = gamma_t_true + noise
        gamma_t_err = np.ones_like(gamma_t) * 0.001
        
        # Fit
        fit = ggl.fit_ggl_profile(
            r_arcmin, gamma_t, gamma_t_err,
            profile_type="nfw",
            observable="shear",
        )
        
        # Check fit results structure
        assert "ks" in fit
        assert "rs" in fit
        assert "chi2" in fit
        assert "redchi2" in fit
        
        # Check parameters recovered (within ~10% for noisy data)
        assert fit["ks"] / ks == pytest.approx(1.0, rel=0.2)
        assert fit["rs"] / rs == pytest.approx(1.0, rel=0.3)
    
    def test_sis_shear_fit(self):
        """Test SIS fitting to synthetic shear profile."""
        # Create synthetic SIS profile
        theta_E = 1.5  # arcsec
        
        profile = SISProfile(theta_E)
        
        r_arcsec = np.logspace(0, 1.5, 8)  # 1 to ~30 arcsec
        gamma_t_true = profile.shear_tangential(r_arcsec)
        
        # Add small noise
        noise = np.random.randn(len(r_arcsec)) * 0.001
        gamma_t = gamma_t_true + noise
        gamma_t_err = np.ones_like(gamma_t) * 0.001
        
        # Fit
        fit = ggl.fit_ggl_profile(
            r_arcsec, gamma_t, gamma_t_err,
            profile_type="sis",
            observable="shear",
        )
        
        # Check fit results structure
        assert "theta_E" in fit
        assert "chi2" in fit
        
        # Check parameter recovered (within ~5% for low noise)
        assert fit["theta_E"] / theta_E == pytest.approx(1.0, rel=0.1)
    
    def test_invalid_profile_type(self):
        """Test error on invalid profile type."""
        r = np.array([0.5, 1.0, 2.0])
        xi = np.array([0.01, 0.005, 0.002])
        xi_err = np.array([0.001, 0.001, 0.001])
        
        with pytest.raises(ValueError, match="Unknown profile_type"):
            ggl.fit_ggl_profile(
                r, xi, xi_err,
                profile_type="invalid"
            )


class TestFitAllGGLProfiles:
    """Tests for fit_all_ggl_profiles function."""
    
    def test_fit_all_profiles(self):
        """Test fitting profiles to all observables."""
        # Create synthetic data for shear
        ks = 0.1
        rs = 5.0  # arcmin
        
        profile = NFWProfile(ks, rs)
        
        # Create mock correlation results
        r_arcmin = np.logspace(0, 2, 8)  # 1 to 100 arcmin
        
        gamma_t = profile.shear_tangential(r_arcmin)
        
        # Add noise
        noise_scale = 0.001
        gamma_t += np.random.randn(len(r_arcmin)) * noise_scale
        
        # Create correlation results structure (just shear for now)
        correlations = {
            "shear": {
                "r": r_arcmin,
                "xi": gamma_t,
                "xi_err": np.ones_like(gamma_t) * noise_scale,
                "xim": np.zeros_like(gamma_t),
                "xim_err": np.ones_like(gamma_t) * noise_scale,
                "npairs": np.ones_like(gamma_t) * 100,
            },
        }
        
        # Fit all profiles
        fits = ggl.fit_all_ggl_profiles(
            correlations,
            profile_type="nfw",
        )
        
        # Check shear fitted
        assert "shear" in fits
        
        # Check fit parameters are reasonable
        fit = fits["shear"]
        assert "ks" in fit
        assert "rs" in fit
        # Parameters should be within factor of 2 of true values
        assert 0.5 < fit["ks"] / ks < 2.0
        assert 0.5 < fit["rs"] / rs < 2.0
    
    def test_skip_insufficient_data(self):
        """Test that observables with insufficient data are skipped."""
        # Create correlation with only 2 valid bins (need at least 3)
        correlations = {
            "shear": {
                "r": np.array([1.0, 2.0, 3.0]),
                "xi": np.array([0.01, 0.005, np.nan]),
                "xi_err": np.array([0.001, 0.001, 0.001]),
                "npairs": np.array([100, 50, 0]),
            }
        }
        
        fits = ggl.fit_all_ggl_profiles(
            correlations,
            profile_type="sis",
        )
        
        # Should skip this observable
        assert "shear" not in fits


class TestGGLIntegration:
    """Integration tests for full GGL workflow."""
    
    def test_full_workflow(self):
        """Test complete GGL workflow from correlation to fitting."""
        # Setup: Create synthetic lensed source catalog
        n_lens = 20
        n_src = 200
        ks = 0.1
        rs = 5.0  # arcmin
        
        # Lens positions (clustered)
        ra_lens = np.random.normal(5, 0.5, n_lens)
        dec_lens = np.random.uniform(4.5, 5.5, n_lens)
        
        # Source positions (uniform)
        ra_src = np.random.uniform(0, 10, n_src)
        dec_src = np.random.uniform(0, 10, n_src)
        
        # Compute distances and lensing signal
        profile = NFWProfile(ks, rs)
        
        # Simple approximation: distance from center
        center_ra, center_dec = ra_lens.mean(), dec_lens.mean()
        dx = (ra_src - center_ra) * np.cos(np.radians(dec_src)) * 60  # arcmin
        dy = (dec_src - center_dec) * 60  # arcmin
        r_arcmin = np.sqrt(dx**2 + dy**2)
        
        # Compute shear
        gamma_t = profile.shear_tangential(r_arcmin)
        phi = np.arctan2(dy, dx)
        gamma1 = -gamma_t * np.cos(2 * phi)
        gamma2 = -gamma_t * np.sin(2 * phi)
        
        # Add noise
        gamma1 += np.random.randn(n_src) * 0.01
        gamma2 += np.random.randn(n_src) * 0.01
        
        # Step 1: Compute correlation
        corr_result = ggl.compute_ggl_correlation(
            ra_lens, dec_lens, ra_src, dec_src,
            gamma1, gamma2,
            min_sep=1, max_sep=100, nbins=8, sep_units="arcmin"
        )
        
        # Check correlation computed
        assert "r" in corr_result
        assert "xi" in corr_result
        valid_bins = corr_result["npairs"] > 0
        assert valid_bins.sum() > 0
        
        # Step 2: Fit profile
        r_arcmin = corr_result["r"][valid_bins]
        
        fit = ggl.fit_ggl_profile(
            r_arcmin,
            corr_result["xi"][valid_bins],
            corr_result["xi_err"][valid_bins],
            profile_type="nfw",
            observable="shear",
        )
        
        # Check fit completed
        assert "ks" in fit
        assert "rs" in fit
        assert fit["chi2"] >= 0
        assert fit["dof"] > 0
        
        # With noisy small sample, parameters should be order-of-magnitude correct
        assert 0.01 < fit["ks"] < 1.0
        assert 1 < fit["rs"] < 50


class TestGGLScaling:
    """Scaling tests for GGL pipeline (marked slow)."""
    
    @pytest.mark.slow
    def test_large_catalog_performance(self):
        """Test performance with realistic catalog sizes."""
        import time
        
        # Realistic sizes
        n_lens = 5000
        n_src = 50000
        
        # Generate data
        ra_lens = np.random.uniform(0, 10, n_lens)
        dec_lens = np.random.uniform(-5, 5, n_lens)
        ra_src = np.random.uniform(0, 10, n_src)
        dec_src = np.random.uniform(-5, 5, n_src)
        
        gamma1 = np.random.randn(n_src) * 0.01
        gamma2 = np.random.randn(n_src) * 0.01
        
        # Time correlation
        start = time.time()
        result = ggl.compute_ggl_correlation(
            ra_lens, dec_lens, ra_src, dec_src,
            gamma1, gamma2,
            min_sep=1, max_sep=100, nbins=10, sep_units="arcmin"
        )
        elapsed = time.time() - start
        
        print(f"\n  Performance test:")
        print(f"    {n_lens} lenses, {n_src} sources")
        print(f"    Correlation time: {elapsed:.2f} seconds")
        print(f"    Pairs computed: {result['npairs'].sum():.0f}")
        
        # Should complete in reasonable time
        assert elapsed < 60, f"Too slow: {elapsed:.1f}s > 60s"
        
        # Should have reasonable number of pairs
        assert result['npairs'].sum() > 1000, "Too few pairs"
    
    @pytest.mark.slow
    def test_memory_efficiency(self):
        """Test memory usage doesn't grow with catalog size."""
        import psutil
        import os
        
        process = psutil.Process(os.getpid())
        mem_before = process.memory_info().rss / 1024**2  # MB
        
        # Process large catalog
        n_src = 100000
        ra_lens = np.random.uniform(0, 10, 1000)
        dec_lens = np.random.uniform(-5, 5, 1000)
        ra_src = np.random.uniform(0, 10, n_src)
        dec_src = np.random.uniform(-5, 5, n_src)
        gamma1 = np.random.randn(n_src) * 0.01
        gamma2 = np.random.randn(n_src) * 0.01
        
        result = ggl.compute_ggl_correlation(
            ra_lens, dec_lens, ra_src, dec_src,
            gamma1, gamma2, nbins=10
        )
        
        mem_after = process.memory_info().rss / 1024**2  # MB
        mem_increase = mem_after - mem_before
        
        print(f"\n  Memory test:")
        print(f"    Before: {mem_before:.1f} MB")
        print(f"    After: {mem_after:.1f} MB")
        print(f"    Increase: {mem_increase:.1f} MB")
        
        # Memory increase should be reasonable (<500 MB for 100k sources)
        assert mem_increase < 500, f"Memory leak? {mem_increase:.1f}MB increase"
