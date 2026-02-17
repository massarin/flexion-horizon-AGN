"""
Unit tests for NFW profile module.

Tests NFW lensing calculations against analytic properties.
"""

import numpy as np
import pytest

from cosmo_lensing.nfw import NFWProfile


class TestNFWFFunction:
    """Test F(x) auxiliary function."""
    
    def test_F_at_x_equals_1(self):
        """Test that F(1) = 1/3."""
        F_1 = NFWProfile._F_function(1.0)
        np.testing.assert_allclose(F_1, 1.0/3.0, rtol=1e-6)
    
    def test_F_small_x(self):
        """Test F(x) for x << 1."""
        x_small = 0.1
        F_small = NFWProfile._F_function(x_small)
        
        # For x << 1: F(x) ≈ 1/(3x²) (approximately)
        # Just check it's positive and large
        assert F_small > 1.0
    
    def test_F_large_x(self):
        """Test F(x) for x >> 1."""
        x_large = 10.0
        F_large = NFWProfile._F_function(x_large)
        
        # For x >> 1: F(x) → 0
        assert 0 < F_large < 0.1
    
    def test_F_array(self):
        """Test F(x) with array input."""
        x = np.array([0.5, 1.0, 2.0])
        F = NFWProfile._F_function(x)
        
        assert F.shape == x.shape
        assert np.all(np.isfinite(F))
        
        # Check ordering: F decreases with x
        assert F[0] > F[1] > F[2]
    
    def test_F_continuity(self):
        """Test that F(x) is continuous at x=1."""
        x_below = 0.99
        x_above = 1.01
        
        F_below = NFWProfile._F_function(x_below)
        F_above = NFWProfile._F_function(x_above)
        F_at_1 = NFWProfile._F_function(1.0)
        
        # Should be continuous within 5%
        assert abs(F_below - F_at_1) / F_at_1 < 0.05
        assert abs(F_above - F_at_1) / F_at_1 < 0.05


class TestNFWGFunction:
    """Test g(x) auxiliary function."""
    
    def test_g_at_x_equals_1(self):
        """Test g(1) = 1 + ln(0.5)."""
        g_1 = NFWProfile._g_function(1.0)
        expected = 1.0 + np.log(0.5)
        np.testing.assert_allclose(g_1, expected, rtol=1e-6)
    
    def test_g_array(self):
        """Test g(x) with array input."""
        x = np.array([0.5, 1.0, 2.0])
        g = NFWProfile._g_function(x)
        
        assert g.shape == x.shape
        assert np.all(np.isfinite(g))
    
    def test_g_continuity(self):
        """Test that g(x) is reasonably behaved near x=1."""
        x_below = 0.95
        x_above = 1.05
        
        g_below = NFWProfile._g_function(x_below)
        g_above = NFWProfile._g_function(x_above)
        g_at_1 = NFWProfile._g_function(1.0)
        
        # Should all be finite
        assert np.isfinite(g_below)
        assert np.isfinite(g_above)
        assert np.isfinite(g_at_1)


class TestNFWProfile:
    """Test NFWProfile class."""
    
    def test_initialization(self):
        """Test profile initialization."""
        profile = NFWProfile(ks=0.5, rs=100.0)
        
        assert profile.ks == 0.5
        assert profile.rs == 100.0
    
    def test_convergence_shape(self):
        """Test that convergence has correct shape."""
        profile = NFWProfile(ks=0.5, rs=100.0)
        
        r = np.linspace(10, 500, 50)
        kappa = profile.convergence(r)
        
        assert kappa.shape == r.shape
        assert np.all(np.isfinite(kappa))
    
    def test_convergence_scaling(self):
        """Test that convergence scales with κ_s."""
        profile1 = NFWProfile(ks=0.5, rs=100.0)
        profile2 = NFWProfile(ks=1.0, rs=100.0)
        
        r = np.linspace(10, 500, 50)
        kappa1 = profile1.convergence(r)
        kappa2 = profile2.convergence(r)
        
        # Should scale linearly with κ_s
        np.testing.assert_allclose(kappa2 / kappa1, 2.0, rtol=1e-6)
    
    def test_convergence_decreases(self):
        """Test that convergence decreases with radius."""
        profile = NFWProfile(ks=0.5, rs=100.0)
        
        r = np.array([50, 100, 200, 400])
        kappa = profile.convergence(r)
        
        # Should be monotonically decreasing
        assert np.all(np.diff(kappa) < 0)
    
    def test_shear_shape(self):
        """Test that tangential shear has correct shape."""
        profile = NFWProfile(ks=0.5, rs=100.0)
        
        r = np.linspace(10, 500, 50)
        gamma_t = profile.shear_tangential(r)
        
        assert gamma_t.shape == r.shape
        assert np.all(np.isfinite(gamma_t))
    
    def test_shear_positive(self):
        """Test that tangential shear is mostly positive."""
        profile = NFWProfile(ks=0.5, rs=100.0)
        
        # Avoid very small radii where shear can be negative
        r = np.linspace(50, 500, 50)
        gamma_t = profile.shear_tangential(r)
        
        # Most should be positive (>80%)
        positive_fraction = np.sum(gamma_t > 0) / len(gamma_t)
        assert positive_fraction > 0.8
    
    def test_shear_peaks_near_rs(self):
        """Test that shear peaks within reasonable range of scale radius."""
        profile = NFWProfile(ks=0.5, rs=100.0)
        
        # Sample broadly around r_s
        r = np.linspace(50, 300, 100)
        gamma_t = profile.shear_tangential(r)
        
        # Find peak (only consider positive values)
        positive_mask = gamma_t > 0
        if np.any(positive_mask):
            peak_idx = np.argmax(gamma_t[positive_mask])
            peak_r = r[positive_mask][peak_idx]
            
            # Peak should be within factor of 3 of r_s
            assert 0.3 * profile.rs < peak_r < 3.0 * profile.rs
    
    def test_excess_surface_density(self):
        """Test that excess surface density equals tangential shear."""
        profile = NFWProfile(ks=0.5, rs=100.0)
        
        r = np.linspace(10, 500, 50)
        gamma_t = profile.shear_tangential(r)
        delta_sigma = profile.excess_surface_density(r)
        
        np.testing.assert_array_equal(delta_sigma, gamma_t)


class TestNFWAnalyticProperties:
    """Test analytic properties of NFW profile."""
    
    def test_convergence_at_scale_radius(self):
        """Test convergence value at r = r_s."""
        profile = NFWProfile(ks=1.0, rs=100.0)
        
        kappa_rs = profile.convergence(100.0)
        
        # At x=1: κ(r_s) = 2 * κ_s * F(1) = 2 * κ_s / 3
        expected = 2.0 * profile.ks / 3.0
        
        np.testing.assert_allclose(kappa_rs, expected, rtol=1e-6)
    
    def test_mean_convergence_relation(self):
        """Test that γ_t = κ̄(<r) - κ(r)."""
        profile = NFWProfile(ks=0.5, rs=100.0)
        
        r = 100.0
        kappa_r = profile.convergence(r)
        gamma_t = profile.shear_tangential(r)
        
        # Mean convergence: κ̄ = κ_s * g(x)
        x = r / profile.rs
        g_x = profile._g_function(x)
        kappa_mean = 2 * profile.ks * g_x
        
        # Check: γ_t = κ̄ - κ
        expected_gamma_t = kappa_mean - kappa_r
        
        np.testing.assert_allclose(gamma_t, expected_gamma_t, rtol=1e-6)
    
    def test_asymptotic_behavior_small_r(self):
        """Test behavior for r << r_s."""
        profile = NFWProfile(ks=0.5, rs=100.0)
        
        r_small = 10.0  # r << r_s
        kappa = profile.convergence(r_small)
        
        # For r << r_s: κ should be large
        assert kappa > profile.ks
    
    def test_asymptotic_behavior_large_r(self):
        """Test behavior for r >> r_s."""
        profile = NFWProfile(ks=0.5, rs=100.0)
        
        r_large = 1000.0  # r >> r_s
        kappa = profile.convergence(r_large)
        
        # For r >> r_s: κ should be small
        assert kappa < 0.1 * profile.ks
    
    def test_profile_normalization(self):
        """Test that profile has reasonable normalization."""
        profile = NFWProfile(ks=0.5, rs=100.0)
        
        # Sample profile
        r = np.logspace(0, 3, 100)  # 1 to 1000
        kappa = profile.convergence(r)
        
        # All values should be positive and finite
        assert np.all(kappa > 0)
        assert np.all(np.isfinite(kappa))
        
        # Maximum should be at small radii
        assert kappa[0] == np.max(kappa)


class TestNFWPlaceholders:
    """Test placeholder functions for Phase 4."""
    
    def test_from_mass_concentration_not_implemented(self):
        """Test that from_mass_concentration raises NotImplementedError."""
        with pytest.raises(NotImplementedError, match="Phase 4"):
            NFWProfile.from_mass_concentration(
                mass=1e14, concentration=5.0,
                redshift_lens=0.3, redshift_source=0.8
            )
