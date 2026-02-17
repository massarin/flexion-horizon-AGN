"""
Unit tests for cosmology module.

Tests cosmology calculations against astropy directly.
"""

import numpy as np
import pytest
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

from cosmo_lensing.cosmology import CosmologyCalculator


class TestCosmologyCalculator:
    """Test CosmologyCalculator class."""
    
    def test_initialization(self):
        """Test cosmology initialization with default parameters."""
        cosmo_calc = CosmologyCalculator()
        assert cosmo_calc.h == pytest.approx(0.6766, rel=1e-4)
        assert cosmo_calc.cosmo.Om0 == pytest.approx(0.3111, rel=1e-4)
    
    def test_initialization_custom(self):
        """Test cosmology initialization with custom parameters."""
        cosmo_calc = CosmologyCalculator(H0=70.0, Om0=0.3, Ob0=0.05)
        assert cosmo_calc.h == pytest.approx(0.7, rel=1e-4)
        assert cosmo_calc.cosmo.Om0 == pytest.approx(0.3, rel=1e-4)
    
    def test_angular_diameter_distance_single_z(self):
        """Test angular diameter distance at single redshift."""
        cosmo_calc = CosmologyCalculator()
        
        # Compare to astropy directly
        z = 0.5
        d_calc = cosmo_calc.angular_diameter_distance(z)
        d_astropy = cosmo_calc.cosmo.angular_diameter_distance(z).value
        
        assert d_calc == pytest.approx(d_astropy, rel=1e-6)
    
    def test_angular_diameter_distance_array(self):
        """Test angular diameter distance for array of redshifts."""
        cosmo_calc = CosmologyCalculator()
        
        z_array = np.array([0.1, 0.5, 1.0, 2.0])
        d_calc = cosmo_calc.angular_diameter_distance(z_array)
        d_astropy = cosmo_calc.cosmo.angular_diameter_distance(z_array).value
        
        np.testing.assert_allclose(d_calc, d_astropy, rtol=1e-6)
    
    def test_angular_diameter_distance_z1z2(self):
        """Test angular diameter distance between two redshifts."""
        cosmo_calc = CosmologyCalculator()
        
        zl, zs = 0.3, 0.8
        d_calc = cosmo_calc.angular_diameter_distance(zl, zs)
        d_astropy = cosmo_calc.cosmo.angular_diameter_distance_z1z2(zl, zs).value
        
        assert d_calc == pytest.approx(d_astropy, rel=1e-6)
    
    def test_comoving_distance(self):
        """Test comoving distance calculation."""
        cosmo_calc = CosmologyCalculator()
        
        z = 1.0
        d_calc = cosmo_calc.comoving_distance(z)
        d_astropy = cosmo_calc.cosmo.comoving_distance(z).value
        
        assert d_calc == pytest.approx(d_astropy, rel=1e-6)
    
    def test_critical_density(self):
        """Test critical density calculation."""
        cosmo_calc = CosmologyCalculator()
        
        z = 0.5
        rho_calc = cosmo_calc.critical_density(z)
        rho_astropy = cosmo_calc.cosmo.critical_density(z).to(u.Msun / u.Mpc**3).value
        
        assert rho_calc == pytest.approx(rho_astropy, rel=1e-6)
    
    def test_critical_density_z0(self):
        """Test critical density at z=0."""
        cosmo_calc = CosmologyCalculator()
        
        rho_0 = cosmo_calc.critical_density(0.0)
        # Planck 2018: rho_crit,0 ≈ 1.27e11 Msun/Mpc³ for h=0.6766
        assert rho_0 == pytest.approx(1.27e11, rel=0.01)
    
    def test_Omega_m(self):
        """Test matter density parameter evolution."""
        cosmo_calc = CosmologyCalculator()
        
        # At z=0, should equal Om0
        Om_0 = cosmo_calc.Omega_m(0.0)
        assert Om_0 == pytest.approx(cosmo_calc.cosmo.Om0, rel=1e-6)
        
        # At high z, should approach 1 (matter dominated)
        Om_high = cosmo_calc.Omega_m(10.0)
        assert Om_high > 0.95
    
    def test_sigma_crit(self):
        """Test critical surface density calculation."""
        cosmo_calc = CosmologyCalculator()
        
        zl, zs = 0.3, 0.8
        sigma_crit = cosmo_calc.sigma_crit(zl, zs)
        
        # Should be positive and reasonable value (1e14 - 1e16 Msun/Mpc²)
        assert sigma_crit > 0
        assert 1e14 < sigma_crit < 1e16
    
    def test_sigma_crit_source_behind_lens(self):
        """Test that sigma_crit raises error if source behind lens."""
        cosmo_calc = CosmologyCalculator()
        
        with pytest.raises(ValueError, match="greater than lens redshift"):
            cosmo_calc.sigma_crit(0.8, 0.3)  # zs < zl
    
    def test_sigma_crit_source_at_lens(self):
        """Test that sigma_crit raises error if source at lens redshift."""
        cosmo_calc = CosmologyCalculator()
        
        with pytest.raises(ValueError, match="greater than lens redshift"):
            cosmo_calc.sigma_crit(0.5, 0.5)  # zs = zl
    
    def test_distance_scaling(self):
        """Test distance scaling factors."""
        cosmo_calc = CosmologyCalculator()
        
        z = 0.5
        arcsec2kpc, sigma_crit_ref = cosmo_calc.distance_scaling(z)
        
        # arcsec2kpc should be reasonable (few kpc)
        assert 1.0 < arcsec2kpc < 10.0
        
        # sigma_crit_ref should be positive
        assert sigma_crit_ref > 0
    
    def test_distance_scaling_array(self):
        """Test distance scaling for array of redshifts."""
        cosmo_calc = CosmologyCalculator()
        
        z_array = np.array([0.2, 0.5, 1.0])
        arcsec2kpc, sigma_crit_ref = cosmo_calc.distance_scaling(z_array)
        
        assert len(arcsec2kpc) == len(z_array)
        assert len(sigma_crit_ref) == len(z_array)
        
        # arcsec2kpc should increase with redshift
        assert np.all(np.diff(arcsec2kpc) > 0)
    
    def test_comparison_with_literature(self):
        """Test cosmological distances against known literature values."""
        # Using Planck 2018 cosmology
        cosmo_calc = CosmologyCalculator()
        
        # At z=1, D_A ≈ 1700 Mpc (rough check)
        d_a_z1 = cosmo_calc.angular_diameter_distance(1.0)
        assert 1600 < d_a_z1 < 1800
        
        # At z=2, D_A ≈ 1750 Mpc (for Planck 2018)
        d_a_z2 = cosmo_calc.angular_diameter_distance(2.0)
        assert 1700 < d_a_z2 < 1800
