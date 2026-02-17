"""
Tests for tangential correlation functions.
"""

import numpy as np
import pytest

from cosmo_lensing.correlations import TangentialCorrelation


class TestTangentialCorrelation:
    """Test suite for TangentialCorrelation class."""
    
    def test_initialization(self):
        """Test basic initialization."""
        corr = TangentialCorrelation(
            rmin=1.0,
            rmax=100.0,
            nbins=10,
            sep_units='arcmin'
        )
        
        assert corr.rmin == 1.0
        assert corr.rmax == 100.0
        assert corr.nbins == 10
        assert corr.sep_units == 'arcmin'
    
    def test_compute_simple_shear(self):
        """Test computation with simple synthetic shear field."""
        corr = TangentialCorrelation(
            rmin=1.0,
            rmax=10.0,
            nbins=5,
            sep_units='arcmin',
            var_method='sample'  # Faster for tests
        )
        
        # Create synthetic data: single lens at origin
        ra_lens = np.array([0.0])
        dec_lens = np.array([0.0])
        
        # Sources in a ring around lens
        n_sources = 100
        theta = np.linspace(0, 2*np.pi, n_sources, endpoint=False)
        radius = 5.0  # arcmin
        ra_src = radius/60.0 * np.cos(theta)  # Convert to degrees
        dec_src = radius/60.0 * np.sin(theta)
        
        # Tangential shear pattern: γ_t = 0.1 (constant)
        # In (g1, g2) components: γ_t = -Re(γ * exp(-2iφ))
        # So g1 = -γ_t * cos(2φ), g2 = -γ_t * sin(2φ)
        gamma_t = 0.1
        g1 = -gamma_t * np.cos(2 * theta)
        g2 = -gamma_t * np.sin(2 * theta)
        
        result = corr.compute(ra_lens, dec_lens, ra_src, dec_src, g1, g2)
        
        # Verify output structure
        assert 'r' in result
        assert 'xi' in result
        assert 'xi_err' in result
        assert 'npairs' in result
        
        assert len(result['r']) == 5
        assert len(result['xi']) == 5
        
        # Check that at least one bin has pairs
        assert np.sum(result['npairs']) > 0
    
    def test_compute_nfw_like_profile(self):
        """Test with NFW-like tangential shear profile."""
        corr = TangentialCorrelation(
            rmin=1.0,
            rmax=50.0,
            nbins=8,
            sep_units='arcmin',
            var_method='sample'
        )
        
        # Single lens
        ra_lens = np.array([0.0])
        dec_lens = np.array([0.0])
        
        # Dense source distribution
        n_sources = 500
        np.random.seed(42)
        
        # Uniform in r (log space), uniform in θ
        log_r_min = np.log(1.0)
        log_r_max = np.log(50.0)
        log_r = np.random.uniform(log_r_min, log_r_max, n_sources)
        r = np.exp(log_r)  # arcmin
        theta = np.random.uniform(0, 2*np.pi, n_sources)
        
        ra_src = r/60.0 * np.cos(theta)  # degrees
        dec_src = r/60.0 * np.sin(theta)
        
        # NFW-like profile: γ_t ∝ 1/r at large r
        r_s = 10.0  # scale radius in arcmin
        gamma_t = 1.0 / (1 + r/r_s)  # Simplified NFW-like
        
        g1 = -gamma_t * np.cos(2 * theta)
        g2 = -gamma_t * np.sin(2 * theta)
        
        result = corr.compute(ra_lens, dec_lens, ra_src, dec_src, g1, g2)
        
        # Should recover declining profile (check magnitude)
        assert np.abs(result['xi'][0]) > np.abs(result['xi'][-1])
        
        # All bins should have pairs
        assert np.all(result['npairs'] > 0)
    
    def test_compute_with_weights(self):
        """Test computation with source weights."""
        corr = TangentialCorrelation(
            rmin=1.0,
            rmax=10.0,
            nbins=5,
            sep_units='arcmin',
            var_method='sample'
        )
        
        ra_lens = np.array([0.0])
        dec_lens = np.array([0.0])
        
        n_sources = 100
        theta = np.linspace(0, 2*np.pi, n_sources, endpoint=False)
        radius = 5.0
        ra_src = radius/60.0 * np.cos(theta)
        dec_src = radius/60.0 * np.sin(theta)
        
        gamma_t = 0.1
        g1 = -gamma_t * np.cos(2 * theta)
        g2 = -gamma_t * np.sin(2 * theta)
        
        # Random weights
        weights = np.random.uniform(0.5, 1.5, n_sources)
        
        result = corr.compute(
            ra_lens, dec_lens, ra_src, dec_src, g1, g2, weights=weights
        )
        
        assert 'weight' in result
        # Some bins might be empty, so check that at least one has weight
        assert np.sum(result['weight']) > 0
    
    def test_compute_zero_shear(self):
        """Test with zero shear field."""
        corr = TangentialCorrelation(
            rmin=1.0,
            rmax=10.0,
            nbins=5,
            sep_units='arcmin',
            var_method='sample'
        )
        
        ra_lens = np.array([0.0])
        dec_lens = np.array([0.0])
        
        n_sources = 100
        theta = np.linspace(0, 2*np.pi, n_sources, endpoint=False)
        radius = 5.0
        ra_src = radius/60.0 * np.cos(theta)
        dec_src = radius/60.0 * np.sin(theta)
        
        g1 = np.zeros(n_sources)
        g2 = np.zeros(n_sources)
        
        result = corr.compute(ra_lens, dec_lens, ra_src, dec_src, g1, g2)
        
        # Should measure near-zero correlation
        assert np.all(np.abs(result['xi']) < 0.01)
    
    def test_compute_multiple_lenses(self):
        """Test with multiple lenses."""
        corr = TangentialCorrelation(
            rmin=1.0,
            rmax=20.0,
            nbins=6,
            sep_units='arcmin',
            var_method='sample'
        )
        
        # Three lenses in a line
        ra_lens = np.array([0.0, 0.5, 1.0])
        dec_lens = np.array([0.0, 0.0, 0.0])
        
        # Many sources
        n_sources = 300
        np.random.seed(123)
        ra_src = np.random.uniform(-1, 2, n_sources)
        dec_src = np.random.uniform(-1, 1, n_sources)
        
        # Random shear
        g1 = np.random.normal(0, 0.1, n_sources)
        g2 = np.random.normal(0, 0.1, n_sources)
        
        result = corr.compute(ra_lens, dec_lens, ra_src, dec_src, g1, g2)
        
        # Should have some pairs
        assert np.sum(result['npairs']) > 30
    
    def test_bin_types(self):
        """Test different binning types."""
        # Logarithmic binning
        corr_log = TangentialCorrelation(
            rmin=1.0,
            rmax=100.0,
            nbins=5,
            bin_type='Log'
        )
        
        # Linear binning
        corr_lin = TangentialCorrelation(
            rmin=1.0,
            rmax=100.0,
            nbins=5,
            bin_type='Linear'
        )
        
        assert corr_log.bin_type == 'Log'
        assert corr_lin.bin_type == 'Linear'


class TestFlexionCorrelation:
    """Test suite for flexion correlations."""
    
    def test_compute_flexion(self):
        """Test basic flexion correlation computation."""
        corr = TangentialCorrelation(
            rmin=1.0,
            rmax=10.0,
            nbins=5,
            sep_units='arcmin',
            var_method='sample'
        )
        
        ra_lens = np.array([0.0])
        dec_lens = np.array([0.0])
        
        n_sources = 100
        theta = np.linspace(0, 2*np.pi, n_sources, endpoint=False)
        radius = 5.0
        ra_src = radius/60.0 * np.cos(theta)
        dec_src = radius/60.0 * np.sin(theta)
        
        # Simple flexion pattern
        F1 = 0.01 * np.cos(theta)
        F2 = 0.01 * np.sin(theta)
        G1 = 0.02 * np.cos(2*theta)
        G2 = 0.02 * np.sin(2*theta)
        
        F_result, G_result = corr.compute_flexion(
            ra_lens, dec_lens, ra_src, dec_src,
            F1, F2, G1, G2
        )
        
        # Both should have proper structure
        assert 'r' in F_result
        assert 'xi' in F_result
        assert 'r' in G_result
        assert 'xi' in G_result
        
        # Should have same radial bins
        assert len(F_result['r']) == len(G_result['r'])
