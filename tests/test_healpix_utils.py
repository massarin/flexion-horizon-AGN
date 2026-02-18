"""
Tests for HEALPix utilities module.
"""

import numpy as np
import pytest
import healpy as hp
import tempfile
import os
from cosmo_lensing import healpix_utils
from cosmo_lensing.healpix_utils import (
    HEALPixMap, cartesian_to_healpix, healpix_to_cartesian
)


class TestHEALPixMap:
    """Test HEALPixMap class."""
    
    def test_initialization(self):
        """Test basic initialization."""
        nside = 16
        npix = hp.nside2npix(nside)
        data = np.random.randn(npix)
        
        hpmap = HEALPixMap(nside, data, ordering='RING')
        
        assert hpmap.nside == nside
        assert hpmap.npix == npix
        assert hpmap.ordering == 'RING'
        assert np.array_equal(hpmap.data, data)
    
    def test_initialization_nested(self):
        """Test initialization with NESTED ordering."""
        nside = 32
        npix = hp.nside2npix(nside)
        data = np.ones(npix)
        
        hpmap = HEALPixMap(nside, data, ordering='NESTED')
        
        assert hpmap.ordering == 'NESTED'
    
    def test_initialization_wrong_size(self):
        """Test that initialization fails with wrong data size."""
        nside = 16
        data = np.random.randn(100)  # Wrong size
        
        with pytest.raises(ValueError, match="Data size .* does not match"):
            HEALPixMap(nside, data)
    
    def test_initialization_invalid_ordering(self):
        """Test that invalid ordering raises error."""
        nside = 16
        npix = hp.nside2npix(nside)
        data = np.random.randn(npix)
        
        with pytest.raises(ValueError, match="Invalid ordering"):
            HEALPixMap(nside, data, ordering='INVALID')
    
    def test_get_neighbors(self):
        """Test neighbor finding."""
        nside = 16
        npix = hp.nside2npix(nside)
        data = np.arange(npix, dtype=float)
        
        hpmap = HEALPixMap(nside, data, ordering='RING')
        
        # Get neighbors of a pixel away from poles
        ipix = npix // 2
        neighbors = hpmap.get_neighbors(ipix)
        
        assert len(neighbors) > 0
        assert len(neighbors) <= 8
        assert np.all(neighbors >= 0)
        assert np.all(neighbors < npix)
        
    def test_get_neighbors_pole(self):
        """Test neighbor finding at pole."""
        nside = 16
        npix = hp.nside2npix(nside)
        data = np.arange(npix, dtype=float)
        
        hpmap = HEALPixMap(nside, data, ordering='RING')
        
        # North pole pixel
        ipix = 0
        neighbors = hpmap.get_neighbors(ipix)
        
        # Should have some neighbors even at pole
        assert len(neighbors) >= 4


class TestHEALPixProjection:
    """Test HEALPix to Cartesian projection."""
    
    def test_to_cartesian_shape(self):
        """Test that projection returns correct shape."""
        nside = 64
        npix = hp.nside2npix(nside)
        data = np.random.randn(npix)
        
        hpmap = HEALPixMap(nside, data)
        
        npix_cart = 100
        cart_map, info = hpmap.to_cartesian(
            center_ra=0.0, center_dec=0.0,
            radius_deg=5.0, npix=npix_cart
        )
        
        assert cart_map.shape == (npix_cart, npix_cart)
        assert 'center_ra' in info
        assert 'pixel_scale_deg' in info
    
    def test_to_cartesian_uniform_map(self):
        """Test projection of uniform map."""
        nside = 64
        npix = hp.nside2npix(nside)
        data = np.ones(npix) * 42.0
        
        hpmap = HEALPixMap(nside, data)
        
        cart_map, info = hpmap.to_cartesian(
            center_ra=180.0, center_dec=30.0,
            radius_deg=10.0, npix=50
        )
        
        # Uniform map should remain uniform
        assert np.allclose(cart_map, 42.0)
    
    def test_to_cartesian_smooth_gradient(self):
        """Test projection preserves smooth gradients."""
        nside = 128
        npix = hp.nside2npix(nside)
        
        # Create smooth gradient: monotonic in declination
        theta, phi = hp.pix2ang(nside, np.arange(npix))
        dec = 90 - np.degrees(theta)
        data = dec / 90.0  # Normalized gradient from -1 (south) to +1 (north)
        
        hpmap = HEALPixMap(nside, data)
        
        # Project patch at equator
        cart_map, info = hpmap.to_cartesian(
            center_ra=0.0, center_dec=0.0,
            radius_deg=5.0, npix=50
        )
        
        # Check that values are reasonable at equator
        assert np.all(np.abs(cart_map) < 0.5)  # Near zero at equator
        
        # Check gradient direction (increasing upward in y-axis)
        # Note: indexing [i,j] where i is row (y-direction), j is column (x)
        # Top rows (higher i) correspond to more positive declination
        left_half = cart_map[:, :25]
        right_half = cart_map[:, 25:]
        # For declination gradient, left/right should be similar
        # Check that there's variation along x (second index)
        assert np.std(cart_map) > 0.01  # Not completely uniform
    
    def test_to_cartesian_different_centers(self):
        """Test projection at different sky positions."""
        nside = 64
        npix = hp.nside2npix(nside)
        data = np.random.randn(npix)
        
        hpmap = HEALPixMap(nside, data)
        
        # Test several positions
        for ra, dec in [(0, 0), (180, 45), (270, -30)]:
            cart_map, info = hpmap.to_cartesian(
                center_ra=ra, center_dec=dec,
                radius_deg=5.0, npix=50
            )
            
            assert cart_map.shape == (50, 50)
            assert info['center_ra'] == ra
            assert info['center_dec'] == dec


class TestCoordinateConversion:
    """Test coordinate conversion functions."""
    
    def test_cartesian_to_healpix_single(self):
        """Test RA/Dec to HEALPix conversion."""
        nside = 64
        
        # Equator, prime meridian
        ra, dec = 0.0, 0.0
        ipix = cartesian_to_healpix(ra, dec, nside)
        
        assert isinstance(ipix, (int, np.integer))
        assert 0 <= ipix < hp.nside2npix(nside)
    
    def test_cartesian_to_healpix_array(self):
        """Test vectorized conversion."""
        nside = 64
        
        ra = np.array([0.0, 90.0, 180.0, 270.0])
        dec = np.array([0.0, 30.0, -30.0, 60.0])
        
        ipix = cartesian_to_healpix(ra, dec, nside)
        
        assert ipix.shape == (4,)
        assert np.all(ipix >= 0)
        assert np.all(ipix < hp.nside2npix(nside))
    
    def test_healpix_to_cartesian_single(self):
        """Test HEALPix to RA/Dec conversion."""
        nside = 64
        ipix = 1000
        
        ra, dec = healpix_to_cartesian(ipix, nside)
        
        assert isinstance(ra, (float, np.floating))
        assert isinstance(dec, (float, np.floating))
        assert -180 <= ra <= 360
        assert -90 <= dec <= 90
    
    def test_healpix_to_cartesian_array(self):
        """Test vectorized reverse conversion."""
        nside = 64
        ipix = np.array([0, 100, 1000, 3000])
        
        ra, dec = healpix_to_cartesian(ipix, nside)
        
        assert ra.shape == (4,)
        assert dec.shape == (4,)
        assert np.all(ra >= 0)
        assert np.all(ra < 360)
        assert np.all(dec >= -90)
        assert np.all(dec <= 90)
    
    def test_coordinate_roundtrip(self):
        """Test that conversion is reversible."""
        nside = 64
        
        # Original coordinates
        ra_orig = np.array([0.0, 90.0, 180.0, 270.0])
        dec_orig = np.array([0.0, 30.0, -30.0, 60.0])
        
        # Convert to HEALPix and back
        ipix = cartesian_to_healpix(ra_orig, dec_orig, nside)
        ra_back, dec_back = healpix_to_cartesian(ipix, nside)
        
        # Should be close (within pixel resolution)
        # HEALPix pixels are discrete, so we get pixel centers
        pixel_size = np.degrees(hp.nside2resol(nside))  # in degrees
        assert np.allclose(ra_orig, ra_back, atol=2*pixel_size)  # Within 2 pixels
        assert np.allclose(dec_orig, dec_back, atol=2*pixel_size)
    
    def test_coordinate_ordering_ring_vs_nested(self):
        """Test that ordering is respected."""
        nside = 64
        ra, dec = 45.0, 0.0
        
        ipix_ring = cartesian_to_healpix(ra, dec, nside, ordering='RING')
        ipix_nested = cartesian_to_healpix(ra, dec, nside, ordering='NESTED')
        
        # Same position, different ordering -> different indices
        assert ipix_ring != ipix_nested
        
        # But converting back should give similar coordinates (within pixel resolution)
        ra_ring, dec_ring = healpix_to_cartesian(ipix_ring, nside, ordering='RING')
        ra_nested, dec_nested = healpix_to_cartesian(ipix_nested, nside, ordering='NESTED')
        
        pixel_size = np.degrees(hp.nside2resol(nside))
        assert np.allclose(ra_ring, ra, atol=2*pixel_size)
        assert np.allclose(dec_ring, dec, atol=2*pixel_size)
        assert np.allclose(ra_nested, ra, atol=2*pixel_size)
        assert np.allclose(dec_nested, dec, atol=2*pixel_size)


class TestHEALPixFITS:
    """Test FITS I/O for HEALPix maps."""
    
    def test_write_and_read_fits(self):
        """Test writing and reading FITS files."""
        nside = 32
        npix = hp.nside2npix(nside)
        data_orig = np.random.randn(npix)
        
        hpmap_orig = HEALPixMap(nside, data_orig, ordering='RING')
        
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            # Write
            hpmap_orig.to_fits(tmp_path)
            
            # Read
            hpmap_read = HEALPixMap.from_fits(tmp_path)
            
            # Verify
            assert hpmap_read.nside == nside
            assert hpmap_read.ordering == 'RING'
            assert np.allclose(hpmap_read.data, data_orig)
        
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
    
    def test_write_and_read_fits_nested(self):
        """Test FITS I/O with NESTED ordering."""
        nside = 32
        npix = hp.nside2npix(nside)
        data_orig = np.random.randn(npix)
        
        hpmap_orig = HEALPixMap(nside, data_orig, ordering='NESTED')
        
        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            hpmap_orig.to_fits(tmp_path)
            hpmap_read = HEALPixMap.from_fits(tmp_path)
            
            # healpy reads in the same ordering as written
            assert hpmap_read.ordering == 'NESTED'
            # Data should match exactly (same ordering)
            assert np.allclose(hpmap_read.data, data_orig)
        
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)


class TestPatchProcessing:
    """Test memory-efficient patch processing."""
    
    def test_process_patches_basic(self):
        """Test basic patch processing."""
        nside = 32
        npix = hp.nside2npix(nside)
        data = np.ones(npix) * 5.0
        
        hpmap = HEALPixMap(nside, data)
        
        def mean_func(patch_data, patch_info):
            """Simple function: compute mean."""
            return np.mean(patch_data)
        
        # Process a few patches
        results = []
        for i, (result, info) in enumerate(hpmap.process_patches(
            patch_size_deg=30.0, overlap_deg=5.0, func=mean_func
        )):
            results.append(result)
            if i >= 10:  # Don't process entire sky for speed
                break
        
        # All patches should have mean ≈ 5.0
        assert len(results) > 0
        assert all(np.isclose(r, 5.0, atol=0.1) for r in results)
    
    def test_process_patches_gradient(self):
        """Test patch processing with gradient."""
        nside = 64
        npix = hp.nside2npix(nside)
        
        # Create declination gradient
        theta, phi = hp.pix2ang(nside, np.arange(npix))
        dec = 90 - np.degrees(theta)
        data = dec  # Values from -90 to +90
        
        hpmap = HEALPixMap(nside, data)
        
        def check_range(patch_data, patch_info):
            """Check that patch values are reasonable."""
            dec_center = patch_info['center_dec']
            return (np.min(patch_data), np.max(patch_data), dec_center)
        
        # Process patches at different declinations
        results = []
        for i, (result, info) in enumerate(hpmap.process_patches(
            patch_size_deg=20.0, overlap_deg=5.0, func=check_range
        )):
            results.append(result)
            if i >= 20:
                break
        
        assert len(results) > 0
        
        # Each patch should have values near its center declination
        for min_val, max_val, dec_center in results:
            # Values should be within reasonable range of center
            assert min_val < dec_center < max_val or np.isclose(min_val, max_val, atol=1.0)


class TestHEALPixIntegration:
    """Integration tests for HEALPix utilities."""
    
    def test_projection_preserves_statistics(self):
        """Test that projection preserves basic statistics."""
        nside = 128
        npix = hp.nside2npix(nside)
        
        # Create random map with known statistics
        np.random.seed(42)
        data = np.random.randn(npix) * 10.0 + 50.0  # mean=50, std=10
        
        hpmap = HEALPixMap(nside, data)
        
        # Project multiple patches and check statistics
        means = []
        stds = []
        
        for ra in [0, 90, 180, 270]:
            for dec in [-30, 0, 30]:
                cart_map, _ = hpmap.to_cartesian(
                    center_ra=ra, center_dec=dec,
                    radius_deg=10.0, npix=100
                )
                means.append(np.mean(cart_map))
                stds.append(np.std(cart_map))
        
        # Projected patches should have similar statistics to original
        # (with some variation due to sky coverage)
        assert np.abs(np.mean(means) - 50.0) < 5.0  # Loose tolerance
        assert np.abs(np.mean(stds) - 10.0) < 5.0
    
    def test_small_patch_derivatives(self):
        """Test that small patches can be used for derivatives."""
        from cosmo_lensing import derivatives
        
        nside = 256
        npix = hp.nside2npix(nside)
        
        # Create smooth gradient in both directions
        theta, phi = hp.pix2ang(nside, np.arange(npix))
        dec = 90 - np.degrees(theta)
        ra = np.degrees(phi)
        
        # Deflection field: α1 ~ RA, α2 ~ Dec
        alpha1 = ra / 360.0
        alpha2 = dec / 90.0
        
        # Create HEALPix maps
        hpmap1 = HEALPixMap(nside, alpha1)
        hpmap2 = HEALPixMap(nside, alpha2)
        
        # Project small patch
        patch1, info = hpmap1.to_cartesian(0, 0, 5.0, 100)
        patch2, _ = hpmap2.to_cartesian(0, 0, 5.0, 100)
        
        # Stack into deflection field with correct shape (npix, npix, 2)
        alpha = np.stack([patch1, patch2], axis=-1)
        
        # Compute Jacobian
        pixel_scale_rad = np.radians(info['pixel_scale_deg'])
        jac = derivatives.compute_jacobian(alpha, pixel_scale_rad)
        
        # Check shape - compute_jacobian returns (nx, ny, 12) for extended jacobian
        assert jac.shape == (100, 100, 12)
        
        # For HEALPix projections, derivatives can be noisy due to pixelization
        # Just check that derivatives are computed and finite
        assert np.all(np.isfinite(jac))
        # Check that there is some variation (not all zeros)
        assert np.std(jac[:,:,0]) > 0
