"""
Tests for new IO loading functions.
"""

import pytest
import numpy as np
from pathlib import Path
from astropy.io import fits
from cosmo_lensing import io
import tempfile


class TestLoadLensingMap:
    """Tests for load_lensing_map function."""

    def test_load_valid_map(self, tmp_path):
        """Test loading a valid lensing map."""
        # Create a simple test FITS file
        filename = tmp_path / "test_map.fits"

        # Create HDUs
        primary = fits.PrimaryHDU()
        gamma1_data = np.random.randn(100, 100).astype(np.float32)
        gamma2_data = np.random.randn(100, 100).astype(np.float32)

        gamma1_hdu = fits.ImageHDU(gamma1_data)
        gamma1_hdu.header["MAP"] = "gamma1"
        gamma1_hdu.header["REDSHIFT"] = 1.0
        gamma1_hdu.header["CDELT1"] = -0.000138889

        gamma2_hdu = fits.ImageHDU(gamma2_data)
        gamma2_hdu.header["MAP"] = "gamma2"

        hdul = fits.HDUList([primary, gamma1_hdu, gamma2_hdu])
        hdul.writeto(filename)

        # Load the map
        data, metadata = io.load_lensing_map(str(filename))

        # Verify structure
        assert "gamma1" in data
        assert "gamma2" in data
        assert data["gamma1"].shape == (100, 100)
        assert metadata["redshift"] == 1.0
        assert metadata["pixscale"] is not None

    def test_load_nonexistent_file(self):
        """Test that loading nonexistent file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            io.load_lensing_map("nonexistent.fits")

    def test_load_empty_fits(self, tmp_path):
        """Test loading FITS with no data extensions."""
        filename = tmp_path / "empty.fits"
        primary = fits.PrimaryHDU()
        hdul = fits.HDUList([primary])
        hdul.writeto(filename)

        with pytest.raises(IOError, match="No data extensions"):
            io.load_lensing_map(str(filename))


class TestLoadGalaxyCatalog:
    """Tests for load_galaxy_catalog function."""

    def test_load_valid_catalog(self, tmp_path):
        """Test loading a valid galaxy catalog."""
        filename = tmp_path / "test_catalog.fits"

        # Create test catalog
        n_gal = 1000
        col1 = fits.Column(name="RA_IMG", format="D", array=np.random.rand(n_gal) * 360)
        col2 = fits.Column(name="DEC_IMG", format="D", array=np.random.rand(n_gal) * 180 - 90)
        col3 = fits.Column(name="z_true", format="D", array=np.random.rand(n_gal) * 3)

        hdu = fits.BinTableHDU.from_columns([col1, col2, col3])
        primary = fits.PrimaryHDU()
        hdul = fits.HDUList([primary, hdu])
        hdul.writeto(filename)

        # Load catalog
        catalog = io.load_galaxy_catalog(str(filename))

        # Verify structure
        assert "ra" in catalog
        assert "dec" in catalog
        assert "z" in catalog
        assert len(catalog["ra"]) == n_gal

    def test_catalog_with_z_cuts(self, tmp_path):
        """Test redshift filtering."""
        filename = tmp_path / "test_catalog.fits"

        n_gal = 1000
        z_values = np.linspace(0, 3, n_gal)
        col1 = fits.Column(name="RA_IMG", format="D", array=np.random.rand(n_gal) * 360)
        col2 = fits.Column(name="DEC_IMG", format="D", array=np.random.rand(n_gal) * 180 - 90)
        col3 = fits.Column(name="z_true", format="D", array=z_values)

        hdu = fits.BinTableHDU.from_columns([col1, col2, col3])
        primary = fits.PrimaryHDU()
        hdul = fits.HDUList([primary, hdu])
        hdul.writeto(filename)

        # Load with redshift cuts
        catalog = io.load_galaxy_catalog(str(filename), z_min=1.0, z_max=2.0)

        # Verify filtering
        assert len(catalog["ra"]) < n_gal
        assert np.all(catalog["z"] >= 1.0)
        assert np.all(catalog["z"] <= 2.0)

    def test_catalog_subsampling(self, tmp_path):
        """Test random subsampling."""
        filename = tmp_path / "test_catalog.fits"

        n_gal = 1000
        col1 = fits.Column(name="RA_IMG", format="D", array=np.random.rand(n_gal) * 360)
        col2 = fits.Column(name="DEC_IMG", format="D", array=np.random.rand(n_gal) * 180 - 90)
        col3 = fits.Column(name="z_true", format="D", array=np.random.rand(n_gal) * 3)

        hdu = fits.BinTableHDU.from_columns([col1, col2, col3])
        primary = fits.PrimaryHDU()
        hdul = fits.HDUList([primary, hdu])
        hdul.writeto(filename)

        # Load with subsampling
        catalog = io.load_galaxy_catalog(str(filename), subsample=100)

        # Verify subsampling
        assert len(catalog["ra"]) == 100

    def test_catalog_missing_columns(self, tmp_path):
        """Test that missing required columns raises ValueError."""
        filename = tmp_path / "bad_catalog.fits"

        # Create catalog without required columns
        col1 = fits.Column(name="OTHER_COL", format="D", array=np.random.rand(100))

        hdu = fits.BinTableHDU.from_columns([col1])
        primary = fits.PrimaryHDU()
        hdul = fits.HDUList([primary, hdu])
        hdul.writeto(filename)

        with pytest.raises(ValueError, match="Required column"):
            io.load_galaxy_catalog(str(filename))

    def test_catalog_nonexistent_file(self):
        """Test that loading nonexistent file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            io.load_galaxy_catalog("nonexistent.fits")
