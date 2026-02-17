"""
End-to-end test of production pipeline scripts.

Creates synthetic deflection data, runs compute_observables.py and
compute_correlations.py scripts to validate the full workflow.
"""

import tempfile
from pathlib import Path
import subprocess
import numpy as np
from astropy.io import fits
from astropy.table import Table

import pytest


def create_test_deflection_binary(output_path: Path, size: int = 100):
    """Create a test deflection field binary file in Fortran format."""
    import struct
    
    # Generate simple deflection (point mass at center)
    y, x = np.mgrid[-size/2:size/2, -size/2:size/2]
    r = np.sqrt(x**2 + y**2) + 1e-6
    
    # Point mass deflection: α ∝ 1/r
    einstein_radius = 10.0
    alpha_x = einstein_radius * x / r**2
    alpha_y = einstein_radius * y / r**2
    
    # Write to binary (Fortran-style format matching io.py expectations)
    with open(output_path, 'wb') as f:
        # Record 1: Size (2 Int32: nx, ny)
        f.write(struct.pack('i', 8))  # record length
        f.write(struct.pack('ii', size, size))
        f.write(struct.pack('i', 8))  # record length
        
        # Record 2: Dummy values (3 Float32)
        f.write(struct.pack('i', 12))
        f.write(struct.pack('fff', 0.0, 0.0, 0.0))
        f.write(struct.pack('i', 12))
        
        # Record 3: Redshift (1 Float64)
        f.write(struct.pack('i', 8))
        f.write(struct.pack('d', 0.5))
        f.write(struct.pack('i', 8))
        
        # Record 4+: Deflection data in chunks
        # For simplicity, write entire array in one chunk
        ntot = size * size
        chunk_data_alpha1 = alpha_x.flatten().astype(np.float32)
        chunk_data_alpha2 = alpha_y.flatten().astype(np.float32)
        
        # Write α₁
        f.write(struct.pack('i', ntot * 4))
        chunk_data_alpha1.tofile(f)
        f.write(struct.pack('i', ntot * 4))
        
        # Write α₂
        f.write(struct.pack('i', ntot * 4))
        chunk_data_alpha2.tofile(f)
        f.write(struct.pack('i', ntot * 4))
    
    deflection = np.stack([alpha_x, alpha_y], axis=-1).astype(np.float32)
    return deflection


def create_test_galaxy_catalog(output_path: Path, n_galaxies: int = 50):
    """Create a test galaxy catalog FITS file."""
    # Random positions within a small field
    ra = np.random.uniform(-0.5, 0.5, n_galaxies)
    dec = np.random.uniform(-0.5, 0.5, n_galaxies)
    
    catalog = Table({
        'RA': ra,
        'DEC': dec,
        'MAG': np.random.uniform(20, 25, n_galaxies)
    })
    
    catalog.write(output_path, format='fits', overwrite=True)
    
    return catalog


@pytest.mark.slow
def test_compute_observables_script():
    """Test compute_observables.py script."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Create test data
        deflection_file = tmpdir / "test_deflection.bin"
        create_test_deflection_binary(deflection_file, size=100)
        
        # Run script
        result = subprocess.run([
            'python', 'scripts/compute_observables.py',
            str(deflection_file),
            '--outdir', str(tmpdir),
            '--observables', 'kappa', 'gamma1', 'gamma2',
            '--no-flexion',
            '--overwrite'
        ], capture_output=True, text=True)
        
        # Check exit code
        assert result.returncode == 0, f"Script failed: {result.stderr}"
        
        # Check outputs exist
        assert (tmpdir / "kappa.fits").exists()
        assert (tmpdir / "gamma1.fits").exists()
        assert (tmpdir / "gamma2.fits").exists()
        
        # Verify FITS content
        with fits.open(tmpdir / "kappa.fits") as hdul:
            kappa = hdul[0].data
            header = hdul[0].header
            
            assert kappa.shape == (100, 100)
            assert 'REDSHIFT' in header
            assert header['REDSHIFT'] == 0.5
            assert 'MAP' in header
            assert header['MAP'] == 'kappa'


@pytest.mark.slow
def test_compute_correlations_script():
    """Test compute_correlations.py script."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Create test shear maps (FITS files)
        size = 100
        gamma1 = np.random.normal(0, 0.01, (size, size)).astype(np.float32)
        gamma2 = np.random.normal(0, 0.01, (size, size)).astype(np.float32)
        kappa = np.random.normal(0, 0.01, (size, size)).astype(np.float32)
        
        # Create headers with WCS
        header = fits.Header()
        header['NAXIS'] = 2
        header['NAXIS1'] = size
        header['NAXIS2'] = size
        header['PIXSCALE'] = 1.0  # arcsec
        header['REDSHIFT'] = 0.5
        
        fits.writeto(tmpdir / "gamma1.fits", gamma1, header, overwrite=True)
        fits.writeto(tmpdir / "gamma2.fits", gamma2, header, overwrite=True)
        fits.writeto(tmpdir / "kappa.fits", kappa, header, overwrite=True)
        
        # Create test catalog
        catalog_file = tmpdir / "test_catalog.fits"
        create_test_galaxy_catalog(catalog_file, n_galaxies=20)
        
        # Run script
        result = subprocess.run([
            'python', 'scripts/compute_correlations.py',
            '--gamma1', str(tmpdir / "gamma1.fits"),
            '--gamma2', str(tmpdir / "gamma2.fits"),
            '--kappa', str(tmpdir / "kappa.fits"),
            '--catalog', str(catalog_file),
            '--outdir', str(tmpdir),
            '--rmin', '0.1',
            '--rmax', '5.0',
            '--nbins', '5',
            '--var-method', 'shot',
            '--sampling', '5',
            '--overwrite'
        ], capture_output=True, text=True)
        
        # Check exit code
        assert result.returncode == 0, f"Script failed: {result.stderr}"
        
        # Check output exists
        output_file = tmpdir / "tangential_shear.npz"
        assert output_file.exists()
        
        # Verify NPZ content
        data = np.load(output_file)
        assert 'r' in data
        assert 'xi' in data
        assert 'xi_err' in data
        assert 'npairs' in data
        
        assert len(data['r']) == 5
        assert np.all(np.isfinite(data['xi']))
        assert np.all(data['xi_err'] >= 0)


@pytest.mark.slow
def test_end_to_end_pipeline():
    """Test full pipeline: deflection → observables → correlations."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Create synthetic data
        deflection_file = tmpdir / "deflection.bin"
        create_test_deflection_binary(deflection_file, size=150)
        
        catalog_file = tmpdir / "catalog.fits"
        create_test_galaxy_catalog(catalog_file, n_galaxies=30)
        
        # Step 1: Compute observables
        result1 = subprocess.run([
            'python', 'scripts/compute_observables.py',
            str(deflection_file),
            '--outdir', str(tmpdir),
            '--observables', 'kappa', 'gamma1', 'gamma2',
            '--no-flexion'
        ], capture_output=True, text=True)
        
        assert result1.returncode == 0, f"Observable computation failed: {result1.stderr}"
        
        # Step 2: Compute correlations
        result2 = subprocess.run([
            'python', 'scripts/compute_correlations.py',
            '--gamma1', str(tmpdir / "gamma1.fits"),
            '--gamma2', str(tmpdir / "gamma2.fits"),
            '--kappa', str(tmpdir / "kappa.fits"),
            '--catalog', str(catalog_file),
            '--outdir', str(tmpdir),
            '--rmin', '0.2',
            '--rmax', '10.0',
            '--nbins', '8',
            '--var-method', 'shot',
            '--sampling', '10'
        ], capture_output=True, text=True)
        
        assert result2.returncode == 0, f"Correlation computation failed: {result2.stderr}"
        
        # Verify final output
        correlation_file = tmpdir / "tangential_shear.npz"
        assert correlation_file.exists()
        
        data = np.load(correlation_file)
        
        # Check that we recovered some signal
        # Point mass should give declining tangential shear
        assert np.sum(data['npairs']) > 0, "No lens-source pairs found"
        
        print(f"\nPipeline test successful!")
        print(f"  Observable maps: {tmpdir / 'kappa.fits'}")
        print(f"  Correlation: {correlation_file}")
        print(f"  Total pairs: {np.sum(data['npairs']):.0f}")


if __name__ == '__main__':
    # Run tests directly
    print("Testing compute_observables.py...")
    test_compute_observables_script()
    print("✓ Passed\n")
    
    print("Testing compute_correlations.py...")
    test_compute_correlations_script()
    print("✓ Passed\n")
    
    print("Testing end-to-end pipeline...")
    test_end_to_end_pipeline()
    print("✓ Passed\n")
    
    print("All pipeline tests passed!")
