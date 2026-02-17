"""
Input/output operations for weak lensing pipeline.

Handles:
- Reading Fortran binary deflection field maps
- Writing FITS maps with WCS headers
- Loading galaxy catalogs
- Explicit error handling and validation
"""

import logging
import struct
from pathlib import Path
from typing import Tuple, Optional

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

logger = logging.getLogger(__name__)


class DeflectionFieldError(Exception):
    """Raised when deflection field data is invalid or corrupted."""
    pass


def read_deflection_field(
    filename: str,
    chunksize: int = 1_000_000,
    info_only: bool = False
) -> Tuple[Optional[np.ndarray], Optional[float]]:
    """
    Read a binary deflection field map plus corresponding source redshift.
    
    The file format is Fortran unformatted binary with the following structure:
    - Record 1: size (2 Int32 values)
    - Record 2: dummy (3 Float32 values)
    - Record 3: source redshift (1 Float64)
    - Records 4+: deflection field data in chunks (Float32)
    
    Args:
        filename: Path to binary deflection field file
        chunksize: Size of chunks for reading (matches Fortran record size)
        info_only: If True, only read header (size and redshift)
    
    Returns:
        Tuple of (deflection_field, redshift) where:
        - deflection_field: (nx, ny, 2) array of α₁ and α₂ components
        - redshift: source redshift value
        
    Raises:
        FileNotFoundError: If file doesn't exist
        DeflectionFieldError: If file is corrupted or has invalid dimensions
    """
    filepath = Path(filename)
    
    if not filepath.exists():
        logger.error(f"Deflection field file not found: {filename}")
        raise FileNotFoundError(f"Cannot find file: {filename}")
    
    # Check file size (minimum viable file should be > 1KB)
    file_size = filepath.stat().st_size
    if file_size < 1000:
        logger.error(f"File {filename} is too small ({file_size} bytes), likely corrupted")
        raise DeflectionFieldError(
            f"File {filename} is too small ({file_size} bytes), likely corrupted"
        )
    
    try:
        with open(filepath, 'rb') as f:
            # Read Fortran record: 4-byte length prefix
            rec_len = struct.unpack('i', f.read(4))[0]
            
            # Read size: 2 Int32 values
            size_bytes = f.read(8)
            nx, ny = struct.unpack('ii', size_bytes)
            
            # Read Fortran record: 4-byte length suffix
            f.read(4)
            
            # Validate dimensions
            if nx <= 0 or ny <= 0:
                raise DeflectionFieldError(
                    f"Invalid dimensions in {filename}: nx={nx}, ny={ny}"
                )
            
            logger.info(f"Reading deflection field: {filename}, size=({nx}, {ny})")
            
            # Read dummy values (3 Float32)
            rec_len = struct.unpack('i', f.read(4))[0]
            f.read(12)  # Skip 3 Float32 values
            f.read(4)
            
            # Read redshift (Float64)
            rec_len = struct.unpack('i', f.read(4))[0]
            zs = struct.unpack('d', f.read(8))[0]
            f.read(4)
            
            if info_only:
                logger.debug(f"Info only: size=({nx}, {ny}), z={zs}")
                return None, zs
            
            # Allocate output array: (nx, ny, 2) for α₁ and α₂
            alpha = np.empty((nx, ny, 2), dtype=np.float32)
            
            # Total number of pixels
            ntot = nx * ny
            
            # Read deflection data in chunks
            i1 = 0
            while i1 < ntot:
                i2 = min(i1 + chunksize, ntot)
                nvalues = i2 - i1
                
                # Read α₁ component
                rec_len = struct.unpack('i', f.read(4))[0]
                expected_bytes = nvalues * 4
                
                if rec_len != expected_bytes:
                    logger.warning(
                        f"Record length mismatch at position {i1}: "
                        f"expected {expected_bytes}, got {rec_len}"
                    )
                
                a1_data = np.frombuffer(f.read(expected_bytes), dtype=np.float32)
                f.read(4)
                
                # Read α₂ component
                rec_len = struct.unpack('i', f.read(4))[0]
                a2_data = np.frombuffer(f.read(expected_bytes), dtype=np.float32)
                f.read(4)
                
                # Fill array (Fortran column-major order)
                alpha.ravel()[i1:i2] = a1_data
                alpha.ravel()[ntot + i1:ntot + i2] = a2_data
                
                i1 += chunksize
            
            # Validate data (check for NaN/Inf)
            if np.any(~np.isfinite(alpha)):
                n_bad = np.sum(~np.isfinite(alpha))
                logger.error(
                    f"Deflection field contains {n_bad} NaN/Inf values"
                )
                raise DeflectionFieldError(
                    f"Deflection field contains {n_bad} NaN/Inf values in {filename}"
                )
            
            logger.info(
                f"Successfully read deflection field: "
                f"{filename}, z={zs:.3f}, "
                f"min={alpha.min():.3e}, max={alpha.max():.3e}"
            )
            
            return alpha, zs
            
    except struct.error as e:
        logger.error(f"Binary format error reading {filename}: {e}")
        raise DeflectionFieldError(f"Binary format error in {filename}: {e}")
    except Exception as e:
        logger.error(f"Unexpected error reading {filename}: {e}")
        raise DeflectionFieldError(f"Error reading {filename}: {e}")


def write_fits_map(
    data: np.ndarray,
    filename: str,
    redshift: float,
    observable: str = "",
    pixel_scale: float = 1.0,
    overwrite: bool = True
) -> None:
    """
    Write a 2D map to FITS format with WCS header.
    
    Args:
        data: 2D array to write
        filename: Output FITS filename
        redshift: Source redshift
        observable: Name of observable (e.g., 'kappa', 'gamma1')
        pixel_scale: Pixel scale in arcsec
        overwrite: Whether to overwrite existing file
        
    Raises:
        ValueError: If data is not 2D
    """
    if data.ndim != 2:
        raise ValueError(f"Data must be 2D, got shape {data.shape}")
    
    # Create WCS header
    w = WCS(naxis=2)
    w.wcs.crpix = [data.shape[1] / 2, data.shape[0] / 2]  # Reference pixel
    w.wcs.cdelt = [pixel_scale / 3600, pixel_scale / 3600]  # Pixel scale in degrees
    w.wcs.crval = [0.0, 0.0]  # Reference coordinates
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    
    # Create FITS header
    header = w.to_header()
    header['REDSHIFT'] = (redshift, 'Source redshift')
    if observable:
        header['OBSRVBLE'] = (observable, 'Lensing observable')
    header['PIXSCALE'] = (pixel_scale, 'Pixel scale (arcsec)')
    header['BUNIT'] = ('dimensionless', 'Physical units')
    
    # Create HDU and write
    hdu = fits.PrimaryHDU(data=data.astype(np.float32), header=header)
    hdu.writeto(filename, overwrite=overwrite)
    
    logger.info(f"Wrote FITS map: {filename}, observable={observable}, z={redshift:.3f}")


def load_galaxy_catalog(filename: str) -> fits.fitsrec.FITS_rec:
    """
    Load galaxy catalog from FITS file.
    
    Args:
        filename: Path to FITS catalog
        
    Returns:
        FITS record array with catalog data
        
    Raises:
        FileNotFoundError: If file doesn't exist
    """
    filepath = Path(filename)
    
    if not filepath.exists():
        logger.error(f"Catalog file not found: {filename}")
        raise FileNotFoundError(f"Cannot find catalog: {filename}")
    
    with fits.open(filepath) as hdul:
        catalog = hdul[1].data
        
    logger.info(f"Loaded catalog: {filename}, {len(catalog)} entries")
    return catalog
