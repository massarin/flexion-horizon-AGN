#!/usr/bin/env python3
"""
Compute lensing observables from deflection fields.

This script reads deflection field binaries, computes the Jacobian matrix
and second derivatives, extracts all lensing observables (convergence, shear,
flexion, rotation), and writes them to FITS files.

Usage:
    python compute_observables.py deflection_50.bin --outdir results/
    python compute_observables.py --field-id 50 --datadir /path/to/data
    python compute_observables.py --config workflow/config.yaml

Outputs:
    - kappa_050.fits: Convergence map
    - gamma1_050.fits: Shear component 1
    - gamma2_050.fits: Shear component 2
    - F1_050.fits: First flexion component 1
    - F2_050.fits: First flexion component 2
    - G1_050.fits: Second flexion component 1
    - G2_050.fits: Second flexion component 2
    - rotation_050.fits: Rotation field (diagnostic)
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional, Dict, Any
import time

import numpy as np
import yaml
from astropy.io import fits

# Add parent directory to path for package import
sys.path.insert(0, str(Path(__file__).parent.parent))

from cosmo_lensing import io as lensing_io
from cosmo_lensing.derivatives import compute_jacobian
from cosmo_lensing.observables import (
    convergence, shear_1, shear_2,
    flexion_F1, flexion_F2,
    flexion_G1, flexion_G2,
    rotation
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def create_wcs_header(
    redshift: float,
    pixel_scale: float,
    shape: tuple,
    observable: str
) -> fits.Header:
    """
    Create FITS header with WCS keywords.
    
    Args:
        redshift: Source plane redshift
        pixel_scale: Pixel scale in arcsec
        shape: Image shape (ny, nx)
        observable: Observable type ('kappa', 'gamma1', etc.)
    
    Returns:
        FITS header with WCS and metadata
    """
    header = fits.Header()
    
    # WCS keywords (simplified - assumes tangent plane projection)
    header['NAXIS'] = 2
    header['NAXIS1'] = shape[1]
    header['NAXIS2'] = shape[0]
    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    header['CRPIX1'] = shape[1] / 2.0
    header['CRPIX2'] = shape[0] / 2.0
    header['CRVAL1'] = 0.0  # RA at reference pixel
    header['CRVAL2'] = 0.0  # Dec at reference pixel
    header['CDELT1'] = -pixel_scale / 3600.0  # degrees
    header['CDELT2'] = pixel_scale / 3600.0
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    
    # Metadata
    header['REDSHIFT'] = (redshift, 'Source plane redshift')
    header['PIXSCALE'] = (pixel_scale, 'Pixel scale [arcsec]')
    header['MAP'] = (observable, 'Observable type')
    header['BUNIT'] = 'dimensionless'
    
    return header


def compute_all_observables(
    deflection: np.ndarray,
    pixel_scale: float,
    compute_flexion: bool = True
) -> Dict[str, np.ndarray]:
    """
    Compute all lensing observables from deflection field.
    
    Args:
        deflection: Deflection field (ny, nx, 2) in arcsec
        pixel_scale: Pixel scale in arcsec
        compute_flexion: Whether to compute flexion (requires second derivatives)
    
    Returns:
        Dictionary with observable arrays:
        - 'kappa': Convergence
        - 'gamma1', 'gamma2': Shear components
        - 'F1', 'F2': First flexion components (if compute_flexion=True)
        - 'G1', 'G2': Second flexion components (if compute_flexion=True)
        - 'rotation': Rotation field (diagnostic)
    """
    logger.info(f"Computing Jacobian for {deflection.shape[:2]} field...")
    start = time.time()
    
    # Compute Jacobian and second derivatives
    jac = compute_jacobian(deflection, pixel_scale, compute_second_derivatives=compute_flexion)
    
    elapsed = time.time() - start
    logger.info(f"Jacobian computed in {elapsed:.2f}s")
    
    # Extract first-order observables
    logger.info("Computing convergence and shear...")
    observables = {
        'kappa': convergence(jac),
        'gamma1': shear_1(jac),
        'gamma2': shear_2(jac),
        'rotation': rotation(jac)
    }
    
    # Extract flexion if requested
    if compute_flexion:
        logger.info("Computing first and second flexion...")
        observables['F1'] = flexion_F1(jac)
        observables['F2'] = flexion_F2(jac)
        observables['G1'] = flexion_G1(jac)
        observables['G2'] = flexion_G2(jac)
    
    # Report statistics
    for name, arr in observables.items():
        logger.info(f"{name}: mean={np.mean(arr):.6f}, std={np.std(arr):.6f}, "
                   f"min={np.min(arr):.6f}, max={np.max(arr):.6f}")
    
    return observables


def save_observable(
    observable: np.ndarray,
    output_path: Path,
    redshift: float,
    pixel_scale: float,
    observable_name: str
):
    """Save single observable to FITS file."""
    header = create_wcs_header(redshift, pixel_scale, observable.shape, observable_name)
    
    hdu = fits.PrimaryHDU(data=observable.astype(np.float32), header=header)
    hdu.writeto(output_path, overwrite=True)
    
    logger.info(f"Saved {observable_name} to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Compute lensing observables from deflection fields',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input specification
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('deflection_file', nargs='?', type=Path,
                            help='Path to deflection field binary file')
    input_group.add_argument('--field-id', type=int,
                            help='Field ID number (requires --datadir)')
    input_group.add_argument('--config', type=Path,
                            help='Path to YAML config file')
    
    # Output control
    parser.add_argument('--outdir', type=Path, default=Path('results'),
                       help='Output directory for FITS files')
    parser.add_argument('--datadir', type=Path,
                       help='Directory containing deflection binaries (for --field-id)')
    parser.add_argument('--prefix', type=str, default='',
                       help='Prefix for output filenames')
    
    # Processing options
    parser.add_argument('--no-flexion', action='store_true',
                       help='Skip flexion computation (faster)')
    parser.add_argument('--observables', nargs='+',
                       default=['kappa', 'gamma1', 'gamma2'],
                       help='Observables to compute')
    
    # Advanced options
    parser.add_argument('--pixel-scale', type=float,
                       help='Override pixel scale (arcsec, read from file if not specified)')
    parser.add_argument('--overwrite', action='store_true',
                       help='Overwrite existing output files')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Determine input file
    if args.config:
        logger.info(f"Loading configuration from {args.config}")
        with open(args.config, 'r') as f:
            config = yaml.safe_load(f)
        # TODO: Parse config to get deflection file path
        raise NotImplementedError("Config file support coming soon")
    
    elif args.field_id is not None:
        if not args.datadir:
            parser.error("--field-id requires --datadir")
        deflection_file = args.datadir / f"deflection_{args.field_id:03d}.bin"
    else:
        deflection_file = args.deflection_file
    
    # Check input exists
    if not deflection_file.exists():
        logger.error(f"Deflection file not found: {deflection_file}")
        sys.exit(1)
    
    # Create output directory
    args.outdir.mkdir(parents=True, exist_ok=True)
    
    # Read deflection field
    logger.info(f"Reading deflection field from {deflection_file}")
    start = time.time()
    
    try:
        deflection, redshift = lensing_io.read_deflection_field(deflection_file)
        pixel_scale = args.pixel_scale if args.pixel_scale else 1.0  # Default pixel scale
        
        elapsed = time.time() - start
        logger.info(f"Read deflection field in {elapsed:.2f}s")
        logger.info(f"Shape: {deflection.shape}, redshift: {redshift}, "
                   f"pixel_scale: {pixel_scale} arcsec")
    
    except Exception as e:
        logger.error(f"Failed to read deflection field: {e}")
        sys.exit(1)
    
    # Compute observables
    compute_flexion = not args.no_flexion
    observables = compute_all_observables(deflection, pixel_scale, compute_flexion)
    
    # Save requested observables
    field_id = args.field_id if args.field_id else "unknown"
    for obs_name in args.observables:
        if obs_name not in observables:
            logger.warning(f"Observable '{obs_name}' not computed, skipping")
            continue
        
        # Create output filename
        filename = f"{args.prefix}{obs_name}_{field_id:03d}.fits" if isinstance(field_id, int) \
                   else f"{args.prefix}{obs_name}.fits"
        output_path = args.outdir / filename
        
        # Check if exists
        if output_path.exists() and not args.overwrite:
            logger.warning(f"Output file exists: {output_path}, skipping (use --overwrite)")
            continue
        
        # Save
        save_observable(
            observables[obs_name],
            output_path,
            redshift,
            pixel_scale,
            obs_name
        )
    
    logger.info("Processing complete!")
    
    # Print summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Input:  {deflection_file}")
    print(f"Output: {args.outdir}/")
    print(f"Redshift: {redshift:.4f}")
    print(f"Pixel scale: {pixel_scale:.4f} arcsec")
    print(f"Shape: {deflection.shape[0]} x {deflection.shape[1]}")
    print(f"Observables computed: {', '.join(args.observables)}")
    print("="*70)


if __name__ == '__main__':
    main()
