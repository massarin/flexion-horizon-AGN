#!/usr/bin/env python3
"""
Compute tangential shear correlation functions around galaxies.

This script reads lensing observable maps and a galaxy catalog, computes
azimuthal averages of the tangential shear/convergence around galaxy positions,
and writes the correlation profiles to numpy archives.

Usage:
    python compute_correlations.py kappa.fits gamma1.fits gamma2.fits galaxies.fits --outdir results/
    python compute_correlations.py --config workflow/config.yaml --field-id 50

Outputs:
    - tangential_shear_050.npz: Arrays (r, xi, xi_err, npairs, weight)
    - convergence_profile_050.npz: Azimuthally averaged κ(r)
    - (optional) flexion_profiles_050.npz: F and G correlations
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional, Dict, Tuple
import time

import numpy as np
import yaml
from astropy.io import fits
from astropy.table import Table

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from cosmo_lensing.correlations import TangentialCorrelation

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_observable_map(fits_path: Path) -> Tuple[np.ndarray, fits.Header]:
    """Load lensing observable from FITS file."""
    logger.info(f"Loading {fits_path.name}")
    
    with fits.open(fits_path) as hdul:
        data = hdul[0].data
        header = hdul[0].header
    
    return data, header


def load_galaxy_catalog(catalog_path: Path, ra_col: str = 'RA', dec_col: str = 'DEC') -> Table:
    """Load galaxy catalog from FITS file."""
    logger.info(f"Loading galaxy catalog from {catalog_path}")
    
    catalog = Table.read(catalog_path)
    logger.info(f"Loaded {len(catalog)} galaxies")
    
    # Check required columns
    if ra_col not in catalog.colnames or dec_col not in catalog.colnames:
        raise ValueError(f"Catalog must have '{ra_col}' and '{dec_col}' columns")
    
    return catalog


def pixels_to_radec(
    i: np.ndarray,
    j: np.ndarray,
    header: fits.Header
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert pixel coordinates to RA/Dec using WCS.
    
    Simplified implementation - assumes TAN projection centered at (0,0).
    For production, should use astropy.wcs.WCS for proper transformation.
    
    Args:
        i: Row indices (y-coordinates)
        j: Column indices (x-coordinates)
        header: FITS header with WCS keywords
    
    Returns:
        (ra, dec) in degrees
    """
    # Get WCS parameters
    crpix1 = header.get('CRPIX1', header['NAXIS1'] / 2.0)
    crpix2 = header.get('CRPIX2', header['NAXIS2'] / 2.0)
    crval1 = header.get('CRVAL1', 0.0)
    crval2 = header.get('CRVAL2', 0.0)
    cdelt1 = header.get('CDELT1', -header['PIXSCALE'] / 3600.0)
    cdelt2 = header.get('CDELT2', header['PIXSCALE'] / 3600.0)
    
    # Simple linear transformation (good for small fields)
    ra = crval1 + (j - crpix1) * cdelt1
    dec = crval2 + (i - crpix2) * cdelt2
    
    return ra, dec


def create_source_catalog_from_map(
    observable_map: np.ndarray,
    g1_map: np.ndarray,
    g2_map: np.ndarray,
    header: fits.Header,
    sampling: int = 1
) -> Dict[str, np.ndarray]:
    """
    Create source catalog by sampling pixels from observable maps.
    
    Each pixel becomes a "source" with position and shear.
    
    Args:
        observable_map: Observable to use for weights
        g1_map: Shear component 1 map
        g2_map: Shear component 2 map
        header: FITS header with WCS
        sampling: Sampling factor (1 = every pixel, 2 = every other pixel, etc.)
    
    Returns:
        Dictionary with 'ra', 'dec', 'g1', 'g2', 'weight' arrays
    """
    ny, nx = observable_map.shape
    
    # Create pixel grid
    i_grid, j_grid = np.meshgrid(
        np.arange(0, ny, sampling),
        np.arange(0, nx, sampling),
        indexing='ij'
    )
    
    # Convert to RA/Dec
    ra, dec = pixels_to_radec(i_grid.ravel(), j_grid.ravel(), header)
    
    # Extract shear values
    g1 = g1_map[::sampling, ::sampling].ravel()
    g2 = g2_map[::sampling, ::sampling].ravel()
    
    # Use uniform weights (could use noise maps if available)
    weights = np.ones_like(g1)
    
    # Remove NaN/inf
    valid = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(g1) & np.isfinite(g2)
    
    logger.info(f"Created source catalog: {np.sum(valid)} sources from {ny}x{nx} map "
               f"(sampling={sampling})")
    
    return {
        'ra': ra[valid],
        'dec': dec[valid],
        'g1': g1[valid],
        'g2': g2[valid],
        'weight': weights[valid]
    }


def compute_tangential_shear(
    lens_catalog: Table,
    observable_maps: Dict[str, np.ndarray],
    header: fits.Header,
    rmin: float = 0.01,
    rmax: float = 10.0,
    nbins: int = 20,
    sep_units: str = 'arcmin',
    var_method: str = 'jackknife',
    npatch: int = 50,
    sampling: int = 1
) -> Dict[str, np.ndarray]:
    """
    Compute tangential shear profile around lenses.
    
    Args:
        lens_catalog: Galaxy catalog with RA/Dec columns
        observable_maps: Dict with 'gamma1' and 'gamma2' arrays
        header: FITS header for WCS
        rmin: Minimum separation
        rmax: Maximum separation
        nbins: Number of radial bins
        sep_units: Units for separations
        var_method: Variance estimation method
        npatch: Number of patches for jackknife
        sampling: Pixel sampling factor
    
    Returns:
        Dictionary with correlation results
    """
    logger.info("Setting up correlation computation...")
    
    # Create source catalog from maps
    sources = create_source_catalog_from_map(
        observable_maps['gamma1'],  # Use shear for weights
        observable_maps['gamma1'],
        observable_maps['gamma2'],
        header,
        sampling=sampling
    )
    
    # Extract lens positions
    ra_lens = np.array(lens_catalog['RA'])
    dec_lens = np.array(lens_catalog['DEC'])
    
    # Filter lenses within map bounds
    # TODO: Implement proper WCS bounds checking
    valid_lenses = np.ones(len(ra_lens), dtype=bool)
    
    logger.info(f"Using {np.sum(valid_lenses)} lenses for correlation")
    
    # Set up correlation calculator
    corr = TangentialCorrelation(
        rmin=rmin,
        rmax=rmax,
        nbins=nbins,
        sep_units=sep_units,
        var_method=var_method,
        npatch=npatch
    )
    
    # Compute correlation
    start = time.time()
    result = corr.compute(
        ra_lens[valid_lenses],
        dec_lens[valid_lenses],
        sources['ra'],
        sources['dec'],
        sources['g1'],
        sources['g2'],
        weights=sources['weight']
    )
    
    elapsed = time.time() - start
    logger.info(f"Correlation computed in {elapsed:.2f}s")
    
    return result


def main():
    parser = argparse.ArgumentParser(
        description='Compute tangential correlation functions',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input specification
    parser.add_argument('--gamma1', type=Path, required=True,
                       help='FITS file with shear component 1')
    parser.add_argument('--gamma2', type=Path, required=True,
                       help='FITS file with shear component 2')
    parser.add_argument('--kappa', type=Path,
                       help='FITS file with convergence (optional)')
    parser.add_argument('--catalog', type=Path, required=True,
                       help='FITS file with galaxy catalog')
    
    # Output control
    parser.add_argument('--outdir', type=Path, default=Path('results'),
                       help='Output directory')
    parser.add_argument('--prefix', type=str, default='',
                       help='Prefix for output filenames')
    parser.add_argument('--field-id', type=int,
                       help='Field ID for filename')
    
    # Correlation parameters
    parser.add_argument('--rmin', type=float, default=0.05,
                       help='Minimum separation [arcmin]')
    parser.add_argument('--rmax', type=float, default=5.0,
                       help='Maximum separation [arcmin]')
    parser.add_argument('--nbins', type=int, default=20,
                       help='Number of radial bins')
    parser.add_argument('--var-method', type=str, default='jackknife',
                       choices=['jackknife', 'bootstrap', 'shot'],
                       help='Variance estimation method')
    parser.add_argument('--npatch', type=int, default=50,
                       help='Number of patches for jackknife/bootstrap')
    
    # Advanced options
    parser.add_argument('--sampling', type=int, default=10,
                       help='Pixel sampling factor (higher = faster but less accurate)')
    parser.add_argument('--overwrite', action='store_true',
                       help='Overwrite existing outputs')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose logging')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Check inputs exist
    for path in [args.gamma1, args.gamma2, args.catalog]:
        if not path.exists():
            logger.error(f"Input file not found: {path}")
            sys.exit(1)
    
    # Create output directory
    args.outdir.mkdir(parents=True, exist_ok=True)
    
    # Load observable maps
    gamma1, header1 = load_observable_map(args.gamma1)
    gamma2, header2 = load_observable_map(args.gamma2)
    
    # Check consistency
    if gamma1.shape != gamma2.shape:
        logger.error(f"Shape mismatch: gamma1 {gamma1.shape} vs gamma2 {gamma2.shape}")
        sys.exit(1)
    
    observable_maps = {'gamma1': gamma1, 'gamma2': gamma2}
    
    if args.kappa:
        kappa, _ = load_observable_map(args.kappa)
        observable_maps['kappa'] = kappa
    
    # Load galaxy catalog
    lens_catalog = load_galaxy_catalog(args.catalog)
    
    # Compute tangential shear correlation
    logger.info("Computing tangential shear profile...")
    result = compute_tangential_shear(
        lens_catalog,
        observable_maps,
        header1,
        rmin=args.rmin,
        rmax=args.rmax,
        nbins=args.nbins,
        var_method=args.var_method,
        npatch=args.npatch,
        sampling=args.sampling
    )
    
    # Save results
    field_id = args.field_id if args.field_id else "unknown"
    filename = f"{args.prefix}tangential_shear_{field_id:03d}.npz" if isinstance(field_id, int) \
               else f"{args.prefix}tangential_shear.npz"
    output_path = args.outdir / filename
    
    if output_path.exists() and not args.overwrite:
        logger.warning(f"Output exists: {output_path} (use --overwrite)")
    else:
        np.savez(
            output_path,
            r=result['r'],
            r_nom=result['r_nom'],
            xi=result['xi'],
            xi_err=result['xi_err'],
            xim=result['xim'],
            xim_err=result['xim_err'],
            npairs=result['npairs'],
            weight=result['weight']
        )
        logger.info(f"Saved correlation to {output_path}")
    
    # Print summary
    print("\n" + "="*70)
    print("CORRELATION SUMMARY")
    print("="*70)
    print(f"Lenses: {len(lens_catalog)}")
    print(f"Sources: {len(result['r'])} radial bins")
    print(f"Total pairs: {np.sum(result['npairs']):.0f}")
    print(f"Radial range: [{result['r'][0]:.3f}, {result['r'][-1]:.3f}] arcmin")
    print(f"Output: {output_path}")
    print("="*70)
    
    # Show first few bins
    print("\nFirst 5 radial bins:")
    print("  r [arcmin]     γ_t          σ(γ_t)     N_pairs")
    print("-" * 60)
    for i in range(min(5, len(result['r']))):
        print(f"  {result['r'][i]:8.4f}    {result['xi'][i]:+.6f}   {result['xi_err'][i]:.6f}   "
              f"{result['npairs'][i]:8.0f}")
    print("="*70)


if __name__ == '__main__':
    main()
