#!/usr/bin/env python3
"""
End-to-end validation script for cosmo_lensing package.

This script validates the complete pipeline on real HAGN lensing data:
1. Load lensing map with all observables
2. Load galaxy catalog
3. Compute tangential shear correlations
4. Generate diagnostic plots
5. Validate outputs

Usage:
    python validate_e2e.py --map lensing_maps/HAGN-lightcone_0250.fits \\
                           --catalog Data/Galaxies_0-6_lensed.v2.0_cut_i27.fits \\
                           --output results/e2e_validation/
"""

import argparse
import sys
import time
from pathlib import Path
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Import cosmo_lensing modules
try:
    from cosmo_lensing import io, correlations, cosmology
except ImportError as e:
    logger.error(f"Failed to import cosmo_lensing: {e}")
    logger.error("Make sure package is installed: pip install -e .")
    sys.exit(1)


# Use consolidated loaders from io module
load_lensing_map = io.load_lensing_map
load_galaxy_catalog = io.load_galaxy_catalog



def compute_tangential_shear(lensing_data, galaxy_cat, metadata, nbins=15):
    """Compute tangential shear correlation around galaxies."""
    logger.info("Computing tangential shear correlations...")
    
    # Extract shear maps
    if 'gamma1' not in lensing_data or 'gamma2' not in lensing_data:
        raise ValueError("Lensing map must contain gamma1 and gamma2")
    
    gamma1 = lensing_data['gamma1']
    gamma2 = lensing_data['gamma2']
    
    ny, nx = gamma1.shape
    logger.info(f"  Shear map size: {nx} × {ny}")
    
    # For simplicity, flatten and assign to "source" positions
    # In real analysis, would cross-match catalog positions with map pixels
    # Here we'll use the catalog as "lenses" and compute correlation
    
    # Extract galaxy positions
    ra_lens = galaxy_cat['ra']
    dec_lens = galaxy_cat['dec']
    n_lens = len(ra_lens)
    
    if n_lens == 0:
        raise ValueError("No galaxies in catalog")
    
    logger.info(f"  Using {n_lens} galaxies as lenses")
    
    # For this demo, create mock "source" positions by sampling the shear field
    # In real analysis, would use actual source galaxy positions
    logger.info("  Creating mock source sample from shear field...")
    n_sources = min(10000, nx * ny // 100)  # Sample 1% of pixels
    
    # Random positions in image plane
    x_sources = np.random.uniform(0, nx, n_sources)
    y_sources = np.random.uniform(0, ny, n_sources)
    
    # Interpolate shear at source positions
    x_idx = x_sources.astype(int)
    y_idx = y_sources.astype(int)
    x_idx = np.clip(x_idx, 0, nx-1)
    y_idx = np.clip(y_idx, 0, ny-1)
    
    g1_sources = gamma1[y_idx, x_idx]
    g2_sources = gamma2[y_idx, x_idx]
    
    # Convert to RA/Dec (simple linear approximation for small fields)
    pixscale_deg = float(metadata['pixscale']) if isinstance(metadata['pixscale'], (int, float, str)) else 1.0 / 3600  # fallback
    if isinstance(pixscale_deg, str):
        try:
            pixscale_deg = float(pixscale_deg)
        except:
            pixscale_deg = 1.0 / 3600
    
    # Get map center from header
    header = metadata.get('header', {})
    ra_center = header.get('CRVAL1', 0.0)
    dec_center = header.get('CRVAL2', 0.0)
    
    # Convert pixel to RA/Dec
    ra_sources = ra_center + (x_sources - nx/2) * pixscale_deg
    dec_sources = dec_center + (y_sources - ny/2) * pixscale_deg
    
    # Subsample lenses for speed
    if n_lens > 1000:
        logger.info(f"  Subsampling lenses to 1000 for speed...")
        idx = np.random.choice(n_lens, 1000, replace=False)
        ra_lens = ra_lens[idx]
        dec_lens = dec_lens[idx]
    
    # Setup correlation
    corr = correlations.TangentialCorrelation(
        rmin=0.01, rmax=1.0, nbins=nbins,
        sep_units='degrees',
        var_method='shot'  # Use shot noise for speed
    )
    
    logger.info(f"  Computing correlation with {len(ra_lens)} lenses and {n_sources} sources...")
    start_time = time.time()
    
    result = corr.compute(
        ra_lens=ra_lens,
        dec_lens=dec_lens,
        ra_src=ra_sources,
        dec_src=dec_sources,
        g1=g1_sources,
        g2=g2_sources
    )
    
    elapsed = time.time() - start_time
    logger.info(f"  Correlation computed in {elapsed:.2f} seconds")
    
    return result


def generate_diagnostic_plots(result, lensing_data, metadata, output_dir):
    """Generate diagnostic plots."""
    logger.info("Generating diagnostic plots...")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Plot 1: Convergence map
    fig, ax = plt.subplots(figsize=(10, 8))
    if 'kappa' in lensing_data:
        im = ax.imshow(lensing_data['kappa'], cmap='RdBu_r', origin='lower')
        plt.colorbar(im, ax=ax, label='κ')
        ax.set_title(f"Convergence Map (z={metadata['redshift']})")
        ax.set_xlabel('Pixel X')
        ax.set_ylabel('Pixel Y')
        plt.savefig(output_dir / 'convergence_map.png', dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"  Saved: {output_dir / 'convergence_map.png'}")
    
    # Plot 2: Shear magnitude
    if 'gamma1' in lensing_data and 'gamma2' in lensing_data:
        gamma_mag = np.sqrt(lensing_data['gamma1']**2 + lensing_data['gamma2']**2)
        fig, ax = plt.subplots(figsize=(10, 8))
        im = ax.imshow(gamma_mag, cmap='viridis', origin='lower')
        plt.colorbar(im, ax=ax, label='|γ|')
        ax.set_title(f"Shear Magnitude (z={metadata['redshift']})")
        ax.set_xlabel('Pixel X')
        ax.set_ylabel('Pixel Y')
        plt.savefig(output_dir / 'shear_magnitude.png', dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"  Saved: {output_dir / 'shear_magnitude.png'}")
    
    # Plot 3: Tangential shear profile
    if result is not None:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.errorbar(result['r'], result['xi'], yerr=result['xi_err'],
                   marker='o', ls='', capsize=3, label='Tangential shear')
        ax.axhline(0, color='k', ls='--', alpha=0.3)
        ax.set_xlabel('Separation [degrees]')
        ax.set_ylabel('Tangential shear ξₜ')
        ax.set_title(f"Tangential Shear Profile (z={metadata['redshift']})")
        ax.set_xscale('log')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.savefig(output_dir / 'tangential_shear.png', dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"  Saved: {output_dir / 'tangential_shear.png'}")
    
    # Plot 4: Observable histograms
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    for i, (key, data) in enumerate(lensing_data.items()):
        if i >= 6:
            break
        ax = axes[i]
        ax.hist(data.flatten(), bins=100, alpha=0.7)
        ax.set_xlabel(key)
        ax.set_ylabel('Count')
        ax.set_title(f'{key} distribution')
        ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'observable_histograms.png', dpi=150, bbox_inches='tight')
    plt.close()
    logger.info(f"  Saved: {output_dir / 'observable_histograms.png'}")


def validate_outputs(lensing_data, result, metadata):
    """Validate outputs are sensible."""
    logger.info("Validating outputs...")
    
    checks_passed = 0
    checks_total = 0
    
    # Check 1: All observables finite
    for key, data in lensing_data.items():
        checks_total += 1
        if np.all(np.isfinite(data)):
            logger.info(f"  ✅ {key}: All values finite")
            checks_passed += 1
        else:
            n_bad = np.sum(~np.isfinite(data))
            logger.warning(f"  ❌ {key}: {n_bad} non-finite values")
    
    # Check 2: Convergence reasonable range
    if 'kappa' in lensing_data:
        checks_total += 1
        kappa = lensing_data['kappa']
        kappa_mean = np.mean(kappa)
        kappa_std = np.std(kappa)
        if -0.5 < kappa_mean < 0.5 and kappa_std < 1.0:
            logger.info(f"  ✅ Convergence: mean={kappa_mean:.4f}, std={kappa_std:.4f}")
            checks_passed += 1
        else:
            logger.warning(f"  ❌ Convergence unusual: mean={kappa_mean:.4f}, std={kappa_std:.4f}")
    
    # Check 3: Shear magnitude reasonable
    if 'gamma1' in lensing_data and 'gamma2' in lensing_data:
        checks_total += 1
        gamma_mag = np.sqrt(lensing_data['gamma1']**2 + lensing_data['gamma2']**2)
        gamma_mean = np.mean(gamma_mag)
        if gamma_mean < 0.5:
            logger.info(f"  ✅ Shear magnitude: mean={gamma_mean:.4f}")
            checks_passed += 1
        else:
            logger.warning(f"  ❌ Shear magnitude high: mean={gamma_mean:.4f}")
    
    # Check 4: Correlation result structure
    if result is not None:
        checks_total += 3
        required_keys = ['r', 'xi', 'xi_err']
        for key in required_keys:
            if key in result and len(result[key]) > 0:
                logger.info(f"  ✅ Correlation '{key}' present with {len(result[key])} bins")
                checks_passed += 1
            else:
                logger.warning(f"  ❌ Correlation '{key}' missing or empty")
    
    logger.info(f"\nValidation: {checks_passed}/{checks_total} checks passed")
    return checks_passed == checks_total


def main():
    parser = argparse.ArgumentParser(description='End-to-end validation of cosmo_lensing package')
    parser.add_argument('--map', default='lensing_maps/HAGN-lightcone_0250.fits',
                       help='Path to lensing map FITS file')
    parser.add_argument('--catalog', default='Data/Galaxies_0-6_lensed.v2.0_cut_i27.fits',
                       help='Path to galaxy catalog FITS file')
    parser.add_argument('--output', default='results/e2e_validation',
                       help='Output directory for results')
    parser.add_argument('--z-min', type=float, default=0.5,
                       help='Minimum redshift for source selection')
    parser.add_argument('--z-max', type=float, default=1.5,
                       help='Maximum redshift for source selection')
    parser.add_argument('--subsample', type=int, default=5000,
                       help='Subsample catalog to N galaxies (for speed)')
    parser.add_argument('--nbins', type=int, default=15,
                       help='Number of radial bins for correlation')
    
    args = parser.parse_args()
    
    logger.info("=" * 70)
    logger.info("End-to-End Validation: cosmo_lensing package")
    logger.info("=" * 70)
    
    try:
        # Step 1: Load lensing map
        lensing_data, metadata = load_lensing_map(args.map)
        
        # Step 2: Load galaxy catalog
        galaxy_cat = load_galaxy_catalog(args.catalog, 
                                        z_min=args.z_min, 
                                        z_max=args.z_max,
                                        subsample=args.subsample)
        
        # Step 3: Compute tangential shear
        result = compute_tangential_shear(lensing_data, galaxy_cat, metadata, nbins=args.nbins)
        
        # Step 4: Generate diagnostic plots
        generate_diagnostic_plots(result, lensing_data, metadata, args.output)
        
        # Step 5: Validate outputs
        all_passed = validate_outputs(lensing_data, result, metadata)
        
        # Save results
        output_dir = Path(args.output)
        np.savez(output_dir / 'correlation_result.npz', **result)
        logger.info(f"\nResults saved to: {output_dir}")
        
        if all_passed:
            logger.info("\n✅ END-TO-END VALIDATION PASSED")
            return 0
        else:
            logger.warning("\n⚠️  END-TO-END VALIDATION: Some checks failed")
            return 1
            
    except Exception as e:
        logger.error(f"\n❌ END-TO-END VALIDATION FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
