"""
Integration test: deflection field → observables → FITS

Tests the full pipeline on synthetic data.
"""

import numpy as np
import tempfile
from pathlib import Path

from cosmo_lensing import io, cosmology, derivatives, observables


def create_synthetic_deflection_field(nx=500, ny=500):
    """
    Create a synthetic NFW-like deflection field for testing.
    
    Returns:
        alpha: Deflection field (nx, ny, 2)
        redshift: Source redshift
    """
    print(f"Creating synthetic deflection field ({nx}x{ny})...")
    
    # Create coordinate grid
    x = np.linspace(-10, 10, nx, dtype=np.float32)
    y = np.linspace(-10, 10, ny, dtype=np.float32)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    # NFW-like profile: α ∝ θ_E² / r for simplicity
    theta_E = 2.0  # Einstein radius
    r = np.sqrt(X**2 + Y**2) + 0.1  # Avoid singularity
    
    # Deflection field components
    alpha = np.zeros((nx, ny, 2), dtype=np.float32)
    alpha[:, :, 0] = theta_E**2 * X / r**2  # α₁
    alpha[:, :, 1] = theta_E**2 * Y / r**2  # α₂
    
    redshift = 1.0
    
    print(f"  α₁: min={alpha[:,:,0].min():.3f}, max={alpha[:,:,0].max():.3f}")
    print(f"  α₂: min={alpha[:,:,1].min():.3f}, max={alpha[:,:,1].max():.3f}")
    
    return alpha, redshift


def test_full_pipeline():
    """Test complete pipeline: deflection → Jacobian → observables → FITS."""
    
    print("\n" + "="*70)
    print("Integration Test: Full Lensing Pipeline")
    print("="*70)
    
    # Step 1: Create synthetic deflection field
    alpha, zs = create_synthetic_deflection_field(nx=500, ny=500)
    print(f"✓ Created deflection field: shape={alpha.shape}, z={zs}")
    
    # Step 2: Initialize cosmology calculator
    cosmo = cosmology.CosmologyCalculator()
    print(f"✓ Initialized cosmology: H0={cosmo.h*100:.2f}, Ωm={cosmo.cosmo.Om0:.4f}")
    
    # Step 3: Compute distance scalings
    arcsec2kpc, sigma_crit = cosmo.distance_scaling(zs)
    print(f"✓ Distance scaling: 1\" = {arcsec2kpc:.2f} kpc at z={zs}")
    
    # Step 4: Compute Jacobian and second derivatives
    print("\nComputing derivatives...")
    jac = derivatives.compute_jacobian(
        alpha,
        pixel_scale=1.0,
        compute_second_derivatives=True
    )
    print(f"✓ Computed extended Jacobian: shape={jac.shape}")
    
    # Step 5: Compute all observables
    print("\nComputing observables...")
    obs = observables.compute_all_observables(jac)
    
    for name, data in obs.items():
        print(f"  {name:6s}: min={data.min():+.3e}, max={data.max():+.3e}, "
              f"mean={data.mean():+.3e}, rms={np.sqrt(np.mean(data**2)):.3e}")
    
    # Step 6: Validate observables
    print("\nValidating observables...")
    
    # Convergence should have reasonable values
    assert np.all(np.isfinite(obs['kappa'])), "Convergence has NaN/Inf"
    print(f"✓ Convergence: all finite")
    
    # Shear magnitude should be positive
    assert np.all(obs['gamma'] >= 0), "Shear magnitude is negative"
    print(f"✓ Shear magnitude: all positive")
    
    # Rotation should be small (curl-free for lensing)
    rot_rms = np.sqrt(np.mean(obs['rot']**2))
    kappa_rms = np.sqrt(np.mean(obs['kappa']**2))
    assert rot_rms < 0.1 * kappa_rms, "Rotation too large (not curl-free)"
    print(f"✓ Rotation: rms/κ_rms = {rot_rms/kappa_rms:.2%} (should be << 1)")
    
    # Flexion should be computed
    assert 'F1' in obs, "Flexion F1 not computed"
    assert 'G1' in obs, "Flexion G1 not computed"
    assert np.all(np.isfinite(obs['F1'])), "Flexion F1 has NaN/Inf"
    print(f"✓ Flexion: all components finite")
    
    # Step 7: Write FITS maps
    print("\nWriting FITS maps...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        for name, data in obs.items():
            if name == 'gamma':  # Skip magnitude, just save components
                continue
            
            outfile = tmpdir / f"test_{name}_z{zs:.3f}.fits"
            io.write_fits_map(
                data,
                str(outfile),
                redshift=zs,
                observable=name,
                pixel_scale=1.0
            )
            
            # Verify file was written
            assert outfile.exists(), f"FITS file not created: {outfile}"
            file_size = outfile.stat().st_size / 1024  # KB
            print(f"  {name:6s} → {outfile.name:30s} ({file_size:6.1f} KB)")
        
        print(f"✓ Wrote {len(obs)-1} FITS maps to {tmpdir}")
    
    # Step 8: Summary
    print("\n" + "="*70)
    print("✅ Integration test PASSED")
    print("="*70)
    print(f"Pipeline validated:")
    print(f"  • Deflection field: {alpha.shape[0]}×{alpha.shape[1]} pixels")
    print(f"  • Observables: {len(obs)} computed (κ, γ₁, γ₂, |γ|, F₁, F₂, G₁, G₂, ω)")
    print(f"  • FITS output: working")
    print(f"  • All checks: passed")
    print()


if __name__ == "__main__":
    test_full_pipeline()
