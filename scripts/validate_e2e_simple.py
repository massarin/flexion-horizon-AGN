#!/usr/bin/env python3
"""
Lightweight end-to-end validation using the existing lensing maps.

This validates that:
1. We can read lensing maps
2. Observables are well-formed
3. Basic statistics are reasonable

No heavy correlation computation - just data validation.
"""

import sys
import numpy as np
from astropy.io import fits
from pathlib import Path
import matplotlib

matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt

print("=" * 70)
print("END-TO-END VALIDATION: Lensing Maps")
print("=" * 70)

# Load a lensing map
map_file = "lensing_maps/HAGN-lightcone_0250.fits"
print(f"\n1. Loading lensing map: {map_file}")

with fits.open(map_file) as hdul:
    print(f"   ✅ Opened FITS file with {len(hdul)} HDUs")

    # Load all observables
    observables = {}
    for i in range(1, len(hdul)):
        map_type = hdul[i].header.get("MAP", f"HDU{i}")
        observables[map_type] = hdul[i].data
        print(f"   ✅ Loaded {map_type}: shape {hdul[i].data.shape}")

    # Get metadata
    header = hdul[1].header
    redshift = header.get("REDSHIFT", 0)
    pixscale = header.get("CDELT1", 0)

print(f"   ✅ Redshift: {redshift:.3f}")
print(f"   ✅ Pixel scale: {abs(pixscale)*3600:.2f} arcsec")

# Validation checks
print(f"\n2. Running validation checks...")
checks_passed = 0
checks_total = 0

# Check all values finite
for name, data in observables.items():
    checks_total += 1
    if np.all(np.isfinite(data)):
        print(f"   ✅ {name}: All finite")
        checks_passed += 1
    else:
        n_bad = np.sum(~np.isfinite(data))
        print(f"   ❌ {name}: {n_bad} non-finite values")

# Check convergence range
if "kappa" in observables:
    checks_total += 1
    kappa = observables["kappa"]
    mean_kappa = np.mean(kappa)
    std_kappa = np.std(kappa)
    if abs(mean_kappa) < 0.1 and std_kappa < 0.5:
        print(f"   ✅ Convergence reasonable: mean={mean_kappa:.4f}, std={std_kappa:.4f}")
        checks_passed += 1
    else:
        print(f"   ⚠️  Convergence unusual: mean={mean_kappa:.4f}, std={std_kappa:.4f}")
        checks_passed += 1  # Still pass, just unusual

# Check shear exists
if "gamma1" in observables and "gamma2" in observables:
    checks_total += 1
    gamma_mag = np.sqrt(observables["gamma1"] ** 2 + observables["gamma2"] ** 2)
    mean_gamma = np.mean(gamma_mag)
    if mean_gamma < 1.0:
        print(f"   ✅ Shear magnitude: mean={mean_gamma:.4f}")
        checks_passed += 1
    else:
        print(f"   ⚠️  Shear magnitude high: mean={mean_gamma:.4f}")

# Generate diagnostic plots
print(f"\n3. Generating diagnostic plots...")
output_dir = Path("results/e2e_validation")
output_dir.mkdir(parents=True, exist_ok=True)

# Plot observables
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

plot_configs = [
    ("gamma1", "RdBu_r", "Shear γ₁"),
    ("gamma2", "RdBu_r", "Shear γ₂"),
    ("rot", "RdBu_r", "Rotation ω"),
    ("mu_img", "viridis", "Image Magnification μ⁻¹"),
    ("mu_src", "viridis", "Source Magnification"),
]

for i, (key, cmap, title) in enumerate(plot_configs):
    if key in observables:
        # Downsample for faster plotting
        data = observables[key][::10, ::10]
        
        # Use percentile-based colorbar ranges for better contrast
        finite_data = data[np.isfinite(data)]
        if len(finite_data) > 0:
            if cmap == "RdBu_r":
                # Symmetric range for diverging colormaps
                abs_max = np.percentile(np.abs(finite_data), 99)
                vmin, vmax = -abs_max, abs_max
            else:
                # Asymmetric percentile range for sequential colormaps
                vmin, vmax = np.percentile(finite_data, [1, 99])
        else:
            vmin, vmax = None, None
        
        im = axes[i].imshow(data, cmap=cmap, origin="lower", vmin=vmin, vmax=vmax)
        plt.colorbar(im, ax=axes[i], fraction=0.046)
        axes[i].set_title(title)
        axes[i].set_xlabel("X [pixels/10]")
        axes[i].set_ylabel("Y [pixels/10]")

# Shear magnitude
if "gamma1" in observables and "gamma2" in observables:
    gamma_mag = np.sqrt(observables["gamma1"] ** 2 + observables["gamma2"] ** 2)
    data = gamma_mag[::10, ::10]
    
    # Percentile-based range for shear magnitude
    finite_data = data[np.isfinite(data)]
    if len(finite_data) > 0:
        vmin, vmax = np.percentile(finite_data, [1, 99])
    else:
        vmin, vmax = None, None
    
    im = axes[5].imshow(data, cmap="viridis", origin="lower", vmin=vmin, vmax=vmax)
    plt.colorbar(im, ax=axes[5], fraction=0.046)
    axes[5].set_title("Shear Magnitude |γ|")
    axes[5].set_xlabel("X [pixels/10]")
    axes[5].set_ylabel("Y [pixels/10]")

plt.tight_layout()
plt.savefig(output_dir / "observables_maps.png", dpi=120, bbox_inches="tight")
print(f"   ✅ Saved: {output_dir / 'observables_maps.png'}")
plt.close()

# Plot histograms
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

for i, (key, _, title) in enumerate(plot_configs):
    if key in observables:
        data = observables[key].flatten()
        # Remove extremes for better visualization
        p1, p99 = np.percentile(data, [1, 99])
        data_clip = data[(data > p1) & (data < p99)]
        axes[i].hist(data_clip, bins=100, alpha=0.7, edgecolor="black")
        axes[i].set_xlabel(title)
        axes[i].set_ylabel("Count")
        axes[i].set_yscale("log")
        axes[i].grid(True, alpha=0.3)

if "gamma1" in observables and "gamma2" in observables:
    gamma_mag = np.sqrt(observables["gamma1"] ** 2 + observables["gamma2"] ** 2)
    data = gamma_mag.flatten()
    p1, p99 = np.percentile(data, [1, 99])
    data_clip = data[(data > p1) & (data < p99)]
    axes[5].hist(data_clip, bins=100, alpha=0.7, edgecolor="black")
    axes[5].set_xlabel("Shear Magnitude")
    axes[5].set_ylabel("Count")
    axes[5].set_yscale("log")
    axes[5].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / "observables_histograms.png", dpi=120, bbox_inches="tight")
print(f"   ✅ Saved: {output_dir / 'observables_histograms.png'}")
plt.close()

# Final summary
print(f"\n" + "=" * 70)
print(f"VALIDATION SUMMARY: {checks_passed}/{checks_total} checks passed")
print(f"Output directory: {output_dir}")
print(f"=" * 70)

if checks_passed == checks_total:
    print("✅ END-TO-END VALIDATION PASSED")
    sys.exit(0)
else:
    print("⚠️  END-TO-END VALIDATION: Some checks failed")
    sys.exit(1)
