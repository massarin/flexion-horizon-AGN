"""
Synthetic data generators for testing and validation.

Provides analytic lensing solutions for:
- Point mass lens
- Singular Isothermal Sphere (SIS)
- Navarro-Frenk-White (NFW) halo profile

These generators create deflection fields with known analytic properties,
enabling rigorous validation of the lensing pipeline against theory.

References:
    Schneider, Ehlers & Falco (1992): "Gravitational Lenses"
    Wright & Brainerd (2000): "Gravitational lensing by NFW halos"
"""

import logging
from typing import Tuple

import numpy as np

logger = logging.getLogger(__name__)


def generate_point_mass_deflection(
    theta_E: float, grid_size: int = 500, box_size: float = 10.0, regularization: float = 0.01
) -> Tuple[np.ndarray, float]:
    """
    Generate deflection field for a point mass lens.

    Analytic solution:
        α(θ) = θ_E² / θ

    where θ is the angular distance from the lens center.

    Convergence:
        κ(θ) = θ_E² / (2θ²)

    Shear:
        γ_t(θ) = θ_E² / (2θ²)

    Args:
        theta_E: Einstein radius (arcsec or dimensionless)
        grid_size: Number of pixels per side
        box_size: Physical size of box (same units as theta_E)
        regularization: Small value to avoid singularity at center

    Returns:
        Tuple of (deflection_field, redshift)
        - deflection_field: (grid_size, grid_size, 2) array
        - redshift: source redshift (set to 1.0 for testing)

    References:
        Schneider, Ehlers & Falco (1992) Eq. 3.52
    """
    logger.info(
        f"Generating point mass deflection: θ_E={theta_E:.3f}, "
        f"grid={grid_size}x{grid_size}, box={box_size:.1f}"
    )

    # Create coordinate grid
    x = np.linspace(-box_size / 2, box_size / 2, grid_size, dtype=np.float32)
    y = np.linspace(-box_size / 2, box_size / 2, grid_size, dtype=np.float32)
    X, Y = np.meshgrid(x, y, indexing="ij")

    # Radial distance from center (with regularization)
    theta = np.sqrt(X**2 + Y**2) + regularization

    # Deflection angle: α = θ_E² / θ * θ_hat
    alpha_magnitude = theta_E**2 / theta

    alpha = np.zeros((grid_size, grid_size, 2), dtype=np.float32)
    alpha[:, :, 0] = alpha_magnitude * X / theta  # α_x
    alpha[:, :, 1] = alpha_magnitude * Y / theta  # α_y

    logger.info(
        f"Point mass deflection: |α|_max={alpha_magnitude.max():.3f}, "
        f"α1: [{alpha[:,:,0].min():.2f}, {alpha[:,:,0].max():.2f}]"
    )

    return alpha, 1.0


def generate_sis_deflection(
    sigma_v: float, grid_size: int = 500, box_size: float = 10.0, regularization: float = 0.01
) -> Tuple[np.ndarray, float]:
    """
    Generate deflection field for a Singular Isothermal Sphere (SIS).

    Analytic solution:
        α(θ) = θ_E
        θ_E = 4π (σ_v/c)² (D_ls/D_s)

    Convergence:
        κ(θ) = θ_E / (2θ)

    Shear:
        γ_t(θ) = θ_E / (2θ)

    For testing, we use simplified θ_E = σ_v (normalized units).

    Args:
        sigma_v: Velocity dispersion (normalized, ~1-2 for typical galaxies)
        grid_size: Number of pixels per side
        box_size: Physical size of box
        regularization: Small value to avoid singularity at center

    Returns:
        Tuple of (deflection_field, redshift)

    References:
        Schneider, Ehlers & Falco (1992) Eq. 3.87
    """
    logger.info(f"Generating SIS deflection: σ_v={sigma_v:.3f}, " f"grid={grid_size}x{grid_size}")

    # Create coordinate grid
    x = np.linspace(-box_size / 2, box_size / 2, grid_size, dtype=np.float32)
    y = np.linspace(-box_size / 2, box_size / 2, grid_size, dtype=np.float32)
    X, Y = np.meshgrid(x, y, indexing="ij")

    # Radial distance from center
    theta = np.sqrt(X**2 + Y**2) + regularization

    # Einstein radius (simplified for testing)
    theta_E = sigma_v

    # Deflection angle: α = θ_E * θ_hat (constant magnitude)
    alpha = np.zeros((grid_size, grid_size, 2), dtype=np.float32)
    alpha[:, :, 0] = theta_E * X / theta  # α_x
    alpha[:, :, 1] = theta_E * Y / theta  # α_y

    logger.info(
        f"SIS deflection: θ_E={theta_E:.3f}, "
        f"α: [{alpha[:,:,0].min():.2f}, {alpha[:,:,0].max():.2f}]"
    )

    return alpha, 1.0


def generate_nfw_deflection(
    mass: float,
    concentration: float,
    redshift: float = 0.5,
    grid_size: int = 500,
    box_size: float = 10.0,
    regularization: float = 0.01,
) -> Tuple[np.ndarray, float]:
    """
    Generate deflection field for an NFW halo profile.

    The NFW profile is:
        ρ(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]

    where:
        r_s = scale radius
        ρ_s = characteristic density

    This is a simplified implementation for testing. Full implementation
    will be in Phase 3 (nfw.py module).

    For now, we approximate with a combination of point mass (inner region)
    and SIS (outer region) to have realistic-looking deflection field.

    Args:
        mass: Halo mass (in normalized units, ~1-10 for testing)
        concentration: Concentration parameter (typically 3-15)
        redshift: Source redshift
        grid_size: Number of pixels per side
        box_size: Physical size of box
        regularization: Small value to avoid singularity

    Returns:
        Tuple of (deflection_field, redshift)

    References:
        Wright & Brainerd (2000) ApJ 534, 34
        Bartelmann (1996) A&A 313, 697
    """
    logger.info(
        f"Generating NFW-like deflection: M={mass:.2f}, c={concentration:.1f}, "
        f"z={redshift:.2f}, grid={grid_size}x{grid_size}"
    )

    # Create coordinate grid
    x = np.linspace(-box_size / 2, box_size / 2, grid_size, dtype=np.float32)
    y = np.linspace(-box_size / 2, box_size / 2, grid_size, dtype=np.float32)
    X, Y = np.meshgrid(x, y, indexing="ij")

    # Radial distance
    theta = np.sqrt(X**2 + Y**2) + regularization

    # NFW approximation: smooth transition from point mass to SIS
    # r_s = scale radius (related to concentration)
    r_s = box_size / (4 * concentration)

    # Deflection magnitude (simplified NFW-like profile)
    # Inner region: ~ point mass, outer: ~ SIS
    x_scaled = theta / r_s

    # F(x) function approximation for NFW (simplified)
    # For testing purposes, use interpolation between point mass and SIS
    inner_weight = np.exp(-x_scaled)
    outer_weight = 1.0 - inner_weight

    # Combine point mass and SIS behavior
    alpha_magnitude = inner_weight * mass / theta + outer_weight * np.sqrt(  # Point mass-like
        mass
    )  # SIS-like

    alpha = np.zeros((grid_size, grid_size, 2), dtype=np.float32)
    alpha[:, :, 0] = alpha_magnitude * X / theta
    alpha[:, :, 1] = alpha_magnitude * Y / theta

    logger.info(
        f"NFW-like deflection: r_s={r_s:.2f}, "
        f"|α|: [{alpha_magnitude.min():.2f}, {alpha_magnitude.max():.2f}]"
    )

    return alpha, redshift


def generate_shear_field(
    shear: float, angle: float = 0.0, grid_size: int = 500, box_size: float = 10.0
) -> Tuple[np.ndarray, float]:
    """
    Generate a pure external shear field.

    This creates a deflection field corresponding to constant external shear,
    useful for testing shear extraction and validation.

    Deflection for pure shear:
        α₁ = γ₁ x + γ₂ y
        α₂ = γ₂ x - γ₁ y

    where:
        γ₁ = γ cos(2φ)
        γ₂ = γ sin(2φ)

    Args:
        shear: Shear magnitude (dimensionless, typically 0.01-0.1)
        angle: Shear orientation angle in radians
        grid_size: Number of pixels per side
        box_size: Physical size of box

    Returns:
        Tuple of (deflection_field, redshift)
    """
    logger.info(
        f"Generating pure shear field: γ={shear:.3f}, "
        f"φ={np.degrees(angle):.1f}°, grid={grid_size}x{grid_size}"
    )

    # Shear components
    gamma1 = shear * np.cos(2 * angle)
    gamma2 = shear * np.sin(2 * angle)

    # Create coordinate grid
    x = np.linspace(-box_size / 2, box_size / 2, grid_size, dtype=np.float32)
    y = np.linspace(-box_size / 2, box_size / 2, grid_size, dtype=np.float32)
    X, Y = np.meshgrid(x, y, indexing="ij")

    # Deflection field for pure shear
    alpha = np.zeros((grid_size, grid_size, 2), dtype=np.float32)
    alpha[:, :, 0] = gamma1 * X + gamma2 * Y
    alpha[:, :, 1] = gamma2 * X - gamma1 * Y

    logger.info(f"Shear field: γ₁={gamma1:.4f}, γ₂={gamma2:.4f}")

    return alpha, 1.0


def generate_convergence_field(
    kappa: float, grid_size: int = 500, box_size: float = 10.0
) -> Tuple[np.ndarray, float]:
    """
    Generate a pure convergence (uniform magnification) field.

    Deflection for pure convergence:
        α = κ * θ

    This creates uniform isotropic magnification, useful for testing
    convergence extraction.

    Args:
        kappa: Convergence value (dimensionless)
        grid_size: Number of pixels per side
        box_size: Physical size of box

    Returns:
        Tuple of (deflection_field, redshift)
    """
    logger.info(
        f"Generating pure convergence field: κ={kappa:.3f}, " f"grid={grid_size}x{grid_size}"
    )

    # Create coordinate grid
    x = np.linspace(-box_size / 2, box_size / 2, grid_size, dtype=np.float32)
    y = np.linspace(-box_size / 2, box_size / 2, grid_size, dtype=np.float32)
    X, Y = np.meshgrid(x, y, indexing="ij")

    # Deflection for pure convergence
    alpha = np.zeros((grid_size, grid_size, 2), dtype=np.float32)
    alpha[:, :, 0] = kappa * X
    alpha[:, :, 1] = kappa * Y

    logger.info(f"Convergence field: κ={kappa:.3f}")

    return alpha, 1.0
