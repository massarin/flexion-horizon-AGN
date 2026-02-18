"""
Lensing observables from deflection field Jacobian.

Computes convergence (κ), shear (γ₁, γ₂), and flexion (F₁, F₂, G₁, G₂)
from the Jacobian matrix and second derivatives.

References:
    Bartelmann & Schneider (2001): "Weak gravitational lensing"
        Physics Reports 340, 291-472
        - Convergence: Eq. 3.47
        - Shear: Eq. 3.48

    Bacon et al. (2006): "Measuring gravitational spin"
        MNRAS 365, 414-428
        - Flexion F (spin-1): Eq. 8
        - Flexion G (spin-3): Eq. 9

    Goldberg & Bacon (2005): "Galaxy distortion by flexion"
        ApJ 619, 741-748
"""

import logging
from typing import Dict

import numpy as np

logger = logging.getLogger(__name__)


def convergence(jac: np.ndarray) -> np.ndarray:
    """
    Compute convergence κ from Jacobian.

    κ = 0.5 * (∂α₁/∂x₁ + ∂α₂/∂x₂)

    The convergence measures the isotropic magnification of images.

    Args:
        jac: Jacobian array from derivatives.compute_jacobian
             Shape (nx, ny, 4) or (nx, ny, 12)

    Returns:
        Convergence map with shape (nx, ny)

    References:
        Bartelmann & Schneider (2001) Eq. 3.47
        raytrace_tk.jl lines 659-668
    """
    kappa = 0.5 * (jac[:, :, 0] + jac[:, :, 3])

    logger.debug(
        f"Computed convergence: min={kappa.min():.3e}, "
        f"max={kappa.max():.3e}, mean={kappa.mean():.3e}"
    )

    return kappa


def shear_1(jac: np.ndarray) -> np.ndarray:
    """
    Compute shear component γ₁ from Jacobian.

    γ₁ = 0.5 * (∂α₁/∂x₁ - ∂α₂/∂x₂)

    The shear measures the anisotropic distortion of images.
    γ₁ corresponds to stretching along x/y axes.

    Args:
        jac: Jacobian array from derivatives.compute_jacobian

    Returns:
        Shear γ₁ map with shape (nx, ny)

    References:
        Bartelmann & Schneider (2001) Eq. 3.48
        raytrace_tk.jl lines 670-679
    """
    gamma1 = 0.5 * (jac[:, :, 0] - jac[:, :, 3])

    logger.debug(
        f"Computed shear γ₁: min={gamma1.min():.3e}, "
        f"max={gamma1.max():.3e}, mean={gamma1.mean():.3e}"
    )

    return gamma1


def shear_2(jac: np.ndarray) -> np.ndarray:
    """
    Compute shear component γ₂ from Jacobian.

    γ₂ = 0.5 * (∂α₁/∂x₂ + ∂α₂/∂x₁)

    The shear measures the anisotropic distortion of images.
    γ₂ corresponds to stretching at 45° angle.

    Args:
        jac: Jacobian array from derivatives.compute_jacobian

    Returns:
        Shear γ₂ map with shape (nx, ny)

    References:
        Bartelmann & Schneider (2001) Eq. 3.48
        raytrace_tk.jl lines 681-690
    """
    gamma2 = 0.5 * (jac[:, :, 1] + jac[:, :, 2])

    logger.debug(
        f"Computed shear γ₂: min={gamma2.min():.3e}, "
        f"max={gamma2.max():.3e}, mean={gamma2.mean():.3e}"
    )

    return gamma2


def shear_magnitude(jac: np.ndarray) -> np.ndarray:
    """
    Compute shear magnitude |γ| = sqrt(γ₁² + γ₂²).

    Args:
        jac: Jacobian array from derivatives.compute_jacobian

    Returns:
        Shear magnitude map with shape (nx, ny)
    """
    g1 = shear_1(jac)
    g2 = shear_2(jac)
    return np.sqrt(g1**2 + g2**2)


def flexion_F1(jac: np.ndarray) -> np.ndarray:
    """
    Compute first flexion component F₁ from extended Jacobian.

    F₁ = 0.5 * (∂²α₁/∂x₁² - ∂²α₂/∂x₁∂x₂) + 0.5 * (∂²α₂/∂x₂² + ∂²α₁/∂x₂∂x₁)

    Flexion measures the third-order lensing distortions (banana-shaped
    image distortions). F is the spin-1 component.

    Args:
        jac: Extended Jacobian array from derivatives.compute_jacobian
             Must have shape (nx, ny, 12) with second derivatives

    Returns:
        Flexion F₁ map with shape (nx, ny)

    Raises:
        ValueError: If jac doesn't contain second derivatives

    References:
        Bacon et al. (2006) MNRAS 365, 414, Eq. 8
        raytrace_tk.jl lines 692-701
    """
    if jac.shape[2] < 12:
        raise ValueError(
            f"Extended Jacobian required for flexion (need 12 components, got {jac.shape[2]})"
        )

    # Indices: 4=d1d1a1, 5=d2d2a1, 6=d1d1a2, 7=d2d2a2, 8=d2d1a1, 9=d2d1a2, 10=d1d2a1, 11=d1d2a2
    # Note: Julia uses 1-based indexing, Python uses 0-based
    # jac[i,j,5] in Julia = jac[:,:,4] in Python (d1d1a1 = ∂²α₁/∂x₁²)

    F1 = 0.5 * (jac[:, :, 4] - jac[:, :, 11]) + 0.5 * (jac[:, :, 6] + jac[:, :, 9])

    logger.debug(
        f"Computed flexion F₁: min={F1.min():.3e}, " f"max={F1.max():.3e}, mean={F1.mean():.3e}"
    )

    return F1


def flexion_F2(jac: np.ndarray) -> np.ndarray:
    """
    Compute second flexion component F₂ from extended Jacobian.

    F₂ = 0.5 * (∂²α₁/∂x₂² + ∂²α₂/∂x₂∂x₁) - 0.5 * (∂²α₁/∂x₁∂x₂ - ∂²α₂/∂x₁²)

    Args:
        jac: Extended Jacobian array with second derivatives

    Returns:
        Flexion F₂ map with shape (nx, ny)

    References:
        Bacon et al. (2006) MNRAS 365, 414, Eq. 8
        raytrace_tk.jl lines 703-712
    """
    if jac.shape[2] < 12:
        raise ValueError(
            f"Extended Jacobian required for flexion (need 12 components, got {jac.shape[2]})"
        )

    F2 = 0.5 * (jac[:, :, 5] + jac[:, :, 10]) - 0.5 * (jac[:, :, 8] - jac[:, :, 7])

    logger.debug(
        f"Computed flexion F₂: min={F2.min():.3e}, " f"max={F2.max():.3e}, mean={F2.mean():.3e}"
    )

    return F2


def flexion_G1(jac: np.ndarray) -> np.ndarray:
    """
    Compute first G-flexion component G₁ from extended Jacobian.

    G₁ = 0.5 * (∂²α₁/∂x₁² - ∂²α₂/∂x₁∂x₂) - 0.5 * (∂²α₁/∂x₂∂x₁ + ∂²α₂/∂x₂²)

    G is the spin-3 flexion component, complementary to F.

    Args:
        jac: Extended Jacobian array with second derivatives

    Returns:
        Flexion G₁ map with shape (nx, ny)

    References:
        Bacon et al. (2006) MNRAS 365, 414, Eq. 9
        raytrace_tk.jl lines 714-723
    """
    if jac.shape[2] < 12:
        raise ValueError(
            f"Extended Jacobian required for flexion (need 12 components, got {jac.shape[2]})"
        )

    G1 = 0.5 * (jac[:, :, 4] - jac[:, :, 11]) - 0.5 * (jac[:, :, 9] + jac[:, :, 6])

    logger.debug(
        f"Computed flexion G₁: min={G1.min():.3e}, " f"max={G1.max():.3e}, mean={G1.mean():.3e}"
    )

    return G1


def flexion_G2(jac: np.ndarray) -> np.ndarray:
    """
    Compute second G-flexion component G₂ from extended Jacobian.

    G₂ = 0.5 * (∂²α₁/∂x₂² + ∂²α₂/∂x₂∂x₁) + 0.5 * (∂²α₁/∂x₁∂x₂ - ∂²α₂/∂x₁²)

    Args:
        jac: Extended Jacobian array with second derivatives

    Returns:
        Flexion G₂ map with shape (nx, ny)

    References:
        Bacon et al. (2006) MNRAS 365, 414, Eq. 9
        raytrace_tk.jl lines 725-734
    """
    if jac.shape[2] < 12:
        raise ValueError(
            f"Extended Jacobian required for flexion (need 12 components, got {jac.shape[2]})"
        )

    G2 = 0.5 * (jac[:, :, 5] + jac[:, :, 10]) + 0.5 * (jac[:, :, 8] - jac[:, :, 7])

    logger.debug(
        f"Computed flexion G₂: min={G2.min():.3e}, " f"max={G2.max():.3e}, mean={G2.mean():.3e}"
    )

    return G2


def rotation(jac: np.ndarray) -> np.ndarray:
    """
    Compute rotation component ω from Jacobian.

    ω = 0.5 * (∂α₁/∂x₂ - ∂α₂/∂x₁)

    Rotation should be zero for pure lensing (no curl).
    Non-zero values indicate numerical errors or non-lensing effects.

    Args:
        jac: Jacobian array from derivatives.compute_jacobian

    Returns:
        Rotation map with shape (nx, ny)

    References:
        raytrace_tk.jl lines 736-745
    """
    omega = 0.5 * (jac[:, :, 1] - jac[:, :, 2])

    logger.debug(
        f"Computed rotation: min={omega.min():.3e}, "
        f"max={omega.max():.3e}, mean={omega.mean():.3e}, "
        f"rms={np.sqrt(np.mean(omega**2)):.3e}"
    )

    return omega


def compute_all_observables(jac: np.ndarray) -> Dict[str, np.ndarray]:
    """
    Compute all lensing observables from Jacobian.

    Args:
        jac: Jacobian or extended Jacobian array

    Returns:
        Dictionary with keys:
        - 'kappa': convergence
        - 'gamma1': shear component 1
        - 'gamma2': shear component 2
        - 'gamma': shear magnitude
        - 'rot': rotation (curl check)

        If jac contains second derivatives (shape[2] >= 12):
        - 'F1': flexion F component 1
        - 'F2': flexion F component 2
        - 'G1': flexion G component 1
        - 'G2': flexion G component 2
    """
    observables = {
        "kappa": convergence(jac),
        "gamma1": shear_1(jac),
        "gamma2": shear_2(jac),
        "gamma": shear_magnitude(jac),
        "rot": rotation(jac),
    }

    if jac.shape[2] >= 12:
        observables.update(
            {
                "F1": flexion_F1(jac),
                "F2": flexion_F2(jac),
                "G1": flexion_G1(jac),
                "G2": flexion_G2(jac),
            }
        )
        logger.info("Computed all observables (including flexion)")
    else:
        logger.info("Computed basic observables (kappa, gamma, rot)")

    return observables
