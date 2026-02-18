"""
Finite-difference derivatives for Jacobian computation.

Computes spatial derivatives of deflection fields using 4th-order
accurate finite differences with special boundary handling.

The Jacobian matrix is:
    J = ∂αᵢ/∂xⱼ = [[∂α₁/∂x₁, ∂α₁/∂x₂],
                    [∂α₂/∂x₁, ∂α₂/∂x₂]]

Extended Jacobian includes second derivatives for flexion:
    - First derivatives: ∂α₁/∂x₁, ∂α₁/∂x₂, ∂α₂/∂x₁, ∂α₂/∂x₂
    - Second derivatives: ∂²α₁/∂x₁², ∂²α₁/∂x₂², ∂²α₂/∂x₁², ∂²α₂/∂x₂²
    - Cross derivatives: ∂²α₁/∂x₁∂x₂, ∂²α₂/∂x₁∂x₂
"""

import logging
from typing import Tuple

import numpy as np

logger = logging.getLogger(__name__)


def derivative_1d(x: np.ndarray, axis: int = 0) -> np.ndarray:
    """
    Compute 1D finite-difference derivative along specified axis.

    Uses 4th-order accurate centered differences in the interior,
    with 3rd-order forward/backward differences at boundaries.

    Interior (i=2 to n-3):
        df/dx[i] = (-f[i-2] + 8*f[i-1] - 8*f[i+1] + f[i+2]) / 12

    Boundaries use 3rd-order one-sided formulas:
        df/dx[0] = (-3*f[0] + 4*f[1] - f[2]) / 2
        df/dx[1] = (-3*f[1] + 4*f[2] - f[3]) / 2
        df/dx[n-2] = (3*f[n-2] - 4*f[n-3] + f[n-4]) / 2
        df/dx[n-1] = (3*f[n-1] - 4*f[n-2] + f[n-3]) / 2

    Args:
        x: Input array (can be multi-dimensional)
        axis: Axis along which to compute derivative

    Returns:
        Derivative array (same shape as input)

    References:
        Fornberg (1988): "Generation of finite difference formulas"
        raytrace_tk.jl lines 146-157
    """
    n = x.shape[axis]
    result = np.zeros_like(x)

    # Move axis to front for easier indexing
    x_rolled = np.moveaxis(x, axis, 0)
    result_rolled = np.moveaxis(result, axis, 0)

    # Boundary: first point (forward difference, 3rd-order)
    result_rolled[0] = -3 * x_rolled[0] + 4 * x_rolled[1] - x_rolled[2]

    # Boundary: second point
    result_rolled[1] = -3 * x_rolled[1] + 4 * x_rolled[2] - x_rolled[3]

    # Interior: 4th-order centered difference
    for i in range(2, n - 2):
        result_rolled[i] = (
            -x_rolled[i - 2] + 8 * x_rolled[i - 1] - 8 * x_rolled[i + 1] + x_rolled[i + 2]
        )

    # Boundary: second-to-last point
    result_rolled[n - 2] = 3 * x_rolled[n - 2] - 4 * x_rolled[n - 3] + x_rolled[n - 4]

    # Boundary: last point (backward difference)
    result_rolled[n - 1] = 3 * x_rolled[n - 1] - 4 * x_rolled[n - 2] + x_rolled[n - 3]

    # Apply global scaling: -1/12 (as in Julia code line 155)
    result_rolled = -result_rolled / 12.0

    # Move axis back
    result = np.moveaxis(result_rolled, 0, axis)

    return result


def compute_jacobian(
    alpha: np.ndarray, pixel_scale: float = 1.0, compute_second_derivatives: bool = True
) -> np.ndarray:
    """
    Compute Jacobian matrix and second derivatives from deflection field.

    Args:
        alpha: Deflection field array with shape (nx, ny, 2)
               alpha[:,:,0] = α₁ (x-component)
               alpha[:,:,1] = α₂ (y-component)
        pixel_scale: Pixel scale for proper derivative units (default 1.0)
        compute_second_derivatives: If True, compute extended Jacobian with
                                    second derivatives for flexion (default True)

    Returns:
        If compute_second_derivatives=False:
            Jacobian array with shape (nx, ny, 4):
            [:,:,0] = ∂α₁/∂x₁
            [:,:,1] = ∂α₁/∂x₂
            [:,:,2] = ∂α₂/∂x₁
            [:,:,3] = ∂α₂/∂x₂

        If compute_second_derivatives=True:
            Extended Jacobian array with shape (nx, ny, 12):
            [:,:,0] = ∂α₁/∂x₁
            [:,:,1] = ∂α₁/∂x₂
            [:,:,2] = ∂α₂/∂x₁
            [:,:,3] = ∂α₂/∂x₂
            [:,:,4] = ∂²α₁/∂x₁²
            [:,:,5] = ∂²α₁/∂x₂²
            [:,:,6] = ∂²α₂/∂x₁²
            [:,:,7] = ∂²α₂/∂x₂²
            [:,:,8] = ∂²α₁/∂x₁∂x₂
            [:,:,9] = ∂²α₂/∂x₁∂x₂
            [:,:,10] = ∂²α₁/∂x₂∂x₁ (duplicate for validation)
            [:,:,11] = ∂²α₂/∂x₂∂x₁ (duplicate for validation)

    Raises:
        ValueError: If alpha doesn't have shape (nx, ny, 2)

    References:
        raytrace_tk.jl lines 218-326 (alpha2jac function)
    """
    if alpha.ndim != 3 or alpha.shape[2] != 2:
        raise ValueError(f"alpha must have shape (nx, ny, 2), got {alpha.shape}")

    nx, ny = alpha.shape[:2]

    # Scaling factors for derivatives
    # (number of pixels / physical scale)
    conv1 = nx / pixel_scale
    conv2 = ny / pixel_scale

    if compute_second_derivatives:
        jac = np.zeros((nx, ny, 12), dtype=np.float32)

        # Compute ∂α₁/∂x₁ and ∂²α₁/∂x₁²
        temp = alpha[:, :, 0] * conv1
        jac[:, :, 0] = derivative_1d(temp, axis=0)
        jac[:, :, 4] = derivative_1d(jac[:, :, 0], axis=0)

        # Compute ∂α₁/∂x₂ and ∂²α₁/∂x₂²
        temp = alpha[:, :, 0] * conv2
        jac[:, :, 1] = derivative_1d(temp, axis=1)
        jac[:, :, 5] = derivative_1d(jac[:, :, 1], axis=1)

        # Compute ∂α₂/∂x₁ and ∂²α₂/∂x₁²
        temp = alpha[:, :, 1] * conv1
        jac[:, :, 2] = derivative_1d(temp, axis=0)
        jac[:, :, 6] = derivative_1d(jac[:, :, 2], axis=0)

        # Compute ∂α₂/∂x₂ and ∂²α₂/∂x₂²
        temp = alpha[:, :, 1] * conv2
        jac[:, :, 3] = derivative_1d(temp, axis=1)
        jac[:, :, 7] = derivative_1d(jac[:, :, 3], axis=1)

        # Compute mixed derivatives ∂²α/∂x₁∂x₂
        # Take ∂/∂x₂ of (∂α₁/∂x₁)
        temp = jac[:, :, 0] / conv1
        jac[:, :, 8] = derivative_1d(temp, axis=1) * conv2

        # Take ∂/∂x₂ of (∂α₂/∂x₁)
        temp = jac[:, :, 2] / conv1
        jac[:, :, 9] = derivative_1d(temp, axis=1) * conv2

        # Compute mixed derivatives ∂²α/∂x₂∂x₁ (should equal above)
        # Take ∂/∂x₁ of (∂α₁/∂x₂)
        temp = jac[:, :, 1] / conv2
        jac[:, :, 10] = derivative_1d(temp, axis=0) * conv1

        # Take ∂/∂x₁ of (∂α₂/∂x₂)
        temp = jac[:, :, 3] / conv2
        jac[:, :, 11] = derivative_1d(temp, axis=0) * conv1

        logger.info(
            f"Computed extended Jacobian: shape={jac.shape}, "
            f"mixed derivative agreement: "
            f"mean(|∂²α₁/∂x₁∂x₂ - ∂²α₁/∂x₂∂x₁|)={np.mean(np.abs(jac[:,:,8] - jac[:,:,10])):.3e}"
        )
    else:
        jac = np.zeros((nx, ny, 4), dtype=np.float32)

        # Compute ∂α₁/∂x₁
        temp = alpha[:, :, 0] * conv1
        jac[:, :, 0] = derivative_1d(temp, axis=0)

        # Compute ∂α₁/∂x₂
        temp = alpha[:, :, 0] * conv2
        jac[:, :, 1] = derivative_1d(temp, axis=1)

        # Compute ∂α₂/∂x₁
        temp = alpha[:, :, 1] * conv1
        jac[:, :, 2] = derivative_1d(temp, axis=0)

        # Compute ∂α₂/∂x₂
        temp = alpha[:, :, 1] * conv2
        jac[:, :, 3] = derivative_1d(temp, axis=1)

        logger.info(f"Computed basic Jacobian: shape={jac.shape}")

    return jac


def jacobian_determinant(jac: np.ndarray) -> np.ndarray:
    """
    Compute determinant of Jacobian matrix.

    det(J) = (∂α₁/∂x₁)(∂α₂/∂x₂) - (∂α₁/∂x₂)(∂α₂/∂x₁)

    For weak lensing (no multiple images), det(J) should be positive everywhere.

    Args:
        jac: Jacobian array from compute_jacobian

    Returns:
        Determinant array with shape (nx, ny)
    """
    det = jac[:, :, 0] * jac[:, :, 3] - jac[:, :, 1] * jac[:, :, 2]

    # Check for negative determinants (multiple images)
    n_negative = np.sum(det < 0)
    if n_negative > 0:
        logger.warning(
            f"Found {n_negative} pixels with negative Jacobian determinant "
            f"(possible multiple images)"
        )

    return det
