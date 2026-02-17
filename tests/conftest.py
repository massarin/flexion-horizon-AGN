"""
Pytest configuration and fixtures for cosmo_lensing tests.
"""

import numpy as np
import pytest


@pytest.fixture
def simple_deflection_field():
    """Small synthetic deflection field for testing."""
    nx, ny = 100, 100
    x = np.linspace(-5, 5, nx)
    y = np.linspace(-5, 5, ny)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    # Point mass deflection: α = θ_E² / θ
    theta_E = 1.0
    r = np.sqrt(X**2 + Y**2) + 1e-10  # Avoid division by zero
    
    alpha = np.zeros((nx, ny, 2), dtype=np.float32)
    alpha[:, :, 0] = theta_E**2 * X / r**2
    alpha[:, :, 1] = theta_E**2 * Y / r**2
    
    return alpha, 1.0  # (deflection_field, redshift)


@pytest.fixture
def test_jacobian():
    """Small synthetic Jacobian for testing observables."""
    nx, ny = 50, 50
    jac = np.zeros((nx, ny, 4), dtype=np.float32)
    
    # Simple linear gradient
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    jac[:, :, 0] = 0.1 * X  # ∂α₁/∂x₁
    jac[:, :, 1] = 0.05 * Y  # ∂α₁/∂x₂
    jac[:, :, 2] = 0.05 * X  # ∂α₂/∂x₁
    jac[:, :, 3] = 0.1 * Y  # ∂α₂/∂x₂
    
    return jac


@pytest.fixture
def extended_jacobian():
    """Jacobian with second derivatives for flexion testing."""
    nx, ny = 50, 50
    jac = np.zeros((nx, ny, 12), dtype=np.float32)
    
    # Fill first derivatives
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    jac[:, :, 0] = 0.1 * X
    jac[:, :, 1] = 0.05 * Y
    jac[:, :, 2] = 0.05 * X
    jac[:, :, 3] = 0.1 * Y
    
    # Fill second derivatives (constant for simple case)
    jac[:, :, 4] = 0.1  # ∂²α₁/∂x₁²
    jac[:, :, 5] = 0.05  # ∂²α₁/∂x₂²
    jac[:, :, 6] = 0.05  # ∂²α₂/∂x₁²
    jac[:, :, 7] = 0.1  # ∂²α₂/∂x₂²
    jac[:, :, 8] = 0.02  # ∂²α₁/∂x₁∂x₂
    jac[:, :, 9] = 0.02  # ∂²α₂/∂x₁∂x₂
    jac[:, :, 10] = 0.02  # ∂²α₁/∂x₂∂x₁
    jac[:, :, 11] = 0.02  # ∂²α₂/∂x₂∂x₁
    
    return jac
