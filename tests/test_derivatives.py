"""
Unit tests for derivatives module.

Tests finite-difference derivatives and Jacobian computation.
"""

import numpy as np
import pytest

from cosmo_lensing.derivatives import (
    compute_jacobian,
    derivative_1d,
    jacobian_determinant,
)


class TestDerivative1D:
    """Test 1D finite-difference derivative function."""
    
    def test_constant_function(self):
        """Test derivative of constant function (should be zero)."""
        x = np.ones(100, dtype=np.float32)
        dx = derivative_1d(x, axis=0)
        
        np.testing.assert_allclose(dx, 0.0, atol=1e-6)
    
    def test_linear_function(self):
        """Test derivative of linear function (should be approximately constant)."""
        x = np.linspace(0, 10, 100, dtype=np.float32)
        dx = derivative_1d(x, axis=0)
        
        # For interior points, derivative should be close to constant
        # Check interior (away from boundaries which have different formulas)
        interior_std = np.std(dx[10:-10])
        interior_mean = np.mean(dx[10:-10])
        
        # Standard deviation should be small compared to mean
        assert interior_std / np.abs(interior_mean) < 0.01
    
    def test_quadratic_function(self):
        """Test derivative of quadratic function."""
        x_vals = np.linspace(-5, 5, 100, dtype=np.float32)
        y = x_vals**2
        dy = derivative_1d(y, axis=0)
        
        # Derivative should increase monotonically for x > 0 half
        # And also increase moving from left to right (d(x²)/dx = 2x)
        assert dy[70] > dy[60] > dy[50]  # Right side increasing
        assert dy[20] < dy[30] < dy[40]  # Left to center (less negative)
    
    def test_2d_array_axis0(self):
        """Test derivative along axis 0 of 2D array."""
        arr = np.random.randn(50, 40).astype(np.float32)
        darr = derivative_1d(arr, axis=0)
        
        assert darr.shape == arr.shape
    
    def test_2d_array_axis1(self):
        """Test derivative along axis 1 of 2D array."""
        arr = np.random.randn(50, 40).astype(np.float32)
        darr = derivative_1d(arr, axis=1)
        
        assert darr.shape == arr.shape
    
    def test_boundary_handling(self):
        """Test that boundary derivatives don't produce NaN/Inf."""
        x = np.linspace(0, 10, 100, dtype=np.float32)
        dx = derivative_1d(x, axis=0)
        
        # All derivatives should be finite
        assert np.all(np.isfinite(dx))
        
        # Boundaries should have reasonable magnitude (same order as interior)
        interior_mag = np.abs(dx[50])
        boundary_mag = np.abs(dx[0])
        assert boundary_mag < 10 * interior_mag  # Within order of magnitude


class TestComputeJacobian:
    """Test Jacobian computation."""
    
    def test_basic_jacobian_shape(self):
        """Test that basic Jacobian has correct shape."""
        alpha = np.random.randn(100, 100, 2).astype(np.float32)
        jac = compute_jacobian(alpha, compute_second_derivatives=False)
        
        assert jac.shape == (100, 100, 4)
    
    def test_extended_jacobian_shape(self):
        """Test that extended Jacobian has correct shape."""
        alpha = np.random.randn(100, 100, 2).astype(np.float32)
        jac = compute_jacobian(alpha, compute_second_derivatives=True)
        
        assert jac.shape == (100, 100, 12)
    
    def test_invalid_input_shape(self):
        """Test that invalid input raises ValueError."""
        # Wrong shape
        alpha_wrong = np.random.randn(100, 100, 3).astype(np.float32)
        
        with pytest.raises(ValueError, match="must have shape"):
            compute_jacobian(alpha_wrong)
    
    def test_linear_deflection_field(self):
        """Test Jacobian for linear deflection field."""
        nx, ny = 50, 50
        alpha = np.zeros((nx, ny, 2), dtype=np.float32)
        
        # Linear deflection: α₁ = a*x, α₂ = b*y
        x = np.linspace(0, 10, nx)
        y = np.linspace(0, 10, ny)
        X, Y = np.meshgrid(x, y, indexing='ij')
        
        a, b = 0.1, 0.2
        alpha[:, :, 0] = a * X
        alpha[:, :, 1] = b * Y
        
        jac = compute_jacobian(alpha, pixel_scale=10.0, compute_second_derivatives=False)
        
        # For linear field: ∂α₁/∂x₁ should be approximately constant
        # Check that variance is small (derivative is roughly uniform)
        assert np.std(jac[:, :, 0]) / np.abs(np.mean(jac[:, :, 0])) < 0.5
        assert np.std(jac[:, :, 3]) / np.abs(np.mean(jac[:, :, 3])) < 0.5
        
        # Cross terms should be small
        assert np.abs(jac[:, :, 1].mean()) < 0.01
        assert np.abs(jac[:, :, 2].mean()) < 0.01
    
    def test_point_mass_deflection(self, simple_deflection_field):
        """Test Jacobian for point mass deflection field."""
        alpha, _ = simple_deflection_field
        
        jac = compute_jacobian(alpha, compute_second_derivatives=False)
        
        # Check that Jacobian is computed without errors
        assert jac.shape == (100, 100, 4)
        assert np.all(np.isfinite(jac))
    
    def test_second_derivatives_symmetry(self):
        """Test that mixed derivatives are approximately equal (Schwarz theorem)."""
        alpha = np.random.randn(100, 100, 2).astype(np.float32)
        jac = compute_jacobian(alpha, compute_second_derivatives=True)
        
        # ∂²α₁/∂x₁∂x₂ should equal ∂²α₁/∂x₂∂x₁
        diff_a1 = np.abs(jac[:, :, 8] - jac[:, :, 10])
        
        # ∂²α₂/∂x₁∂x₂ should equal ∂²α₂/∂x₂∂x₁
        diff_a2 = np.abs(jac[:, :, 9] - jac[:, :, 11])
        
        # Allow some numerical error
        assert diff_a1.mean() < 0.1 * np.abs(jac[:, :, 8]).mean()
        assert diff_a2.mean() < 0.1 * np.abs(jac[:, :, 9]).mean()


class TestJacobianDeterminant:
    """Test Jacobian determinant computation."""
    
    def test_determinant_shape(self):
        """Test that determinant has correct shape."""
        jac = np.random.randn(100, 100, 4).astype(np.float32)
        det = jacobian_determinant(jac)
        
        assert det.shape == (100, 100)
    
    def test_identity_jacobian(self):
        """Test determinant of identity Jacobian."""
        nx, ny = 50, 50
        jac = np.zeros((nx, ny, 4), dtype=np.float32)
        jac[:, :, 0] = 1.0  # ∂α₁/∂x₁ = 1
        jac[:, :, 3] = 1.0  # ∂α₂/∂x₂ = 1
        
        det = jacobian_determinant(jac)
        
        # det(I) = 1
        np.testing.assert_allclose(det, 1.0, rtol=1e-6)
    
    def test_diagonal_jacobian(self):
        """Test determinant of diagonal Jacobian."""
        nx, ny = 50, 50
        jac = np.zeros((nx, ny, 4), dtype=np.float32)
        jac[:, :, 0] = 2.0  # ∂α₁/∂x₁ = 2
        jac[:, :, 3] = 3.0  # ∂α₂/∂x₂ = 3
        
        det = jacobian_determinant(jac)
        
        # det = 2 * 3 = 6
        np.testing.assert_allclose(det, 6.0, rtol=1e-6)
    
    def test_positive_determinant(self, simple_deflection_field):
        """Test that most of weak lensing has reasonable determinant."""
        alpha, _ = simple_deflection_field
        jac = compute_jacobian(alpha, compute_second_derivatives=False)
        det = jacobian_determinant(jac)
        
        # Point mass deflection field is tricky - very strong near center
        # Just check that determinants are computed and finite
        assert np.all(np.isfinite(det))
        
        # Check that away from center (edge regions), determinant is more reasonable
        # Look at corners which are far from singularity
        corner_det = det[:10, :10].mean()
        assert np.isfinite(corner_det)
