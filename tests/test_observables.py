"""
Unit tests for observables module.

Tests lensing observable calculations (convergence, shear, flexion).
"""

import numpy as np
import pytest

from cosmo_lensing.observables import (
    compute_all_observables,
    convergence,
    flexion_F1,
    flexion_F2,
    flexion_G1,
    flexion_G2,
    rotation,
    shear_1,
    shear_2,
    shear_magnitude,
)


class TestConvergence:
    """Test convergence calculation."""
    
    def test_convergence_shape(self, test_jacobian):
        """Test that convergence has correct shape."""
        kappa = convergence(test_jacobian)
        
        assert kappa.shape == test_jacobian.shape[:2]
    
    def test_convergence_formula(self):
        """Test convergence formula κ = 0.5 * (∂α₁/∂x₁ + ∂α₂/∂x₂)."""
        nx, ny = 10, 10
        jac = np.zeros((nx, ny, 4), dtype=np.float32)
        jac[:, :, 0] = 0.1  # ∂α₁/∂x₁
        jac[:, :, 3] = 0.2  # ∂α₂/∂x₂
        
        kappa = convergence(jac)
        expected = 0.5 * (0.1 + 0.2)
        
        np.testing.assert_allclose(kappa, expected, rtol=1e-6)
    
    def test_convergence_from_extended_jacobian(self, extended_jacobian):
        """Test that convergence works with extended Jacobian."""
        kappa = convergence(extended_jacobian)
        
        assert kappa.shape == extended_jacobian.shape[:2]


class TestShear:
    """Test shear calculations."""
    
    def test_shear1_shape(self, test_jacobian):
        """Test that shear_1 has correct shape."""
        gamma1 = shear_1(test_jacobian)
        
        assert gamma1.shape == test_jacobian.shape[:2]
    
    def test_shear2_shape(self, test_jacobian):
        """Test that shear_2 has correct shape."""
        gamma2 = shear_2(test_jacobian)
        
        assert gamma2.shape == test_jacobian.shape[:2]
    
    def test_shear1_formula(self):
        """Test shear_1 formula γ₁ = 0.5 * (∂α₁/∂x₁ - ∂α₂/∂x₂)."""
        nx, ny = 10, 10
        jac = np.zeros((nx, ny, 4), dtype=np.float32)
        jac[:, :, 0] = 0.3  # ∂α₁/∂x₁
        jac[:, :, 3] = 0.1  # ∂α₂/∂x₂
        
        gamma1 = shear_1(jac)
        expected = 0.5 * (0.3 - 0.1)
        
        np.testing.assert_allclose(gamma1, expected, rtol=1e-6)
    
    def test_shear2_formula(self):
        """Test shear_2 formula γ₂ = 0.5 * (∂α₁/∂x₂ + ∂α₂/∂x₁)."""
        nx, ny = 10, 10
        jac = np.zeros((nx, ny, 4), dtype=np.float32)
        jac[:, :, 1] = 0.2  # ∂α₁/∂x₂
        jac[:, :, 2] = 0.4  # ∂α₂/∂x₁
        
        gamma2 = shear_2(jac)
        expected = 0.5 * (0.2 + 0.4)
        
        np.testing.assert_allclose(gamma2, expected, rtol=1e-6)
    
    def test_shear_magnitude(self):
        """Test shear magnitude calculation."""
        nx, ny = 10, 10
        jac = np.zeros((nx, ny, 4), dtype=np.float32)
        jac[:, :, 0] = 0.4  # γ₁ = 0.5*(0.4-0.0) = 0.2
        jac[:, :, 1] = 0.3  # γ₂ = 0.5*(0.3+0.3) = 0.3
        jac[:, :, 2] = 0.3
        
        gamma_mag = shear_magnitude(jac)
        gamma1 = shear_1(jac)
        gamma2 = shear_2(jac)
        expected = np.sqrt(gamma1**2 + gamma2**2)
        
        np.testing.assert_allclose(gamma_mag, expected, rtol=1e-6)


class TestFlexion:
    """Test flexion calculations."""
    
    def test_flexion_requires_extended_jacobian(self, test_jacobian):
        """Test that flexion raises error without second derivatives."""
        with pytest.raises(ValueError, match="Extended Jacobian required"):
            flexion_F1(test_jacobian)
    
    def test_flexion_F1_shape(self, extended_jacobian):
        """Test that flexion F1 has correct shape."""
        F1 = flexion_F1(extended_jacobian)
        
        assert F1.shape == extended_jacobian.shape[:2]
    
    def test_flexion_F2_shape(self, extended_jacobian):
        """Test that flexion F2 has correct shape."""
        F2 = flexion_F2(extended_jacobian)
        
        assert F2.shape == extended_jacobian.shape[:2]
    
    def test_flexion_G1_shape(self, extended_jacobian):
        """Test that flexion G1 has correct shape."""
        G1 = flexion_G1(extended_jacobian)
        
        assert G1.shape == extended_jacobian.shape[:2]
    
    def test_flexion_G2_shape(self, extended_jacobian):
        """Test that flexion G2 has correct shape."""
        G2 = flexion_G2(extended_jacobian)
        
        assert G2.shape == extended_jacobian.shape[:2]
    
    def test_flexion_F1_formula(self):
        """Test flexion F1 formula."""
        nx, ny = 10, 10
        jac = np.zeros((nx, ny, 12), dtype=np.float32)
        
        # F₁ = 0.5 * (∂²α₁/∂x₁² - ∂²α₂/∂x₁∂x₂) + 0.5 * (∂²α₂/∂x₂² + ∂²α₁/∂x₂∂x₁)
        jac[:, :, 4] = 0.1   # ∂²α₁/∂x₁²
        jac[:, :, 11] = 0.05  # ∂²α₂/∂x₁∂x₂ (index 11 in 0-based)
        jac[:, :, 6] = 0.2   # ∂²α₂/∂x₂²
        jac[:, :, 9] = 0.15  # ∂²α₁/∂x₂∂x₁ (index 9 in 0-based)
        
        F1 = flexion_F1(jac)
        expected = 0.5 * (0.1 - 0.05) + 0.5 * (0.2 + 0.15)
        
        np.testing.assert_allclose(F1, expected, rtol=1e-6)
    
    def test_flexion_F2_formula(self):
        """Test flexion F2 formula."""
        nx, ny = 10, 10
        jac = np.zeros((nx, ny, 12), dtype=np.float32)
        
        # F₂ = 0.5 * (∂²α₁/∂x₂² + ∂²α₂/∂x₂∂x₁) - 0.5 * (∂²α₁/∂x₁∂x₂ - ∂²α₂/∂x₁²)
        jac[:, :, 5] = 0.3   # ∂²α₁/∂x₂²
        jac[:, :, 10] = 0.1  # ∂²α₂/∂x₂∂x₁
        jac[:, :, 8] = 0.2   # ∂²α₁/∂x₁∂x₂
        jac[:, :, 7] = 0.15  # ∂²α₂/∂x₁²
        
        F2 = flexion_F2(jac)
        expected = 0.5 * (0.3 + 0.1) - 0.5 * (0.2 - 0.15)
        
        np.testing.assert_allclose(F2, expected, rtol=1e-6)


class TestRotation:
    """Test rotation calculation."""
    
    def test_rotation_shape(self, test_jacobian):
        """Test that rotation has correct shape."""
        omega = rotation(test_jacobian)
        
        assert omega.shape == test_jacobian.shape[:2]
    
    def test_rotation_formula(self):
        """Test rotation formula ω = 0.5 * (∂α₁/∂x₂ - ∂α₂/∂x₁)."""
        nx, ny = 10, 10
        jac = np.zeros((nx, ny, 4), dtype=np.float32)
        jac[:, :, 1] = 0.1  # ∂α₁/∂x₂
        jac[:, :, 2] = 0.3  # ∂α₂/∂x₁
        
        omega = rotation(jac)
        expected = 0.5 * (0.1 - 0.3)
        
        np.testing.assert_allclose(omega, expected, rtol=1e-6)
    
    def test_rotation_zero_for_curl_free(self):
        """Test that rotation is zero for curl-free (lensing) field."""
        nx, ny = 50, 50
        jac = np.zeros((nx, ny, 4), dtype=np.float32)
        
        # For curl-free field: ∂α₁/∂x₂ = ∂α₂/∂x₁
        jac[:, :, 1] = 0.2
        jac[:, :, 2] = 0.2
        
        omega = rotation(jac)
        
        np.testing.assert_allclose(omega, 0.0, atol=1e-6)


class TestComputeAllObservables:
    """Test computing all observables at once."""
    
    def test_basic_observables(self, test_jacobian):
        """Test computing all observables from basic Jacobian."""
        obs = compute_all_observables(test_jacobian)
        
        # Should have basic observables
        assert 'kappa' in obs
        assert 'gamma1' in obs
        assert 'gamma2' in obs
        assert 'gamma' in obs
        assert 'rot' in obs
        
        # Should not have flexion
        assert 'F1' not in obs
        assert 'F2' not in obs
        assert 'G1' not in obs
        assert 'G2' not in obs
        
        # Check shapes
        assert obs['kappa'].shape == test_jacobian.shape[:2]
        assert obs['gamma1'].shape == test_jacobian.shape[:2]
    
    def test_extended_observables(self, extended_jacobian):
        """Test computing all observables including flexion."""
        obs = compute_all_observables(extended_jacobian)
        
        # Should have all observables
        assert 'kappa' in obs
        assert 'gamma1' in obs
        assert 'gamma2' in obs
        assert 'gamma' in obs
        assert 'rot' in obs
        assert 'F1' in obs
        assert 'F2' in obs
        assert 'G1' in obs
        assert 'G2' in obs
        
        # Check shapes
        assert obs['F1'].shape == extended_jacobian.shape[:2]
        assert obs['G2'].shape == extended_jacobian.shape[:2]
    
    def test_observable_consistency(self, test_jacobian):
        """Test that individual functions match compute_all_observables."""
        obs = compute_all_observables(test_jacobian)
        
        np.testing.assert_array_equal(obs['kappa'], convergence(test_jacobian))
        np.testing.assert_array_equal(obs['gamma1'], shear_1(test_jacobian))
        np.testing.assert_array_equal(obs['gamma2'], shear_2(test_jacobian))
        np.testing.assert_array_equal(obs['rot'], rotation(test_jacobian))
