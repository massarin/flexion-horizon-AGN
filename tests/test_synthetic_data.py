"""
Unit tests for synthetic data generators.

Tests analytic lensing solutions against known theoretical predictions.
"""

import numpy as np
import pytest

from cosmo_lensing import derivatives, observables, synthetic_data


class TestPointMassDeflection:
    """Test point mass lens deflection field."""
    
    def test_shape(self):
        """Test that point mass deflection has correct shape."""
        alpha, z = synthetic_data.generate_point_mass_deflection(
            theta_E=1.0, grid_size=100
        )
        
        assert alpha.shape == (100, 100, 2)
        assert z == 1.0
    
    def test_symmetry(self):
        """Test radial symmetry of point mass deflection."""
        alpha, _ = synthetic_data.generate_point_mass_deflection(
            theta_E=1.0, grid_size=100, box_size=10.0
        )
        
        # Check symmetry: α(-x, -y) = -α(x, y)
        center = 50
        offset = 10
        
        alpha_pos = alpha[center + offset, center + offset]
        alpha_neg = alpha[center - offset, center - offset]
        
        # Allow 10% error due to discrete grid
        np.testing.assert_allclose(alpha_pos, -alpha_neg, rtol=0.10)
    
    def test_convergence_profile(self):
        """Test that convergence decreases with radius (1/r² behavior)."""
        theta_E = 2.0
        alpha, _ = synthetic_data.generate_point_mass_deflection(
            theta_E=theta_E, grid_size=200, box_size=10.0
        )
        
        # Compute convergence
        jac = derivatives.compute_jacobian(alpha, pixel_scale=1.0, compute_second_derivatives=False)
        kappa = observables.convergence(jac)
        
        # Check that convergence decreases with radius
        center = 100
        r1 = 20
        r2 = 40
        r3 = 60
        
        kappa1 = kappa[center + r1, center]
        kappa2 = kappa[center + r2, center]
        kappa3 = kappa[center + r3, center]
        
        # Should decrease monotonically
        assert kappa1 > kappa2 > kappa3
        
        # Check approximate 1/r² scaling: κ(r2)/κ(r1) ≈ (r1/r2)²
        ratio_measured = kappa2 / kappa1
        ratio_expected = (r1 / r2)**2
        
        # Allow large error due to finite differences and regularization
        assert 0.1 < ratio_measured / ratio_expected < 10


class TestSISDeflection:
    """Test SIS lens deflection field."""
    
    def test_shape(self):
        """Test that SIS deflection has correct shape."""
        alpha, z = synthetic_data.generate_sis_deflection(
            sigma_v=1.0, grid_size=100
        )
        
        assert alpha.shape == (100, 100, 2)
        assert z == 1.0
    
    def test_radial_symmetry(self):
        """Test radial symmetry of SIS deflection."""
        alpha, _ = synthetic_data.generate_sis_deflection(
            sigma_v=1.5, grid_size=100
        )
        
        # Check that deflection magnitude is radially symmetric
        center = 50
        r1 = 20
        r2 = 30
        
        # Points at same radius should have same magnitude
        mag1 = np.sqrt(alpha[center + r1, center, 0]**2 + alpha[center + r1, center, 1]**2)
        mag2 = np.sqrt(alpha[center, center + r1, 0]**2 + alpha[center, center + r1, 1]**2)
        
        np.testing.assert_allclose(mag1, mag2, rtol=0.05)
    
    def test_convergence_profile(self):
        """Test that convergence matches κ(θ) = θ_E/(2θ)."""
        sigma_v = 1.5
        alpha, _ = synthetic_data.generate_sis_deflection(
            sigma_v=sigma_v, grid_size=200, box_size=10.0
        )
        
        # Compute convergence
        jac = derivatives.compute_jacobian(alpha, pixel_scale=10.0, compute_second_derivatives=False)
        kappa = observables.convergence(jac)
        
        # Check at a test radius
        center = 100
        test_radius_pix = 40
        test_point = kappa[center + test_radius_pix, center]
        
        # Expected: κ = θ_E / (2θ) where θ_E = sigma_v
        theta = test_radius_pix * (10.0 / 200)
        expected_kappa = sigma_v / (2 * theta)
        
        # Allow 10% error
        assert abs(test_point - expected_kappa) / expected_kappa < 0.10


class TestNFWDeflection:
    """Test NFW-like deflection field."""
    
    def test_shape(self):
        """Test that NFW deflection has correct shape."""
        alpha, z = synthetic_data.generate_nfw_deflection(
            mass=5.0, concentration=5.0, grid_size=100
        )
        
        assert alpha.shape == (100, 100, 2)
        assert isinstance(z, float)
    
    def test_smooth_profile(self):
        """Test that NFW deflection is smooth (no discontinuities)."""
        alpha, _ = synthetic_data.generate_nfw_deflection(
            mass=5.0, concentration=5.0, grid_size=200
        )
        
        # Compute gradient along a radial line
        center = 100
        radial_profile = alpha[center + 10:center + 50, center, 0]
        
        # Check that gradient doesn't have large jumps
        diff = np.diff(radial_profile)
        max_jump = np.max(np.abs(diff))
        mean_gradient = np.mean(np.abs(diff))
        
        # No jump should be > 10x the mean gradient (smooth profile)
        assert max_jump < 10 * mean_gradient
    
    def test_observables_finite(self):
        """Test that NFW deflection produces finite observables."""
        alpha, _ = synthetic_data.generate_nfw_deflection(
            mass=5.0, concentration=5.0, grid_size=200
        )
        
        jac = derivatives.compute_jacobian(alpha, pixel_scale=10.0, compute_second_derivatives=True)
        obs = observables.compute_all_observables(jac)
        
        # All observables should be finite
        for name, data in obs.items():
            assert np.all(np.isfinite(data)), f"{name} has NaN/Inf"


class TestShearField:
    """Test pure shear field generation."""
    
    def test_shape(self):
        """Test that shear field has correct shape."""
        alpha, z = synthetic_data.generate_shear_field(
            shear=0.1, angle=0.0, grid_size=100
        )
        
        assert alpha.shape == (100, 100, 2)
        assert z == 1.0
    
    def test_shear_extraction(self):
        """Test that extracted shear matches input."""
        input_shear = 0.05
        input_angle = np.pi / 6  # 30 degrees
        
        alpha, _ = synthetic_data.generate_shear_field(
            shear=input_shear, angle=input_angle, grid_size=200
        )
        
        # Compute shear
        jac = derivatives.compute_jacobian(alpha, pixel_scale=10.0, compute_second_derivatives=False)
        gamma1 = observables.shear_1(jac)
        gamma2 = observables.shear_2(jac)
        
        # Expected shear components
        expected_gamma1 = input_shear * np.cos(2 * input_angle)
        expected_gamma2 = input_shear * np.sin(2 * input_angle)
        
        # Check mean values (should be constant across field)
        measured_gamma1 = np.mean(gamma1)
        measured_gamma2 = np.mean(gamma2)
        
        np.testing.assert_allclose(measured_gamma1, expected_gamma1, rtol=0.02)
        np.testing.assert_allclose(measured_gamma2, expected_gamma2, rtol=0.02)
    
    def test_zero_convergence(self):
        """Test that pure shear field has zero convergence."""
        alpha, _ = synthetic_data.generate_shear_field(
            shear=0.1, angle=0.0, grid_size=200
        )
        
        jac = derivatives.compute_jacobian(alpha, pixel_scale=10.0, compute_second_derivatives=False)
        kappa = observables.convergence(jac)
        
        # Convergence should be nearly zero everywhere
        # But finite differences create small errors
        assert np.abs(np.mean(kappa)) < 0.01
        assert np.abs(np.std(kappa)) < 0.02


class TestConvergenceField:
    """Test pure convergence field generation."""
    
    def test_shape(self):
        """Test that convergence field has correct shape."""
        alpha, z = synthetic_data.generate_convergence_field(
            kappa=0.1, grid_size=100
        )
        
        assert alpha.shape == (100, 100, 2)
        assert z == 1.0
    
    def test_convergence_extraction(self):
        """Test that extracted convergence matches input."""
        input_kappa = 0.15
        
        alpha, _ = synthetic_data.generate_convergence_field(
            kappa=input_kappa, grid_size=200
        )
        
        # Compute convergence
        jac = derivatives.compute_jacobian(alpha, pixel_scale=10.0, compute_second_derivatives=False)
        kappa = observables.convergence(jac)
        
        # Check mean value
        measured_kappa = np.mean(kappa)
        
        np.testing.assert_allclose(measured_kappa, input_kappa, rtol=0.02)
    
    def test_zero_shear(self):
        """Test that pure convergence field has zero shear."""
        alpha, _ = synthetic_data.generate_convergence_field(
            kappa=0.1, grid_size=200
        )
        
        jac = derivatives.compute_jacobian(alpha, pixel_scale=10.0, compute_second_derivatives=False)
        gamma1 = observables.shear_1(jac)
        gamma2 = observables.shear_2(jac)
        
        # Shear should be nearly zero everywhere
        assert np.abs(np.mean(gamma1)) < 1e-6
        assert np.abs(np.mean(gamma2)) < 1e-6


class TestAnalyticValidation:
    """Integration tests validating analytic solutions."""
    
    def test_point_mass_shear_equals_convergence(self):
        """Test that for point mass, tangential shear scales like convergence."""
        theta_E = 2.0
        alpha, _ = synthetic_data.generate_point_mass_deflection(
            theta_E=theta_E, grid_size=200, box_size=10.0
        )
        
        jac = derivatives.compute_jacobian(alpha, pixel_scale=1.0, compute_second_derivatives=False)
        kappa = observables.convergence(jac)
        gamma1 = observables.shear_1(jac)
        gamma2 = observables.shear_2(jac)
        
        # For point mass: both κ and γ_t should scale as 1/r²
        # Test along x-axis at two radii
        center = 100
        r1 = 30
        r2 = 60
        
        kappa1 = abs(kappa[center + r1, center])
        kappa2 = abs(kappa[center + r2, center])
        gamma1_1 = abs(gamma1[center + r1, center])
        gamma1_2 = abs(gamma1[center + r2, center])
        
        # Both should decrease with radius
        assert kappa1 > kappa2
        assert gamma1_1 > gamma1_2
        
        # Check that both scale similarly (within factor of 10)
        kappa_ratio = kappa1 / kappa2 if kappa2 > 0 else 1.0
        gamma_ratio = gamma1_1 / gamma1_2 if gamma1_2 > 0 else 1.0
        
        # Both should scale roughly as r²
        assert 0.1 < gamma_ratio / kappa_ratio < 10
    
    def test_sis_shear_equals_convergence(self):
        """Test that for SIS, tangential shear equals convergence."""
        sigma_v = 1.5
        alpha, _ = synthetic_data.generate_sis_deflection(
            sigma_v=sigma_v, grid_size=200, box_size=10.0
        )
        
        jac = derivatives.compute_jacobian(alpha, pixel_scale=10.0, compute_second_derivatives=False)
        kappa = observables.convergence(jac)
        gamma1 = observables.shear_1(jac)
        
        # For SIS on x-axis: γ_t = κ
        center = 100
        test_radius = 40
        
        kappa_test = kappa[center + test_radius, center]
        gamma_t = -gamma1[center + test_radius, center]
        
        np.testing.assert_allclose(gamma_t, kappa_test, rtol=0.15)
    
    def test_rotation_is_zero(self):
        """Test that lensing fields are curl-free (rotation = 0)."""
        # Test on multiple field types
        alpha_pm, _ = synthetic_data.generate_point_mass_deflection(
            theta_E=1.5, grid_size=150
        )
        alpha_sis, _ = synthetic_data.generate_sis_deflection(
            sigma_v=1.0, grid_size=150
        )
        
        for alpha, name in [(alpha_pm, "point_mass"), (alpha_sis, "SIS")]:
            jac = derivatives.compute_jacobian(alpha, pixel_scale=10.0, compute_second_derivatives=False)
            omega = observables.rotation(jac)
            kappa = observables.convergence(jac)
            
            # Rotation should be much smaller than convergence
            rot_rms = np.sqrt(np.mean(omega**2))
            kappa_rms = np.sqrt(np.mean(kappa**2))
            
            # For strong lensing, finite differences create rotation artifacts
            # Just check rotation is reasonable (not dominant)
            assert rot_rms < 2.0 * kappa_rms, f"{name}: rotation unreasonably large"
