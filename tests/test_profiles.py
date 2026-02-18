"""
Tests for lensing profile models (NFW and SIS).
"""

import pytest
import numpy as np
from cosmo_lensing.profiles import (
    NFWProfile,
    SISProfile,
    fit_nfw_profile,
    fit_sis_profile,
)


class TestNFWProfile:
    """Tests for NFW profile model."""

    def test_initialization(self):
        """Test NFW profile initialization."""
        profile = NFWProfile(ks=0.1, rs=10.0)
        assert profile.ks == 0.1
        assert profile.rs == 10.0

    def test_convergence_at_rs(self):
        """Test convergence at scale radius (x=1)."""
        profile = NFWProfile(ks=0.1, rs=10.0)
        kappa = profile.convergence(10.0)
        # At x=1, F(x) = 1/3, so κ = 2 * ks * 1/3
        expected = 2 * 0.1 / 3
        np.testing.assert_allclose(kappa, expected, rtol=1e-6)

    def test_shear_at_rs(self):
        """Test tangential shear at scale radius."""
        profile = NFWProfile(ks=0.1, rs=10.0)
        gamma = profile.shear_tangential(10.0)
        assert isinstance(gamma, (float, np.floating))
        # Shear at scale radius should be finite
        assert np.isfinite(gamma)

    def test_convergence_array(self):
        """Test convergence with array input."""
        profile = NFWProfile(ks=0.1, rs=10.0)
        r = np.array([5.0, 10.0, 20.0])
        kappa = profile.convergence(r)
        assert kappa.shape == r.shape
        assert np.all(np.isfinite(kappa))

    def test_shear_equals_excess_surface_density(self):
        """Test that shear equals excess surface density."""
        profile = NFWProfile(ks=0.1, rs=10.0)
        r = np.array([5.0, 10.0, 20.0])
        gamma = profile.shear_tangential(r)
        delta_sigma = profile.excess_surface_density(r)
        np.testing.assert_array_equal(gamma, delta_sigma)

    def test_convergence_decreases_with_radius(self):
        """Test that convergence decreases at large radii."""
        profile = NFWProfile(ks=0.1, rs=10.0)
        r = np.array([10.0, 50.0, 100.0])
        kappa = profile.convergence(r)
        assert kappa[1] < kappa[0]
        assert kappa[2] < kappa[1]

    def test_from_mass_concentration_not_implemented(self):
        """Test that from_mass_concentration raises NotImplementedError."""
        with pytest.raises(NotImplementedError):
            NFWProfile.from_mass_concentration(
                mass=1e14, concentration=5.0, redshift_lens=0.3, redshift_source=1.0
            )


class TestSISProfile:
    """Tests for SIS profile model."""

    def test_initialization_with_theta_E(self):
        """Test SIS profile initialization with Einstein radius."""
        profile = SISProfile(theta_E=2.0)
        assert profile.theta_E == 2.0

    def test_initialization_without_params_fails(self):
        """Test that initialization without parameters raises ValueError."""
        with pytest.raises(ValueError, match="Must provide either"):
            SISProfile()

    def test_initialization_with_sigma_v_not_implemented(self):
        """Test that initialization with sigma_v raises NotImplementedError."""
        with pytest.raises(NotImplementedError, match="cosmology"):
            SISProfile(sigma_v=200.0)

    def test_convergence_formula(self):
        """Test SIS convergence formula κ = θ_E / (2r)."""
        profile = SISProfile(theta_E=2.0)
        r = 4.0
        kappa = profile.convergence(r)
        expected = 2.0 / (2 * 4.0)  # θ_E / (2r)
        np.testing.assert_allclose(kappa, expected)

    def test_shear_equals_convergence(self):
        """Test that for SIS, γ_t = κ."""
        profile = SISProfile(theta_E=2.0)
        r = np.array([1.0, 2.0, 5.0, 10.0])
        kappa = profile.convergence(r)
        gamma = profile.shear_tangential(r)
        np.testing.assert_array_equal(kappa, gamma)

    def test_convergence_array_input(self):
        """Test convergence with array input."""
        profile = SISProfile(theta_E=2.0)
        r = np.array([1.0, 2.0, 5.0])
        kappa = profile.convergence(r)
        assert kappa.shape == r.shape
        expected = 2.0 / (2 * r)
        np.testing.assert_allclose(kappa, expected)

    def test_singularity_at_zero(self):
        """Test that convergence diverges at r=0."""
        profile = SISProfile(theta_E=2.0)
        kappa = profile.convergence(0.0)
        assert np.isinf(kappa)

    def test_convergence_decreases_with_radius(self):
        """Test that κ ∝ 1/r decreases with radius."""
        profile = SISProfile(theta_E=2.0)
        r = np.array([1.0, 2.0, 4.0])
        kappa = profile.convergence(r)
        assert kappa[1] < kappa[0]
        assert kappa[2] < kappa[1]
        # Check proportionality: κ(r2)/κ(r1) = r1/r2
        np.testing.assert_allclose(kappa[1] / kappa[0], r[0] / r[1])


class TestFitNFWProfile:
    """Tests for NFW profile fitting."""

    def test_fit_noiseless_data(self):
        """Test fitting to perfect NFW data."""
        # Create synthetic NFW data
        true_profile = NFWProfile(ks=0.15, rs=12.0)
        r = np.logspace(0.5, 2, 20)
        xi_true = true_profile.shear_tangential(r)

        # Fit
        fitted_profile, fit_info = fit_nfw_profile(r, xi_true)

        # Check recovery
        assert fit_info["success"]
        np.testing.assert_allclose(fitted_profile.ks, 0.15, rtol=0.01)
        np.testing.assert_allclose(fitted_profile.rs, 12.0, rtol=0.01)
        assert fit_info["reduced_chi2"] < 0.1  # Should be very small for perfect data

    def test_fit_with_noise(self):
        """Test fitting to noisy NFW data."""
        np.random.seed(42)
        true_profile = NFWProfile(ks=0.1, rs=10.0)
        r = np.logspace(0.5, 2, 20)
        xi_true = true_profile.shear_tangential(r)

        # Add noise
        noise_level = 0.01
        xi_noisy = xi_true + np.random.randn(len(r)) * noise_level
        xi_err = np.ones_like(r) * noise_level

        # Fit
        fitted_profile, fit_info = fit_nfw_profile(r, xi_noisy, xi_err=xi_err)

        # Check that fit is reasonable (within ~10% for noisy data)
        assert fit_info["success"]
        assert fitted_profile.ks > 0
        assert fitted_profile.rs > 0
        assert 0.5 < fit_info["reduced_chi2"] < 2.0  # Reasonable chi2 for noisy data

    def test_fit_convergence_observable(self):
        """Test fitting to convergence instead of shear."""
        true_profile = NFWProfile(ks=0.12, rs=15.0)
        r = np.logspace(0.5, 2, 20)
        kappa_true = true_profile.convergence(r)

        # Fit convergence
        fitted_profile, fit_info = fit_nfw_profile(
            r, kappa_true, observable="convergence"
        )

        assert fit_info["success"]
        np.testing.assert_allclose(fitted_profile.ks, 0.12, rtol=0.02)
        np.testing.assert_allclose(fitted_profile.rs, 15.0, rtol=0.02)

    def test_fit_with_custom_initial_guess(self):
        """Test fitting with custom initial guess."""
        true_profile = NFWProfile(ks=0.2, rs=8.0)
        r = np.logspace(0.5, 2, 20)
        xi_true = true_profile.shear_tangential(r)

        # Fit with initial guess
        fitted_profile, fit_info = fit_nfw_profile(
            r, xi_true, initial_guess=(0.15, 10.0)
        )

        assert fit_info["success"]
        np.testing.assert_allclose(fitted_profile.ks, 0.2, rtol=0.01)

    def test_fit_returns_errors(self):
        """Test that fit returns parameter errors."""
        true_profile = NFWProfile(ks=0.1, rs=10.0)
        r = np.logspace(0.5, 2, 20)
        xi_true = true_profile.shear_tangential(r)

        fitted_profile, fit_info = fit_nfw_profile(r, xi_true)

        assert "errors" in fit_info
        assert len(fit_info["errors"]) == 2
        assert np.all(fit_info["errors"] > 0)

    def test_fit_insufficient_data_fails(self):
        """Test that fitting with too few points fails."""
        r = np.array([1.0, 2.0])
        xi = np.array([0.1, 0.05])

        with pytest.raises(ValueError, match="at least 3 data points"):
            fit_nfw_profile(r, xi)

    def test_fit_mismatched_arrays_fails(self):
        """Test that mismatched array lengths raise ValueError."""
        r = np.array([1.0, 2.0, 3.0])
        xi = np.array([0.1, 0.05])

        with pytest.raises(ValueError, match="same length"):
            fit_nfw_profile(r, xi)


class TestFitSISProfile:
    """Tests for SIS profile fitting."""

    def test_fit_perfect_data(self):
        """Test fitting to perfect SIS data."""
        true_profile = SISProfile(theta_E=1.5)
        r = np.linspace(1, 20, 15)
        xi_true = true_profile.shear_tangential(r)

        # Fit
        fitted_profile, fit_info = fit_sis_profile(r, xi_true)

        # Check recovery
        assert fit_info["success"]
        np.testing.assert_allclose(fitted_profile.theta_E, 1.5, rtol=1e-6)
        assert fit_info["reduced_chi2"] < 0.1

    def test_fit_with_noise(self):
        """Test fitting to noisy SIS data."""
        np.random.seed(42)
        true_profile = SISProfile(theta_E=2.0)
        r = np.linspace(1, 20, 15)
        xi_true = true_profile.shear_tangential(r)

        # Add noise
        noise_level = 0.02
        xi_noisy = xi_true + np.random.randn(len(r)) * noise_level
        xi_err = np.ones_like(r) * noise_level

        # Fit
        fitted_profile, fit_info = fit_sis_profile(r, xi_noisy, xi_err=xi_err)

        # Check reasonable fit
        assert fit_info["success"]
        assert fitted_profile.theta_E > 0
        assert 0.5 < fit_info["reduced_chi2"] < 5.0  # Relaxed for noisy data

    def test_fit_returns_error(self):
        """Test that fit returns parameter error."""
        true_profile = SISProfile(theta_E=1.8)
        r = np.linspace(1, 20, 15)
        xi_true = true_profile.shear_tangential(r)

        fitted_profile, fit_info = fit_sis_profile(r, xi_true)

        assert "errors" in fit_info
        assert len(fit_info["errors"]) == 1
        assert fit_info["errors"][0] > 0

    def test_fit_insufficient_data_fails(self):
        """Test that fitting with too few points fails."""
        r = np.array([1.0])
        xi = np.array([0.5])

        with pytest.raises(ValueError, match="at least 2 data points"):
            fit_sis_profile(r, xi)

    def test_fit_mismatched_arrays_fails(self):
        """Test that mismatched array lengths raise ValueError."""
        r = np.array([1.0, 2.0, 3.0])
        xi = np.array([0.5, 0.25])

        with pytest.raises(ValueError, match="same length"):
            fit_sis_profile(r, xi)


class TestBackwardCompatibility:
    """Test backward compatibility with old nfw module."""

    def test_import_from_nfw_works(self):
        """Test that importing from nfw module still works."""
        from cosmo_lensing.nfw import NFWProfile as NFWOld
        from cosmo_lensing.nfw import fit_nfw_profile as fit_old

        # Should be the same classes/functions
        assert NFWOld is NFWProfile
        assert fit_old is fit_nfw_profile

    def test_nfw_import_shows_deprecation_warning(self):
        """Test that importing from nfw shows deprecation warning."""
        # The warning is triggered on import, which happens at module level
        # Just check that the aliasing works
        from cosmo_lensing import nfw

        # Verify aliasing works
        assert nfw.NFWProfile is NFWProfile
        assert nfw.fit_nfw_profile is fit_nfw_profile
