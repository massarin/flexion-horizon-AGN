"""
Lensing profile models for weak gravitational lensing.

Implements base class and specific profile models:
- NFW (Navarro-Frenk-White): Dark matter halos
- SIS (Singular Isothermal Sphere): Simple isothermal models

References:
    Wright & Brainerd (2000): "Gravitational lensing by NFW halos"
        ApJ 534, 34-40
    Kormann et al. (1994): "Isothermal elliptical gravitational lens models"
        A&A 284, 285-299
"""

import logging
from abc import ABC, abstractmethod
from typing import Tuple, Union

import numpy as np
from scipy.optimize import curve_fit

logger = logging.getLogger(__name__)


class LensingProfile(ABC):
    """
    Abstract base class for lensing profile models.

    All profile models must implement convergence, shear_tangential,
    and fit methods.
    """

    @abstractmethod
    def convergence(self, r: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute convergence κ(r).

        Args:
            r: Radial distances

        Returns:
            Convergence values
        """
        pass

    @abstractmethod
    def shear_tangential(self, r: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute tangential shear γ_t(r).

        Args:
            r: Radial distances

        Returns:
            Tangential shear values
        """
        pass


class NFWProfile(LensingProfile):
    """
    NFW halo profile for lensing calculations.

    The NFW profile is characterized by:
    - κ_s: convergence scale (dimensionless)
    - r_s: scale radius (arcsec or kpc)

    The surface density profile is:
        Σ(R) = Σ_crit * κ(x) * f(x)

    where x = R/r_s and f(x) is given in Wright & Brainerd (2000) Eq. 11.

    Attributes:
        ks: Convergence scale κ_s
        rs: Scale radius r_s

    References:
        Wright & Brainerd (2000) ApJ 534, 34-40
        Bartelmann (1996) A&A 313, 697-702
    """

    def __init__(self, ks: float, rs: float):
        """
        Initialize NFW profile.

        Args:
            ks: Convergence scale κ_s (dimensionless)
            rs: Scale radius r_s (same units as r in lensing calculations)
        """
        self.ks = ks
        self.rs = rs

        logger.info(f"Initialized NFW profile: κ_s={ks:.4f}, r_s={rs:.2f}")

    @staticmethod
    def _F_function(x: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Auxiliary function F(x) for NFW lensing.

        F(x) = { (1 - 2*arctanh(√((1-x)/(1+x)))/√(1-x²))/(x²-1)  for x < 1
               { 1/3                                              for x = 1
               { (1 - 2*arctan(√((x-1)/(1+x)))/√(x²-1))/(x²-1)   for x > 1

        Args:
            x: Dimensionless radius r/r_s

        Returns:
            F(x) values

        References:
            Wright & Brainerd (2000) Eq. 13
        """
        x = np.atleast_1d(np.asarray(x, dtype=np.float64))
        F = np.zeros_like(x)

        # x < 1 case
        mask_lt = x < 0.999
        if np.any(mask_lt):
            x_lt = x[mask_lt]
            arg = np.sqrt((1 - x_lt) / (1 + x_lt))
            F[mask_lt] = (1 - 2 * np.arctanh(arg) / np.sqrt(1 - x_lt**2)) / (x_lt**2 - 1)

        # x = 1 case
        mask_eq = (x >= 0.999) & (x <= 1.001)
        F[mask_eq] = 1.0 / 3.0

        # x > 1 case
        mask_gt = x > 1.001
        if np.any(mask_gt):
            x_gt = x[mask_gt]
            arg = np.sqrt((x_gt - 1) / (1 + x_gt))
            F[mask_gt] = (1 - 2 * np.arctan(arg) / np.sqrt(x_gt**2 - 1)) / (x_gt**2 - 1)

        return F.item() if F.size == 1 else F

    @staticmethod
    def _g_function(x: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Auxiliary function g(x) for mean convergence (shear calculation).

        g(x) = (1/x²) * ∫₀ˣ f(x') x' dx'

        where f(x) is the convergence profile shape.

        References:
            Wright & Brainerd (2000) Eq. 12-14
        """
        x = np.atleast_1d(np.asarray(x, dtype=np.float64))
        g = np.zeros_like(x)

        # x < 1 case
        mask_lt = x < 0.999
        if np.any(mask_lt):
            x_lt = x[mask_lt]
            arg = np.sqrt((1 - x_lt) / (1 + x_lt))
            ln_term = np.log(x_lt / 2)
            g[mask_lt] = (
                2 / (x_lt**2 - 1) * (1 - 2 / np.sqrt(1 - x_lt**2) * np.arctanh(arg))
                + ln_term / x_lt**2
            )

        # x = 1 case
        mask_eq = (x >= 0.999) & (x <= 1.001)
        g[mask_eq] = 1.0 + np.log(0.5)

        # x > 1 case
        mask_gt = x > 1.001
        if np.any(mask_gt):
            x_gt = x[mask_gt]
            arg = np.sqrt((x_gt - 1) / (1 + x_gt))
            ln_term = np.log(x_gt / 2)
            g[mask_gt] = (
                2 / (x_gt**2 - 1) * (1 - 2 / np.sqrt(x_gt**2 - 1) * np.arctan(arg))
                + ln_term / x_gt**2
            )

        return g.item() if g.size == 1 else g

    def convergence(self, r: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute NFW convergence profile κ(r).

        κ(r) = κ_s * f(x)

        where x = r/r_s and f(x) is given by F_function.

        Args:
            r: Radial distances (same units as r_s)

        Returns:
            Convergence values

        References:
            Wright & Brainerd (2000) Eq. 11
        """
        x = r / self.rs
        f_x = self._F_function(x)
        return 2 * self.ks * f_x

    def shear_tangential(self, r: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute NFW tangential shear profile γ_t(r).

        γ_t(r) = κ̄(<r) - κ(r)

        where κ̄(<r) is the mean convergence within radius r.

        Args:
            r: Radial distances (same units as r_s)

        Returns:
            Tangential shear values

        References:
            Wright & Brainerd (2000) Eq. 12
        """
        x = r / self.rs
        g_x = self._g_function(x)
        f_x = self._F_function(x)
        return self.ks * (2 * g_x - 2 * f_x)

    def excess_surface_density(self, r: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute excess surface density ΔΣ(r).

        ΔΣ(r) = Σ̄(<r) - Σ(r) = Σ_crit * γ_t(r)

        This is what's measured in weak lensing.

        Args:
            r: Radial distances

        Returns:
            Excess surface density (same units as Σ_crit)
        """
        return self.shear_tangential(r)

    @classmethod
    def from_mass_concentration(
        cls,
        mass: float,
        concentration: float,
        redshift_lens: float,
        redshift_source: float,
        cosmology=None,
    ) -> "NFWProfile":
        """
        Create NFW profile from mass and concentration.

        This requires cosmology to compute the critical surface density
        and convert M_200 to κ_s and r_s.

        Args:
            mass: Halo mass M_200 (Msun)
            concentration: Concentration parameter c_200
            redshift_lens: Lens redshift
            redshift_source: Source redshift
            cosmology: CosmologyCalculator instance

        Returns:
            NFWProfile instance

        Note:
            Full implementation requires cosmology module.
        """
        raise NotImplementedError(
            "from_mass_concentration() requires cosmology integration"
        )


class SISProfile(LensingProfile):
    """
    Singular Isothermal Sphere (SIS) profile for lensing.

    The SIS is characterized by a constant velocity dispersion σ_v,
    leading to a flat rotation curve and ρ ∝ r⁻².

    The lensing properties are:
        κ(r) = θ_E / (2r)
        γ_t(r) = θ_E / (2r)

    where θ_E is the Einstein radius:
        θ_E = 4π (σ_v/c)² (D_LS/D_S)

    For simplicity, we parameterize directly by θ_E.

    Attributes:
        theta_E: Einstein radius (same units as r)
        sigma_v: Velocity dispersion (km/s), if provided

    References:
        Kormann et al. (1994) A&A 284, 285-299
        Schneider, Ehlers & Falco (1992) "Gravitational Lenses"
    """

    def __init__(self, theta_E: float = None, sigma_v: float = None):
        """
        Initialize SIS profile.

        Must provide either theta_E or sigma_v (converted to theta_E).

        Args:
            theta_E: Einstein radius (arcsec or same units as r)
            sigma_v: Velocity dispersion (km/s), optional

        Raises:
            ValueError: If neither theta_E nor sigma_v provided
        """
        if theta_E is None and sigma_v is None:
            raise ValueError("Must provide either theta_E or sigma_v")

        if theta_E is not None:
            self.theta_E = theta_E
            self.sigma_v = sigma_v
        else:
            # For now, just store sigma_v
            # Full conversion requires cosmology and distances
            self.sigma_v = sigma_v
            self.theta_E = None
            raise NotImplementedError(
                "SIS from sigma_v requires cosmology integration"
            )

        logger.info(f"Initialized SIS profile: θ_E={self.theta_E:.4f}")

    def convergence(self, r: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute SIS convergence profile κ(r).

        κ(r) = θ_E / (2r)

        Args:
            r: Radial distances (same units as theta_E)

        Returns:
            Convergence values

        Note:
            Diverges as r → 0 (singular isothermal sphere).
            In practice, use with r > r_min to avoid singularity.
        """
        r = np.atleast_1d(np.asarray(r, dtype=np.float64))

        # Avoid division by zero
        with np.errstate(divide="ignore", invalid="ignore"):
            kappa = self.theta_E / (2 * r)

        # Handle r=0 case
        kappa[r == 0] = np.inf

        return kappa.item() if kappa.size == 1 else kappa

    def shear_tangential(self, r: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute SIS tangential shear profile γ_t(r).

        For SIS, γ_t(r) = κ(r) = θ_E / (2r)

        This special property (γ_t = κ) simplifies SIS modeling.

        Args:
            r: Radial distances (same units as theta_E)

        Returns:
            Tangential shear values
        """
        # For SIS, shear equals convergence
        return self.convergence(r)


def fit_nfw_profile(
    r: np.ndarray,
    xi: np.ndarray,
    xi_err: np.ndarray = None,
    initial_guess: Tuple[float, float] = None,
    observable: str = "shear",
) -> Tuple[NFWProfile, dict]:
    """
    Fit NFW profile to tangential shear or convergence data.

    Uses scipy.optimize.curve_fit to find best-fit (κ_s, r_s).

    Args:
        r: Radial bins (arcsec or kpc)
        xi: Measured tangential shear or convergence
        xi_err: Errors on xi (optional, assumed uniform if not provided)
        initial_guess: Initial (κ_s, r_s) values, default (0.1, r[len(r)//2])
        observable: Which observable to fit ('shear' or 'convergence')

    Returns:
        Tuple of (best_fit_profile, fit_info)
        - best_fit_profile: NFWProfile with best-fit parameters
        - fit_info: dict with params, errors, covariance, chi2, reduced_chi2, dof

    Raises:
        ValueError: If fitting fails or data is invalid
    """
    # Validate inputs
    if len(r) != len(xi):
        raise ValueError(f"r and xi must have same length: {len(r)} vs {len(xi)}")

    if len(r) < 3:
        raise ValueError("Need at least 3 data points to fit NFW (2 parameters + DOF)")

    # Set initial guess if not provided
    if initial_guess is None:
        # Reasonable defaults: moderate κ_s, scale radius near middle of data
        initial_guess = (0.1, r[len(r) // 2])

    logger.info(f"Fitting NFW to {len(r)} data points, initial guess: {initial_guess}")

    # Define model function
    def model(r_data, ks, rs):
        profile = NFWProfile(ks, rs)
        if observable == "shear":
            return profile.shear_tangential(r_data)
        elif observable == "convergence":
            return profile.convergence(r_data)
        else:
            raise ValueError(f"Unknown observable: {observable}")

    # Fit
    try:
        if xi_err is not None:
            # Weighted fit
            popt, pcov = curve_fit(
                model, r, xi, p0=initial_guess, sigma=xi_err, absolute_sigma=True
            )
        else:
            # Unweighted fit
            popt, pcov = curve_fit(model, r, xi, p0=initial_guess)

        # Extract results
        ks_fit, rs_fit = popt
        errors = np.sqrt(np.diag(pcov))

        # Compute chi-squared
        model_values = model(r, ks_fit, rs_fit)
        if xi_err is not None:
            residuals = (xi - model_values) / xi_err
        else:
            residuals = (xi - model_values) / np.std(xi)

        chi2 = np.sum(residuals**2)
        dof = len(r) - 2  # 2 parameters
        reduced_chi2 = chi2 / dof if dof > 0 else np.inf

        # Create best-fit profile
        best_profile = NFWProfile(ks_fit, rs_fit)

        # Assemble fit info
        fit_info = {
            "params": popt,
            "errors": errors,
            "covariance": pcov,
            "chi2": chi2,
            "reduced_chi2": reduced_chi2,
            "dof": dof,
            "success": True,
        }

        logger.info(
            f"NFW fit successful: κ_s={ks_fit:.4f}±{errors[0]:.4f}, "
            f"r_s={rs_fit:.2f}±{errors[1]:.2f}, χ²/dof={reduced_chi2:.2f}"
        )

        return best_profile, fit_info

    except Exception as e:
        logger.error(f"NFW fitting failed: {e}")
        raise ValueError(f"NFW profile fitting failed: {e}")


def fit_sis_profile(
    r: np.ndarray, xi: np.ndarray, xi_err: np.ndarray = None
) -> Tuple[SISProfile, dict]:
    """
    Fit SIS profile to tangential shear or convergence data.

    Since SIS has only one parameter (θ_E), fitting is simple.
    Uses weighted least squares: θ_E = 2 * Σ(xi * r * w) / Σ(w)
    where w = 1/xi_err² for weighted case.

    Args:
        r: Radial bins (arcsec or kpc)
        xi: Measured tangential shear or convergence
        xi_err: Errors on xi (optional, assumed uniform if not provided)

    Returns:
        Tuple of (best_fit_profile, fit_info)
        - best_fit_profile: SISProfile with best-fit θ_E
        - fit_info: dict with params, errors, chi2, reduced_chi2, dof

    Raises:
        ValueError: If data is invalid
    """
    # Validate inputs
    if len(r) != len(xi):
        raise ValueError(f"r and xi must have same length: {len(r)} vs {len(xi)}")

    if len(r) < 2:
        raise ValueError("Need at least 2 data points to fit SIS")

    logger.info(f"Fitting SIS to {len(r)} data points")

    # SIS model: γ_t(r) = θ_E / (2r)  →  θ_E = 2 * r * γ_t
    # Weighted least squares

    if xi_err is not None:
        weights = 1.0 / xi_err**2
    else:
        weights = np.ones_like(xi)

    # Fit: minimize Σ w * (xi - θ_E/(2r))²
    # Solution: θ_E = 2 * Σ(w * xi * r) / Σ(w)
    theta_E_fit = 2 * np.sum(weights * xi * r) / np.sum(weights)

    # Compute error on θ_E
    if xi_err is not None:
        # Error propagation: σ²(θ_E) = 4 * Σ(r² * σ²_xi) / n
        # For weighted fit: σ²(θ_E) = 1 / Σ(w)
        theta_E_err = np.sqrt(1.0 / np.sum(weights))
    else:
        # Unweighted: use residual variance
        model_values = theta_E_fit / (2 * r)
        residuals = xi - model_values
        sigma_res = np.std(residuals)
        theta_E_err = 2 * sigma_res * np.sqrt(np.sum(r**2)) / len(r)

    # Compute chi-squared
    model_values = theta_E_fit / (2 * r)
    if xi_err is not None:
        chi2 = np.sum(((xi - model_values) / xi_err) ** 2)
    else:
        chi2 = np.sum(((xi - model_values) / np.std(xi)) ** 2)

    dof = len(r) - 1  # 1 parameter
    reduced_chi2 = chi2 / dof if dof > 0 else np.inf

    # Create best-fit profile
    best_profile = SISProfile(theta_E=theta_E_fit)

    # Assemble fit info
    fit_info = {
        "params": np.array([theta_E_fit]),
        "errors": np.array([theta_E_err]),
        "chi2": chi2,
        "reduced_chi2": reduced_chi2,
        "dof": dof,
        "success": True,
    }

    logger.info(
        f"SIS fit successful: θ_E={theta_E_fit:.4f}±{theta_E_err:.4f}, "
        f"χ²/dof={reduced_chi2:.2f}"
    )

    return best_profile, fit_info
