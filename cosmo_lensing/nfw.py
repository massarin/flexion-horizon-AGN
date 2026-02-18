"""
NFW profile modeling for weak lensing.

Implements NFW (Navarro-Frenk-White) halo profile calculations
for convergence and shear.

References:
    Wright & Brainerd (2000): "Gravitational lensing by NFW halos"
        ApJ 534, 34-40
    Bartelmann (1996): "Analytical lensing calculations for NFW halos"
        A&A 313, 697-702
"""

import logging
from typing import Tuple, Union

import numpy as np

logger = logging.getLogger(__name__)


class NFWProfile:
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
            This is a placeholder for Phase 4.
        """
        # TODO: Implement in Phase 4 when integrating with cosmology
        raise NotImplementedError(
            "from_mass_concentration() requires cosmology integration - Phase 4"
        )


def fit_nfw_profile(
    r: np.ndarray,
    xi: np.ndarray,
    xi_err: np.ndarray = None,
    initial_guess: Tuple[float, float] = None,
) -> Tuple[NFWProfile, dict]:
    """
    Fit NFW profile to tangential shear/convergence data.

    Uses scipy.optimize.curve_fit to find best-fit (κ_s, r_s).

    Args:
        r: Radial bins (arcsec or kpc)
        xi: Measured tangential shear or convergence
        xi_err: Errors on xi (optional)
        initial_guess: Initial (κ_s, r_s) values

    Returns:
        Tuple of (best_fit_profile, fit_info)
        - best_fit_profile: NFWProfile with best-fit parameters
        - fit_info: dict with chi2, reduced_chi2, covariance, etc.

    Note:
        Full implementation in Phase 4 with correlation module.
    """
    # TODO: Implement in Phase 4 with scipy.optimize.curve_fit
    raise NotImplementedError("fit_nfw_profile() - Phase 4")
