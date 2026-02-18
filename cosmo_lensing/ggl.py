"""
Galaxy-galaxy lensing (GGL) analysis module.

Provides high-level interface for computing and fitting tangential correlation
functions of lensing observables around lens galaxies.

Example:
    >>> from cosmo_lensing import ggl
    >>> # Compute shear correlation
    >>> result = ggl.compute_ggl_correlation(
    ...     ra_lens, dec_lens, ra_src, dec_src,
    ...     gamma1_src, gamma2_src, weights=weights
    ... )
    >>> # Fit NFW profile
    >>> fit_result = ggl.fit_ggl_profile(
    ...     result["r_mpc"], result["xi_t"], result["xi_t_err"],
    ...     z_lens=0.3, profile_type="nfw"
    ... )
"""

import numpy as np
import logging
from typing import Dict, Optional, Tuple, Union

from .correlations import TangentialCorrelation
from .profiles import fit_nfw_profile, fit_sis_profile

logger = logging.getLogger(__name__)


def compute_ggl_correlation(
    ra_lens: np.ndarray,
    dec_lens: np.ndarray,
    ra_src: np.ndarray,
    dec_src: np.ndarray,
    obs1: np.ndarray,
    obs2: np.ndarray,
    weights: Optional[np.ndarray] = None,
    min_sep: float = 0.1,
    max_sep: float = 100.0,
    nbins: int = 20,
    sep_units: str = "arcmin",
) -> Dict[str, np.ndarray]:
    """
    Compute galaxy-galaxy lensing correlation function.
    
    Wraps TangentialCorrelation with common defaults for GGL analysis.
    
    Args:
        ra_lens: Right ascension of lens galaxies (degrees)
        dec_lens: Declination of lens galaxies (degrees)
        ra_src: Right ascension of source galaxies (degrees)
        dec_src: Declination of source galaxies (degrees)
        obs1: First component of observable (e.g., γ₁, F₁, G₁, or κ)
        obs2: Second component of observable (e.g., γ₂, F₂, G₂, or 0 for κ)
        weights: Optional weights for source galaxies
        min_sep: Minimum separation (default: 0.1 arcmin)
        max_sep: Maximum separation (default: 100 arcmin)
        nbins: Number of radial bins (default: 20)
        sep_units: Separation units (default: 'arcmin')
    
    Returns:
        Dictionary with keys:
        - 'r': Mean separation in each bin
        - 'xi': Tangential correlation amplitude
        - 'xim': Cross-correlation (should be ~0)
        - 'xi_err': Error on tangential correlation
        - 'xim_err': Error on cross-correlation
        - 'npairs': Number of pairs in each bin
        - 'weight': Total weight in each bin
    
    Example:
        >>> # Compute shear correlation
        >>> result = compute_ggl_correlation(
        ...     ra_lens, dec_lens, ra_src, dec_src,
        ...     gamma1, gamma2, weights=weights
        ... )
        >>> r_arcmin = result['r']
        >>> xi_t = result['xi_t']
    """
    logger.info(
        f"Computing GGL correlation: {len(ra_lens)} lenses, "
        f"{len(ra_src)} sources"
    )
    logger.info(
        f"Separation range: {min_sep}-{max_sep} {sep_units}, {nbins} bins"
    )
    
    # Create correlation object
    corr = TangentialCorrelation(
        rmin=min_sep,
        rmax=max_sep,
        nbins=nbins,
        sep_units=sep_units,
    )
    
    # Compute correlation
    result = corr.compute(
        ra_lens, dec_lens, ra_src, dec_src,
        obs1, obs2, weights
    )
    
    # Log diagnostics
    valid = result["npairs"] > 0
    if np.any(valid):
        logger.info(
            f"Correlation computed: {valid.sum()}/{nbins} bins with data"
        )
        logger.info(
            f"Tangential signal: min={result['xi'][valid].min():.3e}, "
            f"max={result['xi'][valid].max():.3e}"
        )
    else:
        logger.warning("No valid bins in correlation function!")
    
    return result


def compute_all_ggl_correlations(
    ra_lens: np.ndarray,
    dec_lens: np.ndarray,
    ra_src: np.ndarray,
    dec_src: np.ndarray,
    gamma1: np.ndarray,
    gamma2: np.ndarray,
    F1: np.ndarray,
    F2: np.ndarray,
    G1: np.ndarray,
    G2: np.ndarray,
    kappa: Optional[np.ndarray] = None,
    weights: Optional[np.ndarray] = None,
    **kwargs,
) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Compute GGL correlations for all lensing observables.
    
    Computes tangential correlations for:
    - Shear (γ)
    - First flexion (F)
    - Second flexion (G)
    - Convergence (κ), if provided
    
    Args:
        ra_lens: Right ascension of lens galaxies (degrees)
        dec_lens: Declination of lens galaxies (degrees)
        ra_src: Right ascension of source galaxies (degrees)
        dec_src: Declination of source galaxies (degrees)
        gamma1, gamma2: Shear components
        F1, F2: First flexion components
        G1, G2: Second flexion components
        kappa: Optional convergence field
        weights: Optional source weights
        **kwargs: Additional arguments passed to compute_ggl_correlation
    
    Returns:
        Dictionary with keys 'shear', 'F', 'G', and optionally 'kappa',
        each containing a correlation result dictionary
    
    Example:
        >>> results = compute_all_ggl_correlations(
        ...     ra_lens, dec_lens, ra_src, dec_src,
        ...     gamma1, gamma2, F1, F2, G1, G2,
        ...     kappa=kappa, weights=weights
        ... )
        >>> shear_result = results['shear']
        >>> F_result = results['F']
    """
    logger.info("Computing all GGL correlations (γ, F, G, κ)")
    
    results = {}
    
    # Shear
    logger.info("Computing shear correlation...")
    results["shear"] = compute_ggl_correlation(
        ra_lens, dec_lens, ra_src, dec_src,
        gamma1, gamma2, weights, **kwargs
    )
    
    # First flexion
    logger.info("Computing first flexion correlation...")
    results["F"] = compute_ggl_correlation(
        ra_lens, dec_lens, ra_src, dec_src,
        F1, F2, weights, **kwargs
    )
    
    # Second flexion
    logger.info("Computing second flexion correlation...")
    results["G"] = compute_ggl_correlation(
        ra_lens, dec_lens, ra_src, dec_src,
        G1, G2, weights, **kwargs
    )
    
    # Convergence (if provided)
    if kappa is not None:
        logger.info("Computing convergence correlation...")
        # For scalar field like κ, use zeros for second component
        results["kappa"] = compute_ggl_correlation(
            ra_lens, dec_lens, ra_src, dec_src,
            kappa, np.zeros_like(kappa), weights, **kwargs
        )
    
    logger.info(f"Completed {len(results)} correlation computations")
    return results


def fit_ggl_profile(
    r: np.ndarray,
    xi: np.ndarray,
    xi_err: np.ndarray,
    profile_type: str = "nfw",
    observable: str = "shear",
    p0: Optional[Union[float, Tuple[float, float]]] = None,
) -> Dict[str, Union[float, np.ndarray]]:
    """
    Fit halo profile to GGL correlation function.
    
    Args:
        r: Separation (same units as used in profile, typically arcmin or kpc)
        xi: Tangential correlation amplitude
        xi_err: Error on tangential correlation
        profile_type: Profile to fit ('nfw' or 'sis')
        observable: Observable type ('shear', 'convergence', 'F', 'G')
        p0: Initial guess for parameters (ks, rs for NFW; theta_E for SIS)
    
    Returns:
        Dictionary with fit results:
        - For NFW: ks, rs, ks_err, rs_err, chi2, redchi2, dof, profile
        - For SIS: theta_E, theta_E_err, chi2, redchi2, dof, profile
    
    Raises:
        ValueError: If profile_type not recognized
    
    Note:
        Returns lensing parameters (κ_s, r_s) not physical (M200, c).
        Conversion to physical parameters requires cosmology integration.
    
    Example:
        >>> # Fit NFW to shear profile
        >>> fit = fit_ggl_profile(
        ...     r_arcmin, xi, xi_err,
        ...     profile_type='nfw', observable='shear'
        ... )
        >>> ks = fit['ks']
        >>> rs = fit['rs']
    """
    logger.info(f"Fitting {profile_type.upper()} profile to {observable}")
    
    if profile_type.lower() == "nfw":
        profile, fit_info = fit_nfw_profile(
            r, xi, xi_err,
            initial_guess=p0,
            observable=observable,
        )
        # Extract results
        result = {
            "ks": fit_info["params"][0],
            "rs": fit_info["params"][1],
            "ks_err": fit_info["errors"][0],
            "rs_err": fit_info["errors"][1],
            "chi2": fit_info["chi2"],
            "redchi2": fit_info["reduced_chi2"],
            "dof": fit_info["dof"],
            "profile": profile,
        }
        logger.info(
            f"NFW fit: κ_s={result['ks']:.4f}, r_s={result['rs']:.2f}, "
            f"chi2/dof={result['redchi2']:.2f}"
        )
    elif profile_type.lower() == "sis":
        profile, fit_info = fit_sis_profile(r, xi, xi_err)
        result = {
            "theta_E": fit_info["params"][0],
            "theta_E_err": fit_info["errors"][0],
            "chi2": fit_info["chi2"],
            "redchi2": fit_info["reduced_chi2"],
            "dof": fit_info["dof"],
            "profile": profile,
        }
        logger.info(
            f"SIS fit: θ_E={result['theta_E']:.3f}, "
            f"chi2/dof={result['redchi2']:.2f}"
        )
    else:
        raise ValueError(
            f"Unknown profile_type: {profile_type}. Use 'nfw' or 'sis'"
        )
    
    return result


def fit_all_ggl_profiles(
    correlations: Dict[str, Dict[str, np.ndarray]],
    profile_type: str = "nfw",
) -> Dict[str, Dict[str, Union[float, np.ndarray]]]:
    """
    Fit halo profiles to all GGL correlations.
    
    Args:
        correlations: Output from compute_all_ggl_correlations()
        profile_type: Profile to fit ('nfw' or 'sis')
    
    Returns:
        Dictionary with keys matching correlations, each containing fit results
    
    Example:
        >>> correlations = compute_all_ggl_correlations(...)
        >>> fits = fit_all_ggl_profiles(
        ...     correlations, profile_type='nfw'
        ... )
        >>> shear_fit = fits['shear']
        >>> F_fit = fits['F']
    """
    logger.info(f"Fitting {profile_type.upper()} profiles to all correlations")
    
    fits = {}
    
    # Map observable names to their types
    observable_map = {
        "shear": "shear",
        "F": "F",
        "G": "G",
        "kappa": "convergence",
    }
    
    for name, corr_result in correlations.items():
        logger.info(f"Fitting {name}...")
        
        # Get separation from correlation
        r = corr_result["r"]
        
        # Filter bins with data
        valid = (corr_result["npairs"] > 0) & np.isfinite(corr_result["xi"])
        
        if valid.sum() < 3:
            logger.warning(
                f"Skipping {name}: only {valid.sum()} valid bins "
                f"(need at least 3)"
            )
            continue
        
        # Fit profile
        try:
            fit_result = fit_ggl_profile(
                r[valid],
                corr_result["xi"][valid],
                corr_result["xi_err"][valid],
                profile_type=profile_type,
                observable=observable_map[name],
            )
            fits[name] = fit_result
        except Exception as e:
            logger.error(f"Failed to fit {name}: {e}")
    
    logger.info(f"Completed {len(fits)}/{len(correlations)} profile fits")
    return fits
