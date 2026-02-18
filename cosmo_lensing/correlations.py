"""
Tangential correlation functions for galaxy-galaxy lensing.

Wraps TreeCorr for computing azimuthal averages of shear/convergence
around lens positions.

References:
    Jarvis et al. (2016) - TreeCorr: correlation functions with tree codes
    https://github.com/rmjarvis/TreeCorr
"""

import logging
from typing import Dict, Optional, Tuple

import numpy as np
import treecorr

logger = logging.getLogger(__name__)


class TangentialCorrelation:
    """
    Compute tangential shear correlation using TreeCorr.

    This class wraps TreeCorr's NGCorrelation for standardized
    tangential shear measurements with proper error estimation.

    The tangential shear profile γ_t(r) is computed using the correlation:
        γ_t = -Re[γ * exp(-2iφ)]
    where φ is the position angle from lens to source.

    TreeCorr handles the azimuthal averaging and coordinate transformations
    automatically.
    """

    def __init__(
        self,
        rmin: float,
        rmax: float,
        nbins: int,
        sep_units: str = "arcmin",
        bin_type: str = "Log",
        var_method: str = "jackknife",
        npatch: int = 10,
        metric: str = "Euclidean",
    ):
        """
        Initialize correlation calculator.

        Args:
            rmin: Minimum separation
            rmax: Maximum separation
            nbins: Number of radial bins
            sep_units: Units for separations ('arcmin', 'degrees', etc.)
            bin_type: Binning type ('Log' or 'Linear')
            var_method: Variance estimation method ('jackknife' or 'bootstrap')
            npatch: Number of patches for jackknife/bootstrap
            metric: Distance metric ('Euclidean', 'Arc', 'Periodic')
        """
        self.rmin = rmin
        self.rmax = rmax
        self.nbins = nbins
        self.sep_units = sep_units
        self.bin_type = bin_type
        self.var_method = var_method
        self.npatch = npatch
        self.metric = metric

        logger.info(
            f"Initialized TangentialCorrelation: "
            f"r=[{rmin}, {rmax}] {sep_units}, nbins={nbins}, "
            f"var_method={var_method}"
        )

    def compute(
        self,
        ra_lens: np.ndarray,
        dec_lens: np.ndarray,
        ra_src: np.ndarray,
        dec_src: np.ndarray,
        g1: np.ndarray,
        g2: np.ndarray,
        weights: Optional[np.ndarray] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Compute tangential shear correlation.

        Args:
            ra_lens: Right ascension of lenses (degrees)
            dec_lens: Declination of lenses (degrees)
            ra_src: Right ascension of sources (degrees)
            dec_src: Declination of sources (degrees)
            g1: Shear component 1
            g2: Shear component 2
            weights: Optional weights for sources

        Returns:
            Dictionary with keys:
            - 'r': Radial bin centers (in sep_units)
            - 'r_nom': Nominal bin centers (log or linear)
            - 'xi': Tangential shear correlation γ_t
            - 'xi_err': Error on correlation
            - 'xim': Cross component (should be ~0)
            - 'xim_err': Error on cross component
            - 'npairs': Number of pairs per bin
            - 'weight': Sum of weights per bin
        """
        logger.info(f"Computing correlation for {len(ra_lens)} lenses, " f"{len(ra_src)} sources")

        # Adjust npatch for small samples (must be same for both catalogs)
        npatch = min(self.npatch, len(ra_lens), len(ra_src))

        # Adjust var_method if not enough patches
        var_method = self.var_method if npatch > 1 else "shot"

        # Create TreeCorr catalogs
        lens_cat = treecorr.Catalog(
            ra=ra_lens, dec=dec_lens, ra_units="deg", dec_units="deg", npatch=npatch
        )

        if weights is None:
            weights = np.ones_like(g1)

        src_cat = treecorr.Catalog(
            ra=ra_src,
            dec=dec_src,
            g1=g1,
            g2=g2,
            w=weights,
            ra_units="deg",
            dec_units="deg",
            npatch=npatch,
        )

        # Configure correlation
        ng = treecorr.NGCorrelation(
            min_sep=self.rmin,
            max_sep=self.rmax,
            nbins=self.nbins,
            sep_units=self.sep_units,
            bin_type=self.bin_type,
            metric=self.metric,
            var_method=var_method,
        )

        # Compute correlation
        ng.process(lens_cat, src_cat)

        # Extract results
        result = {
            "r": np.exp(ng.meanlogr),  # Mean separation per bin
            "r_nom": ng.rnom,  # Nominal bin centers
            "xi": ng.xi,  # Tangential shear <γ_t>
            "xi_err": np.sqrt(ng.varxi),  # Error on γ_t
            "xim": ng.xi_im,  # Cross component (should be ~0)
            "xim_err": np.sqrt(ng.varxi),  # Symmetric error
            "npairs": ng.npairs,  # Number of pairs
            "weight": ng.weight,  # Sum of weights
        }

        logger.info(
            f"Correlation computed: {len(result['r'])} bins, "
            f"{np.sum(result['npairs'])} total pairs"
        )

        return result

    def compute_flexion(
        self,
        ra_lens: np.ndarray,
        dec_lens: np.ndarray,
        ra_src: np.ndarray,
        dec_src: np.ndarray,
        F1: np.ndarray,
        F2: np.ndarray,
        G1: np.ndarray,
        G2: np.ndarray,
        weights: Optional[np.ndarray] = None,
    ) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
        """
        Compute first and second flexion correlations.

        Flexion is a third-order lensing observable:
        - F (first flexion): gradient of convergence
        - G (second flexion): gradient of shear

        Args:
            ra_lens: Right ascension of lenses (degrees)
            dec_lens: Declination of lenses (degrees)
            ra_src: Right ascension of sources (degrees)
            dec_src: Declination of sources (degrees)
            F1, F2: First flexion components
            G1, G2: Second flexion components
            weights: Optional weights for sources

        Returns:
            (F_result, G_result): Tuple of dictionaries with correlation results
        """
        logger.info("Computing flexion correlations")

        # F correlation (treat like shear)
        F_result = self.compute(ra_lens, dec_lens, ra_src, dec_src, F1, F2, weights)

        # G correlation (treat like shear)
        G_result = self.compute(ra_lens, dec_lens, ra_src, dec_src, G1, G2, weights)

        return F_result, G_result
