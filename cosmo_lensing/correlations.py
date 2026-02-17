"""
Tangential correlation functions for galaxy-galaxy lensing.

Wraps TreeCorr for computing azimuthal averages of shear/convergence
around lens positions.
"""

import logging
from typing import Dict, Optional

import numpy as np

logger = logging.getLogger(__name__)


class TangentialCorrelation:
    """
    Compute tangential shear correlation using TreeCorr.
    
    This class will wrap TreeCorr's NGCorrelation for standardized
    tangential shear measurements.
    
    TODO: Full implementation in Phase 3
    """
    
    def __init__(
        self,
        rmin: float,
        rmax: float,
        nbins: int,
        sep_units: str = 'arcmin',
        var_method: str = 'jackknife'
    ):
        """
        Initialize correlation calculator.
        
        Args:
            rmin: Minimum separation
            rmax: Maximum separation
            nbins: Number of radial bins
            sep_units: Units for separations ('arcmin', 'degrees', etc.)
            var_method: Variance estimation method ('jackknife' or 'bootstrap')
        """
        self.rmin = rmin
        self.rmax = rmax
        self.nbins = nbins
        self.sep_units = sep_units
        self.var_method = var_method
        
        logger.info(
            f"Initialized TangentialCorrelation: "
            f"r=[{rmin}, {rmax}] {sep_units}, nbins={nbins}"
        )
    
    def compute(
        self,
        ra_lens: np.ndarray,
        dec_lens: np.ndarray,
        ra_src: np.ndarray,
        dec_src: np.ndarray,
        g1: np.ndarray,
        g2: np.ndarray,
        weights: Optional[np.ndarray] = None
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
            - 'r': Radial bin centers
            - 'xi': Tangential correlation
            - 'xi_err': Error on correlation
            - 'npairs': Number of pairs per bin
        """
        # TODO: Implement full TreeCorr wrapper in Phase 3
        raise NotImplementedError("TangentialCorrelation.compute() - Phase 3")
