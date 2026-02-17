"""
NFW profile modeling for weak lensing.

Implements NFW (Navarro-Frenk-White) halo profile calculations
for convergence and shear.
"""

import logging
from typing import Tuple

import numpy as np

logger = logging.getLogger(__name__)


class NFWProfile:
    """
    NFW halo profile for lensing calculations.
    
    The NFW profile is characterized by:
    - κ_s: convergence scale
    - r_s: scale radius
    
    References:
        Wright & Brainerd (2000): "Gravitational lensing by NFW halos"
            ApJ 534, 34-40
    """
    
    def __init__(self, ks: float, rs: float):
        """
        Initialize NFW profile.
        
        Args:
            ks: Convergence scale (dimensionless)
            rs: Scale radius (same units as r in lensing calculations)
        """
        self.ks = ks
        self.rs = rs
        
        logger.info(f"Initialized NFW profile: κ_s={ks:.4f}, r_s={rs:.2f}")
    
    @staticmethod
    def _F_function(x: np.ndarray) -> np.ndarray:
        """
        Auxiliary function F(x) for NFW lensing.
        
        F(x) = { (1 - 2*arctanh(sqrt((1-x)/(1+x)))/sqrt(1-x²))/(x²-1)  for x < 1
               { 1/3                                                    for x = 1
               { (1 - 2*arctan(sqrt((x-1)/(1+x)))/sqrt(x²-1))/(x²-1)   for x > 1
        
        Args:
            x: Dimensionless radius r/r_s
            
        Returns:
            F(x) values
            
        References:
            Wright & Brainerd (2000) Eq. 13
        """
        # TODO: Implement F function in Phase 3
        raise NotImplementedError("NFWProfile._F_function() - Phase 3")
    
    def convergence(self, r: np.ndarray) -> np.ndarray:
        """
        Compute NFW convergence profile κ(r).
        
        Args:
            r: Radial distances (same units as r_s)
            
        Returns:
            Convergence values
            
        References:
            Wright & Brainerd (2000) Eq. 11
        """
        # TODO: Implement convergence profile in Phase 3
        raise NotImplementedError("NFWProfile.convergence() - Phase 3")
    
    def shear_tangential(self, r: np.ndarray) -> np.ndarray:
        """
        Compute NFW tangential shear profile γ_t(r).
        
        Args:
            r: Radial distances (same units as r_s)
            
        Returns:
            Tangential shear values
            
        References:
            Wright & Brainerd (2000) Eq. 12
        """
        # TODO: Implement tangential shear in Phase 3
        raise NotImplementedError("NFWProfile.shear_tangential() - Phase 3")
