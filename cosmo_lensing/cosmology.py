"""
Cosmological calculations for weak lensing.

Replaces external C library (libsoftlens.so) with astropy.cosmology.
All cosmological quantities computed using FlatLambdaCDM model.

References:
    - Hogg (1999): "Distance measures in cosmology" (arXiv:astro-ph/9905116)
    - Wright & Brainerd (2000): NFW lensing formalism
"""

import logging
from typing import Union, Tuple

import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

logger = logging.getLogger(__name__)


class CosmologyCalculator:
    """
    Wrapper for cosmological calculations using astropy.
    
    Default cosmology: Planck 2018 (TT,TE,EE+lowE+lensing)
    - H0 = 67.66 km/s/Mpc
    - Om0 = 0.3111
    - Flat ΛCDM (Ωk = 0)
    
    Attributes:
        cosmo: astropy FlatLambdaCDM instance
        h: Reduced Hubble constant (H0 / 100)
    """
    
    def __init__(
        self,
        H0: float = 67.66,
        Om0: float = 0.3111,
        Ob0: float = 0.0490,
        Tcmb0: float = 2.7255
    ):
        """
        Initialize cosmology calculator.
        
        Args:
            H0: Hubble constant in km/s/Mpc
            Om0: Matter density parameter at z=0
            Ob0: Baryon density parameter at z=0
            Tcmb0: CMB temperature at z=0 in Kelvin
        """
        self.cosmo = FlatLambdaCDM(
            H0=H0 * u.km / u.s / u.Mpc,
            Om0=Om0,
            Ob0=Ob0,
            Tcmb0=Tcmb0 * u.K
        )
        self.h = H0 / 100.0
        
        logger.info(
            f"Initialized cosmology: H0={H0:.2f}, Om0={Om0:.4f}, "
            f"Ob0={Ob0:.4f}, h={self.h:.4f}"
        )
    
    def angular_diameter_distance(
        self,
        z1: Union[float, np.ndarray],
        z2: Union[float, np.ndarray] = 0.0
    ) -> Union[float, np.ndarray]:
        """
        Compute angular diameter distance D_A(z1, z2).
        
        For z2 = 0, returns D_A(z1) = D_C(z1) / (1 + z1)
        For z2 > z1, returns D_A(z1, z2) = D_C(z1, z2) / (1 + z2)
        
        Replaces jsl_dda from libsoftlens.so
        
        Args:
            z1: Observer/lens redshift
            z2: Source redshift (default 0)
            
        Returns:
            Angular diameter distance in Mpc
            
        References:
            Hogg (1999) Eq. 18
        """
        # Handle z2 = 0 case (distance from observer)
        if np.any(z2 == 0):
            result = self.cosmo.angular_diameter_distance(z1).value
        else:
            # Distance between two redshifts
            result = self.cosmo.angular_diameter_distance_z1z2(z1, z2).value
        
        return result
    
    def comoving_distance(
        self,
        z: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """
        Compute comoving distance D_C(z).
        
        Args:
            z: Redshift
            
        Returns:
            Comoving distance in Mpc
        """
        return self.cosmo.comoving_distance(z).value
    
    def critical_density(
        self,
        z: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """
        Compute critical density ρ_crit(z).
        
        ρ_crit(z) = 3 H(z)² / (8π G)
        
        Args:
            z: Redshift
            
        Returns:
            Critical density in Msun/Mpc³
        """
        rho_crit = self.cosmo.critical_density(z)
        # Convert to Msun/Mpc³
        rho_crit_mpc3 = rho_crit.to(u.Msun / u.Mpc**3).value
        return rho_crit_mpc3
    
    def Omega_m(
        self,
        z: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """
        Compute matter density parameter Ωm(z).
        
        Ωm(z) = Ωm0 (1 + z)³ / E²(z)
        where E(z) = H(z) / H0
        
        Args:
            z: Redshift
            
        Returns:
            Matter density parameter (dimensionless)
        """
        return self.cosmo.Om(z)
    
    def sigma_crit(
        self,
        zl: Union[float, np.ndarray],
        zs: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """
        Compute critical surface density Σ_crit for lensing.
        
        Σ_crit = (c² / 4πG) * (D_s / (D_l * D_ls))
        
        where:
        - D_l = angular diameter distance to lens
        - D_s = angular diameter distance to source
        - D_ls = angular diameter distance from lens to source
        
        Args:
            zl: Lens redshift
            zs: Source redshift
            
        Returns:
            Critical surface density in Msun/Mpc²
            
        Raises:
            ValueError: If zs <= zl (source behind lens)
            
        References:
            Bartelmann & Schneider (2001) Eq. 3.34
        """
        # Validate redshifts
        if np.any(zs <= zl):
            raise ValueError(
                f"Source redshift ({zs}) must be greater than lens redshift ({zl})"
            )
        
        # Angular diameter distances in Mpc
        D_l = self.angular_diameter_distance(zl)
        D_s = self.angular_diameter_distance(zs)
        D_ls = self.angular_diameter_distance(zl, zs)
        
        # Constants
        c = 299792.458  # km/s (speed of light)
        G_Mpc_Msun = 4.30091e-9  # G in Mpc³ / (Msun * s²)
        
        # Σ_crit in Msun/Mpc²
        sigma_crit = (c**2 / (4 * np.pi * G_Mpc_Msun)) * (D_s / (D_l * D_ls))
        
        return sigma_crit
    
    def distance_scaling(
        self,
        z: Union[float, np.ndarray]
    ) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
        """
        Compute useful distance scaling factors for lensing.
        
        Returns:
        - arcsec2kpc: Conversion factor from arcsec to physical kpc at redshift z
        - sigma_crit: Critical surface density for z_lens=0, z_source=z
        
        Args:
            z: Redshift
            
        Returns:
            Tuple of (arcsec2kpc, sigma_crit_ref)
            - arcsec2kpc: Physical size of 1 arcsec in kpc at redshift z
            - sigma_crit_ref: Reference Σ_crit in Msun/Mpc²
        """
        # Angular diameter distance in Mpc
        D_A = self.angular_diameter_distance(z)
        
        # Conversion: 1 arcsec = D_A * (1 arcsec in radians) = D_A * (1/206265) Mpc
        # Convert to kpc: * 1000
        arcsec2kpc = D_A * (1.0 / 206265.0) * 1000.0  # kpc
        
        # Reference critical surface density (for z_lens ≈ 0)
        # Use small z_lens = 0.01 to avoid singularity
        try:
            sigma_crit_ref = self.sigma_crit(0.01, z)
        except ValueError:
            # If z < 0.01, just return large value
            sigma_crit_ref = 1e20
        
        return arcsec2kpc, sigma_crit_ref
