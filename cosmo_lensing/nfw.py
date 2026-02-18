"""
NFW profile modeling for weak lensing (DEPRECATED).

This module is deprecated. Use cosmo_lensing.profiles instead.

For backward compatibility, NFWProfile and fit_nfw_profile are aliased here.
"""

import warnings
from cosmo_lensing.profiles import NFWProfile, fit_nfw_profile

warnings.warn(
    "cosmo_lensing.nfw is deprecated. Use cosmo_lensing.profiles instead.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["NFWProfile", "fit_nfw_profile"]

