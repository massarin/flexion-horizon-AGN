"""
Weak gravitational lensing analysis pipeline for cosmological N-body simulations.

This package provides tools for:
- Reading deflection fields from Fortran binary files
- Computing lensing observables (convergence, shear, flexion)
- Measuring tangential correlations (galaxy-galaxy lensing)
- Fitting NFW and SIS halo profiles
- Generating publication-quality plots

Modules:
    io: Input/output operations for deflection fields and catalogs
    cosmology: Cosmological calculations (distances, critical density)
    derivatives: Finite-difference derivatives for Jacobian computation
    observables: Lensing observables (κ, γ, F, G)
    correlations: Tangential correlation functions
    profiles: NFW and SIS profile models (replaces nfw)
    nfw: Deprecated, use profiles instead
    healpix_utils: HEALPix sky map utilities
    synthetic_data: Test field generation
"""

__version__ = "1.0.0"
__author__ = "Laurent Magri-Stella"

from . import (
    io,
    cosmology,
    derivatives,
    observables,
    correlations,
    profiles,
    nfw,  # Deprecated
    healpix_utils,
    synthetic_data,
)

__all__ = [
    "io",
    "cosmology",
    "derivatives",
    "observables",
    "correlations",
    "profiles",
    "nfw",
    "healpix_utils",
    "synthetic_data",
]
