"""
Weak gravitational lensing analysis pipeline for cosmological N-body simulations.

This package provides tools for:
- Reading deflection fields from Fortran binary files
- Computing lensing observables (convergence, shear, flexion)
- Measuring tangential correlations (galaxy-galaxy lensing)
- Fitting NFW halo profiles
- Generating publication-quality plots

Modules:
    io: Input/output operations for deflection fields and catalogs
    cosmology: Cosmological calculations (distances, critical density)
    derivatives: Finite-difference derivatives for Jacobian computation
    observables: Lensing observables (κ, γ, F, G)
    correlations: Tangential correlation functions
    nfw: NFW profile modeling
"""

__version__ = "0.1.0"
__author__ = "Laurent Magri-Stella"

from . import io, cosmology, derivatives, observables, correlations, nfw

__all__ = [
    "io",
    "cosmology",
    "derivatives",
    "observables",
    "correlations",
    "nfw",
]
