"""
HEALPix utilities for weak lensing analysis.

This module provides tools for working with HEALPix-formatted sky maps,
including projections to Cartesian patches for derivative calculations
and memory-efficient processing of large maps.
"""

import numpy as np
import healpy as hp
from typing import Tuple, Optional, Callable, Iterator
from astropy.io import fits
import logging

logger = logging.getLogger(__name__)


class HEALPixMap:
    """
    Wrapper for HEALPix sky maps with lensing-specific utilities.

    Parameters
    ----------
    nside : int
        HEALPix resolution parameter (npix = 12 * nside**2)
    data : np.ndarray
        Map data in HEALPix pixelization
    ordering : {'RING', 'NESTED'}, optional
        Pixel ordering scheme (default: 'RING')

    Attributes
    ----------
    nside : int
        Resolution parameter
    npix : int
        Total number of pixels
    data : np.ndarray
        Map data
    ordering : str
        Pixel ordering scheme
    """

    def __init__(self, nside: int, data: np.ndarray, ordering: str = "RING"):
        self.nside = nside
        self.npix = hp.nside2npix(nside)
        self.ordering = ordering.upper()

        if data.shape[0] != self.npix:
            raise ValueError(
                f"Data size {data.shape[0]} does not match " f"npix={self.npix} for nside={nside}"
            )

        if self.ordering not in ["RING", "NESTED"]:
            raise ValueError(f"Invalid ordering: {self.ordering}")

        self.data = np.asarray(data, dtype=np.float64)

    @classmethod
    def from_fits(
        cls, filename: str, field: int = 0, ordering: Optional[str] = None
    ) -> "HEALPixMap":
        """
        Read HEALPix map from FITS file.

        Parameters
        ----------
        filename : str
            Path to FITS file
        field : int, optional
            Field number if multiple maps in file (default: 0)
        ordering : str, optional
            Override ordering from FITS header

        Returns
        -------
        HEALPixMap
            Loaded map
        """
        # First read header to determine ordering
        _, header = hp.read_map(filename, field=field, h=True, verbose=False, nest=False)
        header_dict = dict(header)

        # Determine ordering from header if not specified
        if ordering is None:
            ordering = header_dict.get("ORDERING", "RING")

        ordering = ordering.upper()

        # Read with correct nest parameter
        nest_flag = ordering == "NESTED"
        data = hp.read_map(filename, field=field, verbose=False, nest=nest_flag)

        nside = hp.npix2nside(len(data))

        logger.info(f"Loaded HEALPix map from {filename}: " f"nside={nside}, ordering={ordering}")

        return cls(nside, data, ordering=ordering)

    def to_fits(
        self, filename: str, coord: str = "C", column_name: str = "SIGNAL", overwrite: bool = True
    ):
        """
        Write HEALPix map to FITS file.

        Parameters
        ----------
        filename : str
            Output FITS file path
        coord : str, optional
            Coordinate system ('C'=Celestial, 'G'=Galactic, 'E'=Ecliptic)
        column_name : str, optional
            Name for the map column (default: 'SIGNAL')
        overwrite : bool, optional
            Overwrite existing file (default: True)

        Note
        ----
        healpy.write_map behavior:
        - nest=False: writes data as-is in RING ordering
        - nest=True: converts data from RING to NESTED before writing

        To write NESTED data: we must provide it in RING form with nest=True,
        OR write it as-is with nest=False but this writes RING format.

        Simpler: Always write in the format our data is in, without conversion.
        """
        # Write directly without reordering
        # nest parameter tells healpy what ordering the data is in
        nest_flag = self.ordering == "NESTED"

        hp.write_map(
            filename,
            self.data,
            nest=nest_flag,
            coord=coord,
            column_names=[column_name],
            overwrite=overwrite,
            dtype=np.float64,
        )
        logger.info(f"Wrote HEALPix map to {filename} (ordering={self.ordering})")

    def get_neighbors(self, ipix: int) -> np.ndarray:
        """
        Get neighboring pixels for finite difference calculations.

        Parameters
        ----------
        ipix : int
            Pixel index

        Returns
        -------
        np.ndarray
            Array of up to 8 neighbor pixel indices
        """
        nest = self.ordering == "NESTED"
        neighbors = hp.get_all_neighbours(self.nside, ipix, nest=nest)
        # Filter out -1 (no neighbor exists, e.g., at poles)
        return neighbors[neighbors >= 0]

    def to_cartesian(
        self, center_ra: float, center_dec: float, radius_deg: float, npix: int
    ) -> Tuple[np.ndarray, dict]:
        """
        Project HEALPix map patch to Cartesian grid.

        Uses gnomonic (tangent plane) projection centered at (RA, Dec).
        Suitable for small patches where flat-sky approximation is valid.

        Parameters
        ----------
        center_ra : float
            Right ascension of patch center (degrees)
        center_dec : float
            Declination of patch center (degrees)
        radius_deg : float
            Radius of patch (degrees)
        npix : int
            Number of pixels along each axis of output grid

        Returns
        -------
        cartesian_map : np.ndarray
            Projected map on Cartesian grid, shape (npix, npix)
        info : dict
            Projection metadata (center, pixel_scale, etc.)
        """
        # Convert center to colatitude/longitude for healpy
        theta_c = np.radians(90.0 - center_dec)  # colatitude
        phi_c = np.radians(center_ra)

        # Create Cartesian grid in gnomonic projection
        # Range: [-radius, +radius] in degrees
        x = np.linspace(-radius_deg, radius_deg, npix)
        y = np.linspace(-radius_deg, radius_deg, npix)
        xx, yy = np.meshgrid(x, y, indexing="ij")

        # Convert Cartesian offsets to angular coordinates
        # Gnomonic projection: x = tan(θ)cos(φ), y = tan(θ)sin(φ)
        # where θ, φ are offsets from patch center
        r = np.sqrt(xx**2 + yy**2)

        # For small angles: approximate with flat-sky
        # This is valid for radius_deg << 10 degrees
        delta_theta = -yy * np.pi / 180  # offset in colatitude (radians)
        delta_phi = xx * np.pi / 180 / np.sin(theta_c)  # offset in longitude

        theta = theta_c + delta_theta
        phi = phi_c + delta_phi

        # Clip to valid range
        theta = np.clip(theta, 0, np.pi)
        phi = np.fmod(phi, 2 * np.pi)
        phi[phi < 0] += 2 * np.pi

        # Convert to HEALPix pixel indices
        nest = self.ordering == "NESTED"
        ipix = hp.ang2pix(self.nside, theta, phi, nest=nest)

        # Extract values
        cartesian_map = self.data[ipix]

        # Metadata
        pixel_scale_deg = 2 * radius_deg / npix
        info = {
            "center_ra": center_ra,
            "center_dec": center_dec,
            "radius_deg": radius_deg,
            "npix": npix,
            "pixel_scale_deg": pixel_scale_deg,
            "projection": "gnomonic",
            "nside": self.nside,
        }

        logger.info(
            f"Projected {npix}×{npix} patch centered at "
            f"(RA={center_ra:.2f}°, Dec={center_dec:.2f}°), "
            f"radius={radius_deg:.2f}°"
        )

        return cartesian_map, info

    def process_patches(
        self,
        patch_size_deg: float,
        overlap_deg: float,
        func: Callable[[np.ndarray, dict], np.ndarray],
    ) -> Iterator[Tuple[np.ndarray, dict]]:
        """
        Process map in overlapping patches (memory-efficient iterator).

        Divides the full sky into patches with specified size and overlap,
        applies a function to each patch, and yields results.

        Parameters
        ----------
        patch_size_deg : float
            Size of each patch (degrees)
        overlap_deg : float
            Overlap between adjacent patches (degrees)
        func : callable
            Function to apply to each patch: func(data, info) -> result

        Yields
        ------
        result : np.ndarray
            Result of func applied to patch
        info : dict
            Patch metadata (center, bounds, etc.)
        """
        # Simple all-sky tiling strategy
        # For production: use hierarchical pixelization

        # Number of patches along declination
        n_dec = int(180 / (patch_size_deg - overlap_deg)) + 1

        for i_dec in range(n_dec):
            # Declination range for this strip
            dec_min = -90 + i_dec * (patch_size_deg - overlap_deg)
            dec_max = dec_min + patch_size_deg
            dec_center = (dec_min + dec_max) / 2

            # Number of patches along RA (fewer near poles)
            cos_dec = np.cos(np.radians(dec_center))
            if cos_dec > 0.1:  # Not too close to pole
                n_ra = int(360 * cos_dec / (patch_size_deg - overlap_deg)) + 1
            else:
                n_ra = 1  # Single patch at pole

            for i_ra in range(n_ra):
                ra_center = i_ra * 360 / n_ra

                # Extract and process patch
                try:
                    patch_data, patch_info = self.to_cartesian(
                        ra_center,
                        dec_center,
                        radius_deg=patch_size_deg / 2,
                        npix=128,  # Fixed resolution for efficiency
                    )

                    result = func(patch_data, patch_info)
                    yield result, patch_info

                except Exception as e:
                    logger.warning(
                        f"Failed to process patch at "
                        f"(RA={ra_center:.1f}°, Dec={dec_center:.1f}°): {e}"
                    )
                    continue


def cartesian_to_healpix(
    ra: np.ndarray, dec: np.ndarray, nside: int, ordering: str = "RING"
) -> np.ndarray:
    """
    Convert RA/Dec coordinates to HEALPix pixel indices.

    Parameters
    ----------
    ra : np.ndarray
        Right ascension (degrees)
    dec : np.ndarray
        Declination (degrees)
    nside : int
        HEALPix resolution parameter
    ordering : {'RING', 'NESTED'}, optional
        Pixel ordering scheme (default: 'RING')

    Returns
    -------
    np.ndarray
        HEALPix pixel indices
    """
    theta = np.radians(90.0 - dec)  # colatitude
    phi = np.radians(ra)
    nest = ordering.upper() == "NESTED"

    return hp.ang2pix(nside, theta, phi, nest=nest)


def healpix_to_cartesian(
    ipix: np.ndarray, nside: int, ordering: str = "RING"
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert HEALPix pixel indices to RA/Dec coordinates.

    Parameters
    ----------
    ipix : np.ndarray
        HEALPix pixel indices
    nside : int
        HEALPix resolution parameter
    ordering : {'RING', 'NESTED'}, optional
        Pixel ordering scheme (default: 'RING')

    Returns
    -------
    ra : np.ndarray
        Right ascension (degrees)
    dec : np.ndarray
        Declination (degrees)
    """
    nest = ordering.upper() == "NESTED"
    theta, phi = hp.pix2ang(nside, ipix, nest=nest)

    dec = 90.0 - np.degrees(theta)  # Convert colatitude to declination
    ra = np.degrees(phi)

    return ra, dec
