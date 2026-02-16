"""
Hale Telescope optical model using HCIpy.

The 200-inch Hale Telescope at Palomar Observatory:
- 5.08m primary mirror
- 1.83m central obscuration (secondary)
- Four-vane spider support structure
"""

import numpy as np
import hcipy as hp
from typing import Tuple, Optional

from config import (
    HALE_DIAMETER,
    HALE_OBSCURATION_RATIO,
    WAVELENGTH
)


# Spider vane parameters for Hale telescope
SPIDER_WIDTH = 0.02  # Relative width (fraction of diameter)
SPIDER_ANGLES = [0, 90, 180, 270]  # Four-vane configuration (degrees)


def create_hale_pupil(
    pupil_grid: hp.Grid,
    include_spider: bool = True,
    spider_width: float = SPIDER_WIDTH
) -> hp.Field:
    """
    Create the Hale telescope pupil function.
    
    Parameters
    ----------
    pupil_grid : hp.Grid
        The pupil plane grid (should be large enough for pupil)
    include_spider : bool
        Whether to include spider vane obscuration
    spider_width : float
        Spider vane width as fraction of diameter
    
    Returns
    -------
    hp.Field
        Complex pupil function (1 inside, 0 outside)
    """
    # Primary aperture (circular)
    primary = hp.circular_aperture(HALE_DIAMETER)(pupil_grid)
    
    # Secondary obscuration (central)
    secondary_diameter = HALE_DIAMETER * HALE_OBSCURATION_RATIO
    secondary = hp.circular_aperture(secondary_diameter)(pupil_grid)
    
    # Pupil = primary - secondary
    pupil = primary - secondary
    
    # Add spider vanes if requested
    if include_spider:
        spider = create_spider_vanes(pupil_grid, spider_width)
        pupil = pupil * (1 - spider)
    
    # Ensure binary (0 or 1)
    pupil = hp.Field(np.clip(pupil, 0, 1), pupil_grid)
    
    return pupil


def create_spider_vanes(
    pupil_grid: hp.Grid,
    width: float = SPIDER_WIDTH
) -> hp.Field:
    """
    Create spider vane obscuration pattern.
    
    Parameters
    ----------
    pupil_grid : hp.Grid
        The pupil plane grid
    width : float
        Vane width as fraction of diameter
    
    Returns
    -------
    hp.Field
        Spider vane pattern (1 where obscured)
    """
    spider = hp.Field(np.zeros(pupil_grid.size), pupil_grid)
    
    # Get coordinates
    x, y = pupil_grid.coords
    
    # Vane width in meters
    vane_width = width * HALE_DIAMETER
    
    # Create four vanes at 0, 90, 180, 270 degrees
    for angle in SPIDER_ANGLES:
        theta = np.radians(angle)
        
        # Rotate coordinates
        x_rot = x * np.cos(theta) + y * np.sin(theta)
        y_rot = -x * np.sin(theta) + y * np.cos(theta)
        
        # Vane is where |y_rot| < width/2 and x_rot > 0
        # (from center to edge)
        vane_mask = (np.abs(y_rot) < vane_width / 2) & (x_rot > 0)
        
        spider += vane_mask
    
    # Clip to binary
    spider = hp.Field(np.clip(spider, 0, 1), pupil_grid)
    
    return spider


def create_pupil_grid(
    num_pixels: int = 256,
    oversampling: float = 1.5
) -> hp.Grid:
    """
    Create an appropriate pupil plane grid.
    
    Parameters
    ----------
    num_pixels : int
        Number of pixels across the grid
    oversampling : float
        Factor by which grid extends beyond pupil diameter
    
    Returns
    -------
    hp.Grid
        Cartesian grid for the pupil plane
    """
    grid_size = HALE_DIAMETER * oversampling
    pupil_grid = hp.make_pupil_grid(num_pixels, grid_size)
    return pupil_grid


def create_focal_grid(
    pupil_grid: hp.Grid,
    wavelength: float = WAVELENGTH,
    q: int = 4,
    num_airy: int = 16
) -> hp.Grid:
    """
    Create focal plane grid matched to pupil grid.
    
    Parameters
    ----------
    pupil_grid : hp.Grid
        The pupil plane grid (in physical units, meters)
    wavelength : float
        Observing wavelength in meters
    q : int
        Oversampling factor (pixels per lambda/D)
    num_airy : int
        Number of Airy rings to include (field of view radius in lambda/D)
    
    Returns
    -------
    hp.Grid
        Focal plane grid in physical units (radians)
    """
    # Use physical units: pupil_diameter, reference_wavelength, focal_length
    # focal_length=1 means focal plane coordinates are in radians
    focal_grid = hp.make_focal_grid(
        q, 
        num_airy,
        pupil_diameter=HALE_DIAMETER,
        reference_wavelength=wavelength,
        focal_length=1.0
    )
    return focal_grid


def create_propagator(
    pupil_grid: hp.Grid,
    focal_grid: hp.Grid,
    wavelength: float = WAVELENGTH
) -> hp.FraunhoferPropagator:
    """
    Create Fraunhofer propagator from pupil to focal plane.
    
    Parameters
    ----------
    pupil_grid : hp.Grid
        Pupil plane grid (in meters)
    focal_grid : hp.Grid
        Focal plane grid (in radians, with focal_length=1)
    wavelength : float
        Wavelength in meters (not used directly, but kept for API)
    
    Returns
    -------
    hp.FraunhoferPropagator
        Propagator object
    """
    # focal_length=1.0 means focal plane is in radians
    propagator = hp.FraunhoferPropagator(pupil_grid, focal_grid, focal_length=1.0)
    return propagator


def create_wavefront(
    pupil_grid: hp.Grid,
    wavelength: float = WAVELENGTH
) -> hp.Wavefront:
    """
    Create a plane wavefront at the telescope pupil.
    
    Parameters
    ----------
    pupil_grid : hp.Grid
        Pupil plane grid
    wavelength : float
        Wavelength in meters
    
    Returns
    -------
    hp.Wavefront
        Plane wavefront with unit amplitude
    """
    # Create uniform electric field
    E_field = hp.Field(np.ones(pupil_grid.size, dtype=complex), pupil_grid)
    
    # Create wavefront
    wf = hp.Wavefront(E_field, wavelength)
    
    return wf


class HaleTelescope:
    """
    Complete Hale telescope optical model.
    
    Combines pupil, propagator, and detector simulation.
    """
    
    def __init__(
        self,
        num_pixels: int = 256,
        wavelength: float = WAVELENGTH,
        focal_sampling: int = 4,
        num_airy: int = 16,
        include_spider: bool = True
    ):
        """
        Initialize Hale telescope model.
        
        Parameters
        ----------
        num_pixels : int
            Pupil grid resolution
        wavelength : float
            Observing wavelength
        focal_sampling : int
            Focal plane oversampling (pixels per lambda/D)
        num_airy : int
            Field of view in Airy rings
        include_spider : bool
            Include spider vane obscuration
        """
        self.wavelength = wavelength
        self.num_pixels = num_pixels
        
        # Create grids
        self.pupil_grid = create_pupil_grid(num_pixels)
        self.focal_grid = create_focal_grid(
            self.pupil_grid, wavelength, focal_sampling, num_airy
        )
        
        # Create pupil
        self.pupil = create_hale_pupil(self.pupil_grid, include_spider)
        
        # Create propagator
        self.propagator = create_propagator(
            self.pupil_grid, self.focal_grid, wavelength
        )
        
        # Compute reference PSF (diffraction-limited)
        self._compute_reference_psf()
    
    def _compute_reference_psf(self):
        """Compute diffraction-limited PSF for normalization."""
        wf = create_wavefront(self.pupil_grid, self.wavelength)
        wf.electric_field *= self.pupil
        
        focal_wf = self.propagator(wf)
        self.reference_psf = focal_wf.intensity
        self.strehl_normalization = self.reference_psf.max()
    
    def propagate(self, wavefront: hp.Wavefront) -> hp.Wavefront:
        """
        Propagate wavefront through telescope.
        
        Parameters
        ----------
        wavefront : hp.Wavefront
            Input wavefront (at pupil)
        
        Returns
        -------
        hp.Wavefront
            Output wavefront (at focal plane)
        """
        # Apply pupil
        wf = wavefront.copy()
        wf.electric_field *= self.pupil
        
        # Propagate to focal plane
        focal_wf = self.propagator(wf)
        
        return focal_wf
    
    def get_psf(self, phase_screen: Optional[hp.Field] = None) -> hp.Field:
        """
        Compute PSF with optional atmospheric phase.
        
        Parameters
        ----------
        phase_screen : hp.Field, optional
            Phase screen to apply (radians)
        
        Returns
        -------
        hp.Field
            PSF intensity
        """
        wf = create_wavefront(self.pupil_grid, self.wavelength)
        
        # Apply phase screen if provided
        if phase_screen is not None:
            wf.electric_field *= np.exp(1j * phase_screen)
        
        # Propagate through telescope
        focal_wf = self.propagate(wf)
        
        return focal_wf.intensity
    
    def compute_strehl(self, psf: hp.Field) -> float:
        """
        Compute Strehl ratio of a PSF.
        
        Parameters
        ----------
        psf : hp.Field
            PSF intensity
        
        Returns
        -------
        float
            Strehl ratio (0 to 1)
        """
        return psf.max() / self.strehl_normalization
    
    @property
    def plate_scale(self) -> float:
        """
        Focal plane plate scale in arcsec/pixel.
        
        Returns
        -------
        float
            Plate scale
        """
        # With physical units, focal_grid.delta is in radians
        pixel_size_rad = self.focal_grid.delta[0]
        pixel_size_arcsec = np.degrees(pixel_size_rad) * 3600
        
        return pixel_size_arcsec


def compute_diffraction_limit(wavelength: float = WAVELENGTH) -> dict:
    """
    Compute diffraction-limited performance metrics.
    
    Parameters
    ----------
    wavelength : float
        Wavelength in meters
    
    Returns
    -------
    dict
        Dictionary with lambda/D, FWHM, first null, etc.
    """
    lambda_D_rad = wavelength / HALE_DIAMETER
    lambda_D_arcsec = np.degrees(lambda_D_rad) * 3600
    lambda_D_mas = lambda_D_arcsec * 1000
    
    return {
        'lambda_D_rad': lambda_D_rad,
        'lambda_D_arcsec': lambda_D_arcsec,
        'lambda_D_mas': lambda_D_mas,
        'airy_fwhm_arcsec': 1.028 * lambda_D_arcsec,
        'first_null_arcsec': 1.22 * lambda_D_arcsec,
        'second_null_arcsec': 2.23 * lambda_D_arcsec
    }


if __name__ == "__main__":
    print("Hale Telescope Optical Model")
    print("=" * 50)
    
    print(f"\nPrimary diameter: {HALE_DIAMETER} m")
    print(f"Central obscuration: {HALE_OBSCURATION_RATIO*100:.0f}%")
    print(f"Wavelength: {WAVELENGTH*1e9:.0f} nm")
    
    limits = compute_diffraction_limit()
    print(f"\nDiffraction limit:")
    print(f"  λ/D = {limits['lambda_D_arcsec']*1000:.1f} mas")
    print(f"  Airy FWHM = {limits['airy_fwhm_arcsec']*1000:.1f} mas")
    print(f"  First null = {limits['first_null_arcsec']*1000:.1f} mas")
    
    # Test pupil creation
    print("\nCreating telescope model...")
    tel = HaleTelescope(num_pixels=256)
    
    # Compute collecting area
    pupil_sum = tel.pupil.sum() * tel.pupil_grid.delta[0]**2
    print(f"Effective collecting area: {pupil_sum:.2f} m^2")
    print(f"Plate scale: {tel.plate_scale*1000:.2f} mas/pixel")
