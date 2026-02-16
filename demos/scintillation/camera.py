"""
Speckle camera simulation using HCIpy.

Implements short-exposure imaging with:
- Photon noise (Poisson statistics)
- Detector read noise
- Atmospheric speckle evolution
"""

import numpy as np
import hcipy as hp
from typing import Optional, Tuple, List

from config import (
    WAVELENGTH,
    EXPOSURE_TIME,
    photons_per_exposure
)


class SpeckleCamera:
    """
    Speckle imaging camera simulator.
    
    Models short-exposure imaging through atmospheric turbulence
    with realistic noise characteristics.
    """
    
    def __init__(
        self,
        focal_grid: hp.Grid,
        exposure_time: float = EXPOSURE_TIME,
        read_noise: float = 3.0,
        dark_current: float = 0.01,
        quantum_efficiency: float = 0.9,
        well_depth: int = 100000,
        gain: float = 1.0
    ):
        """
        Initialize speckle camera.
        
        Parameters
        ----------
        focal_grid : hp.Grid
            Focal plane grid from telescope
        exposure_time : float
            Single exposure time in seconds
        read_noise : float
            Read noise in electrons RMS
        dark_current : float
            Dark current in electrons/pixel/second
        quantum_efficiency : float
            Detector QE (0-1)
        well_depth : int
            Full well depth in electrons
        gain : float
            Detector gain (electrons per ADU)
        """
        self.focal_grid = focal_grid
        self.exposure_time = exposure_time
        self.read_noise = read_noise
        self.dark_current = dark_current
        self.quantum_efficiency = quantum_efficiency
        self.well_depth = well_depth
        self.gain = gain
        
        # Number of pixels
        self.num_pixels = focal_grid.size
    
    def expose(
        self,
        psf: hp.Field,
        total_photons: float,
        add_noise: bool = True,
        seed: Optional[int] = None
    ) -> hp.Field:
        """
        Simulate a single exposure.
        
        Parameters
        ----------
        psf : hp.Field
            Normalized PSF (should sum to 1)
        total_photons : float
            Total photons incident on detector
        add_noise : bool
            Whether to add photon and read noise
        seed : int, optional
            Random seed for reproducibility
        
        Returns
        -------
        hp.Field
            Detector image in electrons
        """
        if seed is not None:
            np.random.seed(seed)
        
        # Normalize PSF to sum to 1
        psf_normalized = psf / psf.sum()
        
        # Expected electrons per pixel
        expected_electrons = psf_normalized * total_photons * self.quantum_efficiency
        
        if add_noise:
            # Add dark current
            expected_electrons += self.dark_current * self.exposure_time
            
            # Photon noise (Poisson)
            electrons = np.random.poisson(expected_electrons)
            
            # Read noise (Gaussian)
            electrons = electrons + np.random.normal(0, self.read_noise, self.num_pixels)
            
            # Clip to well depth
            electrons = np.clip(electrons, 0, self.well_depth)
        else:
            electrons = expected_electrons
        
        return hp.Field(electrons, self.focal_grid)
    
    def expose_sequence(
        self,
        psfs: List[hp.Field],
        total_photons_per_frame: float,
        add_noise: bool = True,
        seed: Optional[int] = None
    ) -> List[hp.Field]:
        """
        Simulate a sequence of exposures.
        
        Parameters
        ----------
        psfs : list of hp.Field
            List of PSFs for each frame
        total_photons_per_frame : float
            Photons per exposure
        add_noise : bool
            Whether to add noise
        seed : int, optional
            Random seed
        
        Returns
        -------
        list of hp.Field
            List of detector images
        """
        if seed is not None:
            np.random.seed(seed)
        
        images = []
        for psf in psfs:
            img = self.expose(psf, total_photons_per_frame, add_noise)
            images.append(img)
        
        return images


class SpeckleDataCube:
    """
    Container for speckle image data cubes.
    
    Stores and processes sequences of short-exposure images.
    """
    
    def __init__(
        self,
        images: List[hp.Field],
        focal_grid: hp.Grid,
        exposure_time: float = EXPOSURE_TIME
    ):
        """
        Initialize data cube.
        
        Parameters
        ----------
        images : list of hp.Field
            List of speckle images
        focal_grid : hp.Grid
            Focal plane grid
        exposure_time : float
            Single frame exposure time
        """
        self.focal_grid = focal_grid
        self.exposure_time = exposure_time
        self.num_frames = len(images)
        
        # Convert to 3D numpy array
        self.data = np.array([img.shaped for img in images])
    
    def get_frame(self, index: int) -> hp.Field:
        """Get single frame as Field."""
        return hp.Field(self.data[index].ravel(), self.focal_grid)
    
    def mean_image(self) -> hp.Field:
        """Compute mean (long-exposure) image."""
        mean = self.data.mean(axis=0)
        return hp.Field(mean.ravel(), self.focal_grid)
    
    def power_spectrum_average(self) -> np.ndarray:
        """
        Compute average power spectrum (for speckle imaging).
        
        Returns
        -------
        np.ndarray
            Average of |FT(image)|^2
        """
        ps_sum = np.zeros(self.data[0].shape)
        
        for frame in self.data:
            ft = np.fft.fft2(frame)
            ps_sum += np.abs(ft)**2
        
        return np.fft.fftshift(ps_sum / self.num_frames)
    
    def bispectrum_average(self) -> np.ndarray:
        """
        Compute average bispectrum (for phase recovery).
        
        Returns a representative slice of the bispectrum.
        
        Returns
        -------
        np.ndarray
            Average bispectrum
        """
        # Full bispectrum is 4D, compute representative slice
        ny, nx = self.data[0].shape
        
        # Initialize bispectrum accumulator
        bispec_sum = np.zeros((ny, nx), dtype=complex)
        
        for frame in self.data:
            ft = np.fft.fft2(frame)
            
            # Compute bispectrum slice: B(u, v) = F(u)F(v)F*(-u-v)
            # Using u = (ux, 0), v = (0, vy) slice
            for ux in range(nx//4, 3*nx//4):
                for vy in range(ny//4, 3*ny//4):
                    bispec_sum[vy, ux] += (
                        ft[0, ux] * 
                        ft[vy, 0] * 
                        np.conj(ft[vy, ux])
                    )
        
        return bispec_sum / self.num_frames
    
    def autocorrelation_average(self) -> np.ndarray:
        """
        Compute average autocorrelation (Knox-Thompson method).
        
        Returns
        -------
        np.ndarray
            Average autocorrelation
        """
        ac_sum = np.zeros(self.data[0].shape)
        
        for frame in self.data:
            ft = np.fft.fft2(frame)
            ac = np.fft.ifft2(np.abs(ft)**2)
            ac_sum += np.real(np.fft.fftshift(ac))
        
        return ac_sum / self.num_frames


def compute_snr(
    magnitude: float,
    exposure_time: float = EXPOSURE_TIME,
    read_noise: float = 3.0,
    sky_background: float = 0.0,
    npix: int = 100
) -> dict:
    """
    Compute signal-to-noise ratio for speckle imaging.
    
    Parameters
    ----------
    magnitude : float
        Stellar magnitude
    exposure_time : float
        Exposure time in seconds
    read_noise : float
        Read noise in electrons
    sky_background : float
        Sky background in electrons/pixel/second
    npix : int
        Number of pixels in aperture
    
    Returns
    -------
    dict
        SNR metrics
    """
    photons = photons_per_exposure(magnitude, exposure=exposure_time)
    
    # Signal (electrons)
    signal = photons * 0.9  # Assuming 90% QE
    
    # Noise sources (variance)
    photon_noise_var = signal  # Poisson
    read_noise_var = npix * read_noise**2
    sky_noise_var = npix * sky_background * exposure_time
    
    total_noise = np.sqrt(photon_noise_var + read_noise_var + sky_noise_var)
    
    snr = signal / total_noise
    
    return {
        'signal_electrons': signal,
        'photon_noise': np.sqrt(photon_noise_var),
        'read_noise': np.sqrt(read_noise_var),
        'sky_noise': np.sqrt(sky_noise_var),
        'total_noise': total_noise,
        'snr': snr,
        'snr_db': 20 * np.log10(snr) if snr > 0 else -np.inf
    }


def estimate_limiting_magnitude(
    snr_threshold: float = 5.0,
    exposure_time: float = EXPOSURE_TIME,
    num_frames: int = 1000,
    read_noise: float = 3.0
) -> float:
    """
    Estimate limiting magnitude for speckle detection.
    
    Parameters
    ----------
    snr_threshold : float
        Required SNR per frame
    exposure_time : float
        Single exposure time
    num_frames : int
        Number of frames to average
    read_noise : float
        Detector read noise
    
    Returns
    -------
    float
        Limiting magnitude
    """
    # SNR improves as sqrt(N) for averaged power spectra
    # Search for magnitude giving threshold SNR
    
    from scipy.optimize import brentq
    
    def snr_minus_threshold(mag):
        snr_dict = compute_snr(mag, exposure_time, read_noise)
        # Account for averaging benefit
        effective_snr = snr_dict['snr'] * np.sqrt(num_frames)
        return effective_snr - snr_threshold
    
    try:
        limiting_mag = brentq(snr_minus_threshold, -5, 20)
    except ValueError:
        limiting_mag = np.nan
    
    return limiting_mag


if __name__ == "__main__":
    from config import WAVELENGTH
    
    print("Speckle Camera Configuration")
    print("=" * 50)
    
    print(f"\nExposure time: {EXPOSURE_TIME*1e3:.1f} ms")
    
    for mag in [0, 5, 10, 12]:
        snr = compute_snr(mag)
        print(f"\nMagnitude {mag}:")
        print(f"  Signal: {snr['signal_electrons']:.0f} e-")
        print(f"  Photon noise: {snr['photon_noise']:.1f} e-")
        print(f"  Read noise: {snr['read_noise']:.1f} e-")
        print(f"  SNR: {snr['snr']:.1f}")
    
    print(f"\nLimiting magnitude (1000 frames): {estimate_limiting_magnitude():.1f}")
