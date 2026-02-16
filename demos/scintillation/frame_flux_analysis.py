"""
Per-Frame Flux Analysis

Simulates 100 speckle frames, measures flux in each frame,
plots the flux time series, and saves the image cube as FITS.
"""

import numpy as np
import hcipy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
from pathlib import Path
from typing import Optional

from config import (
    WAVELENGTH,
    EXPOSURE_TIME,
    HALE_DIAMETER,
    TOTAL_R0_ZENITH,
    photons_per_exposure,
    zenith_correction
)
from atmosphere import create_frozen_flow_atmosphere
from telescope import HaleTelescope, create_wavefront
from camera import SpeckleCamera


def create_atmosphere_with_scintillation(
    pupil_grid: hp.Grid,
    zenith_angle_deg: float = 0.0,
    outer_scale: float = 25.0,
    seed: Optional[int] = None
) -> hp.MultiLayerAtmosphere:
    """Create three-layer atmosphere with scintillation enabled."""
    return create_frozen_flow_atmosphere(
        pupil_grid,
        zenith_angle_deg=zenith_angle_deg,
        outer_scale=outer_scale,
        scintillation=True,
        seed=seed
    )


def aperture_photometry(image: hp.Field, aperture_radius_pix: int) -> float:
    """Simple circular aperture photometry."""
    img_2d = image.shaped
    ny, nx = img_2d.shape
    cy, cx = ny // 2, nx // 2
    
    y, x = np.ogrid[:ny, :nx]
    r = np.sqrt((x - cx)**2 + (y - cy)**2)
    mask = r <= aperture_radius_pix
    
    return img_2d[mask].sum()


def simulate_frame_sequence(
    magnitude: float = 5.0,
    zenith_angle_deg: float = 30.0,
    num_frames: int = 100,
    seed: Optional[int] = 42
) -> tuple:
    """
    Simulate a sequence of speckle frames.
    
    Returns
    -------
    tuple
        (image_cube, fluxes, times, telescope, focal_grid)
    """
    # Setup
    telescope = HaleTelescope(
        num_pixels=256,
        wavelength=WAVELENGTH,
        focal_sampling=4,
        include_spider=True
    )
    
    camera = SpeckleCamera(
        focal_grid=telescope.focal_grid,
        exposure_time=EXPOSURE_TIME
    )
    
    atmosphere = create_atmosphere_with_scintillation(
        telescope.pupil_grid,
        zenith_angle_deg=zenith_angle_deg,
        seed=seed
    )
    
    photons_per_frame = photons_per_exposure(magnitude)
    
    # Aperture: 3 lambda/D
    # With physical units: focal_grid.delta is in radians
    lambda_over_D = WAVELENGTH / HALE_DIAMETER  # radians
    pixel_scale_rad = telescope.focal_grid.delta[0]  # radians per pixel
    pixels_per_lambda_D = lambda_over_D / pixel_scale_rad
    aperture_radius_pix = int(3.0 * pixels_per_lambda_D)
    
    # Get image dimensions
    sample_img = telescope.reference_psf.shaped
    ny, nx = sample_img.shape
    
    # Storage
    image_cube = np.zeros((num_frames, ny, nx))
    fluxes = np.zeros(num_frames)
    times = np.zeros(num_frames)
    
    print(f"Simulating {num_frames} frames...")
    print(f"  Magnitude: {magnitude}")
    print(f"  Zenith angle: {zenith_angle_deg}°")
    print(f"  Photons/frame: {photons_per_frame:.2e}")
    
    for i in range(num_frames):
        t = i * EXPOSURE_TIME
        atmosphere.evolve_until(t)
        
        wf = create_wavefront(telescope.pupil_grid, WAVELENGTH)
        wf = atmosphere(wf)
        focal_wf = telescope.propagate(wf)
        psf = focal_wf.intensity
        
        image = camera.expose(psf, photons_per_frame, add_noise=True)
        
        image_cube[i] = image.shaped
        fluxes[i] = aperture_photometry(image, aperture_radius_pix)
        times[i] = t
        
        if (i + 1) % 20 == 0:
            print(f"  Frame {i+1}/{num_frames}")
    
    return image_cube, fluxes, times, telescope, aperture_radius_pix


def save_fits(image_cube: np.ndarray, 
              fluxes: np.ndarray,
              magnitude: float,
              zenith_angle_deg: float,
              output_path: str):
    """Save image cube and fluxes to FITS file."""
    
    # Primary HDU: image cube
    primary_hdu = fits.PrimaryHDU(image_cube.astype(np.float32))
    primary_hdu.header['NFRAMES'] = image_cube.shape[0]
    primary_hdu.header['EXPTIME'] = EXPOSURE_TIME
    primary_hdu.header['WAVELEN'] = WAVELENGTH
    primary_hdu.header['MAG'] = magnitude
    primary_hdu.header['ZENITH'] = zenith_angle_deg
    primary_hdu.header['R0'] = zenith_correction(TOTAL_R0_ZENITH, zenith_angle_deg)
    primary_hdu.header['TELESCOP'] = 'Hale 5.08m'
    primary_hdu.header['BUNIT'] = 'electrons'
    
    # Extension: flux measurements
    flux_hdu = fits.ImageHDU(fluxes.astype(np.float32), name='FLUXES')
    flux_hdu.header['TTYPE1'] = 'FLUX'
    flux_hdu.header['TUNIT1'] = 'electrons'
    
    hdul = fits.HDUList([primary_hdu, flux_hdu])
    hdul.writeto(output_path, overwrite=True)
    print(f"Saved FITS: {output_path}")


def plot_flux_timeseries(fluxes: np.ndarray, 
                         times: np.ndarray,
                         magnitude: float,
                         zenith_angle_deg: float,
                         output_path: str = None):
    """Plot flux vs frame number."""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Panel 1: Flux time series
    ax = axes[0, 0]
    frames = np.arange(len(fluxes))
    ax.plot(frames, fluxes, 'b-', linewidth=0.8, alpha=0.8)
    ax.axhline(np.mean(fluxes), color='r', linestyle='--', label=f'Mean: {np.mean(fluxes):.2e}')
    ax.fill_between(frames, 
                    np.mean(fluxes) - np.std(fluxes),
                    np.mean(fluxes) + np.std(fluxes),
                    alpha=0.2, color='r', label=f'±1σ: {np.std(fluxes):.2e}')
    ax.set_xlabel('Frame', fontsize=12)
    ax.set_ylabel('Flux (electrons)', fontsize=12)
    ax.set_title('Per-Frame Aperture Flux', fontsize=14)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    
    # Panel 2: Normalized flux (relative fluctuations)
    ax = axes[0, 1]
    rel_flux = fluxes / np.mean(fluxes)
    ax.plot(frames, rel_flux, 'b-', linewidth=0.8, alpha=0.8)
    ax.axhline(1.0, color='r', linestyle='--')
    ax.set_xlabel('Frame', fontsize=12)
    ax.set_ylabel('Relative Flux', fontsize=12)
    ax.set_title(f'Normalized Flux (σ/μ = {np.std(fluxes)/np.mean(fluxes)*100:.2f}%)', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    # Panel 3: Histogram
    ax = axes[1, 0]
    ax.hist(fluxes, bins=20, edgecolor='black', alpha=0.7)
    ax.axvline(np.mean(fluxes), color='r', linestyle='--', linewidth=2, label='Mean')
    ax.axvline(np.median(fluxes), color='g', linestyle=':', linewidth=2, label='Median')
    ax.set_xlabel('Flux (electrons)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Flux Distribution', fontsize=14)
    ax.legend()
    
    # Panel 4: Statistics
    ax = axes[1, 1]
    ax.axis('off')
    
    r0 = zenith_correction(TOTAL_R0_ZENITH, zenith_angle_deg)
    stats_text = f"""
    Simulation Parameters
    ─────────────────────
    Magnitude:        {magnitude}
    Zenith angle:     {zenith_angle_deg}°
    Wavelength:       {WAVELENGTH*1e9:.0f} nm
    Exposure time:    {EXPOSURE_TIME*1e3:.1f} ms
    r₀:               {r0*100:.1f} cm
    
    Flux Statistics
    ─────────────────────
    Frames:           {len(fluxes)}
    Mean flux:        {np.mean(fluxes):.4e} e⁻
    Std flux:         {np.std(fluxes):.4e} e⁻
    RMS variation:    {np.std(fluxes)/np.mean(fluxes)*100:.2f}%
    Min flux:         {np.min(fluxes):.4e} e⁻
    Max flux:         {np.max(fluxes):.4e} e⁻
    Peak-to-peak:     {(np.max(fluxes)-np.min(fluxes))/np.mean(fluxes)*100:.1f}%
    """
    ax.text(0.1, 0.9, stats_text, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.suptitle(f'Speckle Frame Flux Analysis (M={magnitude}, z={zenith_angle_deg}°)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved plot: {output_path}")
    
    plt.show()


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Per-frame flux analysis')
    parser.add_argument('-m', '--magnitude', type=float, default=5.0,
                        help='Stellar magnitude (default: 5)')
    parser.add_argument('-z', '--zenith', type=float, default=30.0,
                        help='Zenith angle in degrees (default: 30)')
    parser.add_argument('-n', '--nframes', type=int, default=100,
                        help='Number of frames (default: 100)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed (default: 42)')
    parser.add_argument('-o', '--output', type=str, default='speckle_sequence',
                        help='Output basename (default: speckle_sequence)')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("Per-Frame Flux Analysis")
    print("=" * 60)
    
    # Run simulation
    image_cube, fluxes, times, telescope, aperture_pix = simulate_frame_sequence(
        magnitude=args.magnitude,
        zenith_angle_deg=args.zenith,
        num_frames=args.nframes,
        seed=args.seed
    )
    
    # Save FITS
    fits_path = f"{args.output}_M{args.magnitude:.0f}_z{args.zenith:.0f}.fits"
    save_fits(image_cube, fluxes, args.magnitude, args.zenith, fits_path)
    
    # Plot
    plot_path = f"{args.output}_M{args.magnitude:.0f}_z{args.zenith:.0f}.png"
    plot_flux_timeseries(fluxes, times, args.magnitude, args.zenith, plot_path)
    
    print("\nDone!")


if __name__ == "__main__":
    main()
