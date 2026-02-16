"""
Photometric Accuracy vs Magnitude Analysis

Simulates speckle imaging photometry:
- 100-frame stacks averaged to get flux estimate
- Scintillation modeled by HCIpy (scintillation=True)
- Simple circular aperture photometry
- 20 realizations per magnitude to measure scatter
"""

import numpy as np
import hcipy as hp
import matplotlib.pyplot as plt
from typing import Optional, List
from dataclasses import dataclass
from tqdm import tqdm

from config import (
    WAVELENGTH,
    EXPOSURE_TIME,
    HALE_DIAMETER,
    TOTAL_R0_ZENITH,
    photons_per_exposure,
    zenith_correction
)
from atmosphere import (
    LAYER_HEIGHTS,
    LAYER_R0_ZENITH,
    LAYER_WIND_SPEEDS,
    LAYER_WIND_DIRECTIONS,
)
from telescope import HaleTelescope, create_wavefront
from camera import SpeckleCamera


@dataclass
class PhotometryResult:
    """Results for one magnitude."""
    magnitude: float
    flux_mean: float
    flux_std: float
    error_percent: float


def create_atmosphere_with_scintillation(
    pupil_grid: hp.Grid,
    zenith_angle_deg: float = 0.0,
    outer_scale: float = 25.0,
    seed: Optional[int] = None
) -> hp.MultiLayerAtmosphere:
    """
    Create three-layer atmosphere with scintillation enabled.
    """
    if seed is not None:
        np.random.seed(seed)
    
    sec_z = 1.0 / np.cos(np.radians(zenith_angle_deg))
    layers = []
    
    for height, r0_zenith, v_wind, theta_wind in zip(
            LAYER_HEIGHTS, LAYER_R0_ZENITH, LAYER_WIND_SPEEDS, LAYER_WIND_DIRECTIONS):
        
        r0 = zenith_correction(r0_zenith, zenith_angle_deg)
        velocity = v_wind * np.array([np.cos(theta_wind), np.sin(theta_wind)])
        eff_height = height * sec_z
        
        layer = hp.InfiniteAtmosphericLayer(
            pupil_grid,
            Cn_squared=1.0,
            L0=outer_scale,
            velocity=velocity,
            height=eff_height,
            use_interpolation=True
        )
        
        k = 2 * np.pi / WAVELENGTH
        layer.Cn_squared = r0**(-5/3) / (0.423 * k**2)
        layers.append(layer)
    
    # Scintillation ON
    return hp.MultiLayerAtmosphere(layers, scintillation=True)


def aperture_photometry(image: hp.Field, aperture_radius_pix: int) -> float:
    """
    Simple circular aperture photometry.
    
    Returns total flux (sum of electrons) within aperture.
    """
    img_2d = image.shaped
    ny, nx = img_2d.shape
    cy, cx = ny // 2, nx // 2
    
    y, x = np.ogrid[:ny, :nx]
    r = np.sqrt((x - cx)**2 + (y - cy)**2)
    mask = r <= aperture_radius_pix
    
    return img_2d[mask].sum()


def measure_photometry_single_realization(
    telescope: HaleTelescope,
    camera: SpeckleCamera,
    zenith_angle_deg: float,
    magnitude: float,
    num_frames: int,
    aperture_radius_pix: int,
    seed: int
) -> float:
    """
    Run one realization: simulate frames, average, measure flux.
    """
    # Fresh atmosphere
    atmosphere = create_atmosphere_with_scintillation(
        telescope.pupil_grid,
        zenith_angle_deg=zenith_angle_deg,
        seed=seed
    )
    
    photons_per_frame = photons_per_exposure(magnitude)
    
    # Accumulate frames
    stack = np.zeros(telescope.focal_grid.size)
    
    for i in range(num_frames):
        atmosphere.evolve_until(i * EXPOSURE_TIME)
        
        wf = create_wavefront(telescope.pupil_grid, WAVELENGTH)
        wf = atmosphere(wf)
        focal_wf = telescope.propagate(wf)
        psf = focal_wf.intensity
        
        image = camera.expose(psf, photons_per_frame, add_noise=True)
        stack += image
    
    # Mean image
    mean_image = hp.Field(stack / num_frames, telescope.focal_grid)
    
    # Aperture photometry
    flux = aperture_photometry(mean_image, aperture_radius_pix)
    
    return flux


def run_magnitude_sweep(
    magnitudes: List[float],
    zenith_angle_deg: float = 30.0,
    num_frames: int = 100,
    num_realizations: int = 20,
    base_seed: int = 42
) -> List[PhotometryResult]:
    """
    Run photometry analysis across magnitudes.
    """
    # Setup telescope and camera once
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
    
    # Aperture radius: 3 lambda/D
    pixel_scale = telescope.focal_grid.delta[0]
    aperture_radius_pix = int(3.0 / pixel_scale)
    
    print(f"Aperture radius: {aperture_radius_pix} pixels (~3 λ/D)")
    print(f"Frames per stack: {num_frames}")
    print(f"Realizations per magnitude: {num_realizations}")
    print()
    
    results = []
    
    for mag in tqdm(magnitudes, desc="Magnitudes"):
        fluxes = []
        
        for r in range(num_realizations):
            seed = base_seed + int(mag * 1000) + r
            flux = measure_photometry_single_realization(
                telescope, camera, zenith_angle_deg, mag,
                num_frames, aperture_radius_pix, seed
            )
            fluxes.append(flux)
        
        fluxes = np.array(fluxes)
        flux_mean = np.mean(fluxes)
        flux_std = np.std(fluxes)
        error_percent = 100 * flux_std / flux_mean if flux_mean > 0 else np.inf
        
        results.append(PhotometryResult(
            magnitude=mag,
            flux_mean=flux_mean,
            flux_std=flux_std,
            error_percent=error_percent
        ))
        
        tqdm.write(f"  M={mag:.1f}: flux={flux_mean:.2e} ± {flux_std:.2e}, error={error_percent:.2f}%")
    
    return results


def plot_results(results: List[PhotometryResult], output_path: str = None):
    """
    Plot photometric accuracy vs magnitude.
    """
    mags = [r.magnitude for r in results]
    errors = [r.error_percent for r in results]
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.semilogy(mags, errors, 'o-', color='C0', markersize=10, linewidth=2)
    
    # Reference lines
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.7, label='1%')
    ax.axhline(0.1, color='gray', linestyle=':', alpha=0.7, label='0.1%')
    
    ax.set_xlabel('Magnitude', fontsize=14)
    ax.set_ylabel('Photometric Error (%)', fontsize=14)
    ax.set_title('Photometric Accuracy vs Magnitude\n'
                 '(100-frame stack, 20 realizations, scintillation ON)', fontsize=14)
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(loc='upper left')
    
    # Info box
    r0 = zenith_correction(TOTAL_R0_ZENITH, 30.0)
    info = (f"Hale 5.08m, λ=500nm\n"
            f"z=30°, r₀={r0*100:.1f}cm\n"
            f"Exposure=5ms")
    ax.text(0.97, 0.03, info, transform=ax.transAxes, fontsize=10,
            ha='right', va='bottom', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"\nSaved: {output_path}")
    
    plt.show()


def main():
    print("=" * 60)
    print("Photometric Accuracy vs Magnitude")
    print("=" * 60)
    
    magnitudes = [0, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    
    results = run_magnitude_sweep(
        magnitudes=magnitudes,
        zenith_angle_deg=30.0,
        num_frames=100,
        num_realizations=20,
        base_seed=42
    )
    
    print("\n" + "=" * 60)
    print(f"{'Mag':>5} {'Flux':>12} {'Std':>12} {'Error%':>10}")
    print("-" * 60)
    for r in results:
        print(f"{r.magnitude:>5.1f} {r.flux_mean:>12.2e} {r.flux_std:>12.2e} {r.error_percent:>10.3f}")
    
    plot_results(results, output_path='photometric_accuracy.png')


if __name__ == "__main__":
    main()
