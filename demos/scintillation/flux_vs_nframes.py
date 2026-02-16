"""
Flux Estimate Accuracy vs Number of Frames

Compares measured aperture flux to ground truth as a function of
the number of frames averaged.

Ground truth = photons_per_exposure × QE × aperture_encircled_fraction
"""

import numpy as np
import hcipy as hp
import matplotlib.pyplot as plt
from typing import Optional, List
from dataclasses import dataclass

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
class AccuracyResult:
    """Results for one nframes value."""
    nframes: int
    measured_flux: float
    ground_truth: float
    bias_percent: float
    
    
def create_atmosphere_with_scintillation(
    pupil_grid: hp.Grid,
    zenith_angle_deg: float = 0.0,
    outer_scale: float = 25.0,
    seed: Optional[int] = None
) -> hp.MultiLayerAtmosphere:
    """Create three-layer atmosphere with scintillation enabled."""
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
    
    return hp.MultiLayerAtmosphere(layers, scintillation=True)


def aperture_photometry(image: hp.Field, aperture_radius_pix: int) -> float:
    """Simple circular aperture photometry."""
    img_2d = image.shaped
    ny, nx = img_2d.shape
    cy, cx = ny // 2, nx // 2
    
    y, x = np.ogrid[:ny, :nx]
    r = np.sqrt((x - cx)**2 + (y - cy)**2)
    mask = r <= aperture_radius_pix
    
    return img_2d[mask].sum()


def compute_ground_truth(
    telescope: HaleTelescope,
    aperture_radius_pix: int,
    magnitude: float,
    nframes: int
) -> float:
    """
    Compute ground truth flux.
    
    Ground truth = photons × QE × encircled_energy_fraction × nframes
    
    We measure encircled energy from the diffraction-limited PSF.
    """
    # Photons per frame
    photons = photons_per_exposure(magnitude)
    
    # QE
    qe = 0.9
    
    # Encircled energy in aperture (from reference PSF)
    ref_psf = telescope.reference_psf
    ref_psf_norm = ref_psf / ref_psf.sum()
    
    img_2d = ref_psf_norm.shaped
    ny, nx = img_2d.shape
    cy, cx = ny // 2, nx // 2
    
    y, x = np.ogrid[:ny, :nx]
    r = np.sqrt((x - cx)**2 + (y - cy)**2)
    mask = r <= aperture_radius_pix
    
    encircled_fraction = img_2d[mask].sum()
    
    # Ground truth per frame
    ground_truth_per_frame = photons * qe * encircled_fraction
    
    return ground_truth_per_frame * nframes


def simulate_and_measure(
    telescope: HaleTelescope,
    camera: SpeckleCamera,
    aperture_radius_pix: int,
    magnitude: float,
    zenith_angle_deg: float,
    nframes: int,
    seed: int
) -> float:
    """
    Simulate nframes, average, measure total flux.
    
    Returns measured flux (sum over mean image × nframes to get total).
    """
    atmosphere = create_atmosphere_with_scintillation(
        telescope.pupil_grid,
        zenith_angle_deg=zenith_angle_deg,
        seed=seed
    )
    
    photons_per_frame = photons_per_exposure(magnitude)
    
    # Accumulate
    stack = np.zeros(telescope.focal_grid.size)
    
    for i in range(nframes):
        atmosphere.evolve_until(i * EXPOSURE_TIME)
        
        wf = create_wavefront(telescope.pupil_grid, WAVELENGTH)
        wf = atmosphere(wf)
        focal_wf = telescope.propagate(wf)
        psf = focal_wf.intensity
        
        image = camera.expose(psf, photons_per_frame, add_noise=True)
        stack += image
    
    # Mean image
    mean_image = hp.Field(stack / nframes, telescope.focal_grid)
    
    # Measure flux on mean image, scale back to total
    flux_per_frame = aperture_photometry(mean_image, aperture_radius_pix)
    total_flux = flux_per_frame * nframes
    
    return total_flux


def run_nframes_analysis(
    magnitude: float = 5.0,
    zenith_angle_deg: float = 30.0,
    nframes_list: List[int] = [10, 50, 100, 200, 500, 1000],
    num_realizations: int = 10,
    base_seed: int = 42
) -> List[dict]:
    """
    Run analysis for different nframes values.
    
    Returns list of dicts with results for each nframes.
    """
    # Setup telescope and camera
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
    
    # Aperture: 3 lambda/D
    # Aperture: 3 lambda/D
    # With physical units: focal_grid.delta is in radians
    lambda_over_D = WAVELENGTH / HALE_DIAMETER  # radians
    pixel_scale_rad = telescope.focal_grid.delta[0]  # radians per pixel
    pixels_per_lambda_D = lambda_over_D / pixel_scale_rad
    aperture_radius_pix = int(3.0 * pixels_per_lambda_D)
    
    print(f"Aperture radius: {aperture_radius_pix} pixels (~3 λ/D)")
    print(f"Magnitude: {magnitude}")
    print(f"Zenith: {zenith_angle_deg}°")
    print(f"Realizations per nframes: {num_realizations}")
    print()
    
    results = []
    
    for nframes in nframes_list:
        print(f"nframes = {nframes}...")
        
        # Ground truth
        ground_truth = compute_ground_truth(
            telescope, aperture_radius_pix, magnitude, nframes
        )
        
        # Multiple realizations
        measured_fluxes = []
        for r in range(num_realizations):
            seed = base_seed + nframes * 100 + r
            flux = simulate_and_measure(
                telescope, camera, aperture_radius_pix,
                magnitude, zenith_angle_deg, nframes, seed
            )
            measured_fluxes.append(flux)
        
        measured_fluxes = np.array(measured_fluxes)
        mean_measured = np.mean(measured_fluxes)
        std_measured = np.std(measured_fluxes)
        
        bias = (mean_measured - ground_truth) / ground_truth * 100
        scatter = std_measured / ground_truth * 100
        
        results.append({
            'nframes': nframes,
            'ground_truth': ground_truth,
            'measured_mean': mean_measured,
            'measured_std': std_measured,
            'bias_percent': bias,
            'scatter_percent': scatter,
            'all_measurements': measured_fluxes
        })
        
        print(f"  Ground truth: {ground_truth:.4e}")
        print(f"  Measured:     {mean_measured:.4e} ± {std_measured:.4e}")
        print(f"  Bias:         {bias:+.2f}%")
        print(f"  Scatter:      {scatter:.2f}%")
        print()
    
    return results


def plot_results(results: List[dict], 
                 magnitude: float,
                 zenith_angle_deg: float,
                 output_path: str = None):
    """Plot accuracy vs nframes."""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    nframes = [r['nframes'] for r in results]
    ground_truth = [r['ground_truth'] for r in results]
    measured = [r['measured_mean'] for r in results]
    measured_std = [r['measured_std'] for r in results]
    bias = [r['bias_percent'] for r in results]
    scatter = [r['scatter_percent'] for r in results]
    
    # Panel 1: Measured vs Ground Truth
    ax = axes[0, 0]
    ax.errorbar(nframes, measured, yerr=measured_std, fmt='o-', 
                capsize=5, markersize=8, label='Measured', color='C0')
    ax.plot(nframes, ground_truth, 's--', markersize=8, 
            label='Ground Truth', color='C1')
    ax.set_xlabel('Number of Frames', fontsize=12)
    ax.set_ylabel('Total Flux (electrons)', fontsize=12)
    ax.set_title('Measured vs Ground Truth Flux', fontsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Panel 2: Bias
    ax = axes[0, 1]
    ax.axhline(0, color='gray', linestyle='--', alpha=0.7)
    ax.plot(nframes, bias, 'o-', markersize=10, color='C2', linewidth=2)
    ax.fill_between(nframes, 
                    [b - s for b, s in zip(bias, scatter)],
                    [b + s for b, s in zip(bias, scatter)],
                    alpha=0.2, color='C2')
    ax.set_xlabel('Number of Frames', fontsize=12)
    ax.set_ylabel('Bias (%)', fontsize=12)
    ax.set_title('Flux Estimate Bias', fontsize=14)
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    
    # Panel 3: Scatter (precision)
    ax = axes[1, 0]
    ax.plot(nframes, scatter, 'o-', markersize=10, color='C3', linewidth=2)
    
    # Theoretical sqrt(N) scaling
    scatter_0 = scatter[0] * np.sqrt(nframes[0])
    theoretical = [scatter_0 / np.sqrt(n) for n in nframes]
    ax.plot(nframes, theoretical, '--', color='gray', 
            label=r'$\propto 1/\sqrt{N}$', linewidth=2)
    
    ax.set_xlabel('Number of Frames', fontsize=12)
    ax.set_ylabel('Scatter (%)', fontsize=12)
    ax.set_title('Measurement Precision', fontsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Panel 4: Ratio measured/truth
    ax = axes[1, 1]
    ratio = [m / g for m, g in zip(measured, ground_truth)]
    ratio_err = [s / g for s, g in zip(measured_std, ground_truth)]
    
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.7)
    ax.errorbar(nframes, ratio, yerr=ratio_err, fmt='o-',
                capsize=5, markersize=10, color='C0', linewidth=2)
    ax.set_xlabel('Number of Frames', fontsize=12)
    ax.set_ylabel('Measured / Ground Truth', fontsize=12)
    ax.set_title('Flux Recovery Ratio', fontsize=14)
    ax.set_xscale('log')
    ax.set_ylim(0.8, 1.2)
    ax.grid(True, alpha=0.3)
    
    # Add annotation with encircled energy info
    # Recompute for reference
    r0 = zenith_correction(TOTAL_R0_ZENITH, zenith_angle_deg)
    info = (f"M={magnitude}, z={zenith_angle_deg}°\n"
            f"r₀={r0*100:.1f}cm, λ=500nm\n"
            f"Aperture: 3λ/D")
    ax.text(0.95, 0.05, info, transform=ax.transAxes, fontsize=10,
            ha='right', va='bottom', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.suptitle(f'Flux Estimate Accuracy vs Number of Frames',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}")
    
    plt.show()


def print_summary_table(results: List[dict]):
    """Print summary table."""
    print("\n" + "=" * 80)
    print("Summary: Measured Flux vs Ground Truth")
    print("=" * 80)
    print(f"{'N frames':>10} {'Ground Truth':>14} {'Measured':>14} {'Bias':>10} {'Scatter':>10}")
    print("-" * 80)
    for r in results:
        print(f"{r['nframes']:>10} {r['ground_truth']:>14.4e} {r['measured_mean']:>14.4e} "
              f"{r['bias_percent']:>+10.2f}% {r['scatter_percent']:>9.2f}%")
    print("=" * 80)


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Flux accuracy vs number of frames'
    )
    parser.add_argument('-m', '--magnitude', type=float, default=5.0,
                        help='Stellar magnitude (default: 5)')
    parser.add_argument('-z', '--zenith', type=float, default=30.0,
                        help='Zenith angle (default: 30)')
    parser.add_argument('-r', '--realizations', type=int, default=10,
                        help='Realizations per nframes (default: 10)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed (default: 42)')
    parser.add_argument('-o', '--output', type=str, default='flux_vs_nframes.png',
                        help='Output plot path')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("Flux Accuracy vs Number of Frames")
    print("=" * 60)
    
    nframes_list = [10, 50, 100, 200, 500, 1000]
    
    results = run_nframes_analysis(
        magnitude=args.magnitude,
        zenith_angle_deg=args.zenith,
        nframes_list=nframes_list,
        num_realizations=args.realizations,
        base_seed=args.seed
    )
    
    print_summary_table(results)
    
    plot_results(results, args.magnitude, args.zenith, args.output)


if __name__ == "__main__":
    main()
