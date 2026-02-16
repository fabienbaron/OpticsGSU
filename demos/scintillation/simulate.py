"""
Hale Telescope Speckle Imaging Simulation

Main simulation script combining:
- Three-layer HV-5/7 atmospheric model with Taylor frozen flow
- Hale telescope optical model
- Speckle camera with 5ms exposures

Usage:
    python simulate.py --magnitude 5 --zenith 30 --nframes 100

Author: Fabien
"""

import numpy as np
import hcipy as hp
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Optional, Tuple
import argparse

from config import (
    WAVELENGTH,
    EXPOSURE_TIME,
    HALE_DIAMETER,
    TOTAL_R0_ZENITH,
    photons_per_exposure
)
from atmosphere import (
    create_frozen_flow_atmosphere,
    MultiLayerAtmosphere,
    compute_seeing,
    compute_coherence_time,
    get_effective_wind_speed,
    zenith_correction,
    LAYER_HEIGHTS,
    LAYER_R0_ZENITH,
    LAYER_WIND_SPEEDS,
    LAYER_WIND_DIRECTIONS
)
from telescope import HaleTelescope, create_wavefront
from camera import SpeckleCamera, SpeckleDataCube, compute_snr


def create_simulation(
    magnitude: float = 5.0,
    zenith_angle_deg: float = 0.0,
    num_pixels: int = 256,
    focal_sampling: int = 4,
    seed: Optional[int] = None
) -> dict:
    """
    Create all simulation components.
    
    Parameters
    ----------
    magnitude : float
        Stellar magnitude (M)
    zenith_angle_deg : float
        Zenith angle (z) in degrees
    num_pixels : int
        Pupil grid resolution
    focal_sampling : int
        Pixels per lambda/D
    seed : int, optional
        Random seed
    
    Returns
    -------
    dict
        Dictionary with all simulation components
    """
    if seed is not None:
        np.random.seed(seed)
    
    print(f"Creating simulation components...")
    print(f"  Magnitude: {magnitude}")
    print(f"  Zenith angle: {zenith_angle_deg}°")
    
    # Create telescope
    telescope = HaleTelescope(
        num_pixels=num_pixels,
        wavelength=WAVELENGTH,
        focal_sampling=focal_sampling,
        include_spider=True
    )
    
    # Create atmosphere
    atmosphere = create_three_layer_atmosphere(
        telescope.pupil_grid,
        zenith_angle_deg=zenith_angle_deg,
        seed=seed
    )
    
    # Create camera
    camera = SpeckleCamera(
        focal_grid=telescope.focal_grid,
        exposure_time=EXPOSURE_TIME
    )
    
    # Compute photons per exposure
    photons = photons_per_exposure(magnitude)
    
    # Corrected r0
    r0 = zenith_correction(TOTAL_R0_ZENITH, zenith_angle_deg)
    
    return {
        'telescope': telescope,
        'atmosphere': atmosphere,
        'camera': camera,
        'magnitude': magnitude,
        'zenith_angle': zenith_angle_deg,
        'photons_per_exposure': photons,
        'r0': r0,
        'seeing': compute_seeing(r0),
        'coherence_time': compute_coherence_time(r0, get_effective_wind_speed())
    }


def create_three_layer_atmosphere(
    pupil_grid: hp.Grid,
    zenith_angle_deg: float = 0.0,
    outer_scale: float = 25.0,
    seed: Optional[int] = None
) -> hp.MultiLayerAtmosphere:
    """
    Create three-layer atmosphere using HCIpy.
    
    Parameters
    ----------
    pupil_grid : hp.Grid
        Pupil plane grid
    zenith_angle_deg : float
        Zenith angle in degrees
    outer_scale : float
        Outer scale (von Karman) in meters
    seed : int, optional
        Random seed
    
    Returns
    -------
    hp.MultiLayerAtmosphere
        HCIpy multi-layer atmosphere
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Airmass correction
    sec_z = 1.0 / np.cos(np.radians(zenith_angle_deg))
    
    # Create layers
    layers = []
    
    for i, (height, r0_zenith, v_wind, theta_wind) in enumerate(zip(
            LAYER_HEIGHTS, LAYER_R0_ZENITH, LAYER_WIND_SPEEDS, LAYER_WIND_DIRECTIONS)):
        
        # Correct r0 for zenith angle
        r0 = zenith_correction(r0_zenith, zenith_angle_deg)
        
        # Velocity vector
        velocity = v_wind * np.array([np.cos(theta_wind), np.sin(theta_wind)])
        
        # Effective height
        eff_height = height * sec_z
        
        # Create infinite atmospheric layer
        layer = hp.InfiniteAtmosphericLayer(
            pupil_grid,
            Cn_squared=1.0,  # Will be overridden
            L0=outer_scale,
            velocity=velocity,
            height=eff_height,
            use_interpolation=True
        )
        
        # Set proper Cn2 from r0
        # For a single layer with unit path length:
        # r0 = (0.423 * k^2 * Cn2)^(-3/5)
        # Cn2 = r0^(-5/3) / (0.423 * k^2)
        k = 2 * np.pi / WAVELENGTH
        cn2 = r0**(-5/3) / (0.423 * k**2)
        layer.Cn_squared = cn2
        
        layers.append(layer)
        
        print(f"  Layer {i}: h={height}m, r0={r0*100:.1f}cm, v={v_wind:.1f}m/s")
    
    # Create multi-layer atmosphere
    atmosphere = hp.MultiLayerAtmosphere(layers, scintillation=False)
    
    return atmosphere


def simulate_speckle_sequence(
    sim: dict,
    num_frames: int = 100,
    add_noise: bool = True,
    seed: Optional[int] = None
) -> SpeckleDataCube:
    """
    Simulate a sequence of speckle images.
    
    Parameters
    ----------
    sim : dict
        Simulation components from create_simulation()
    num_frames : int
        Number of frames to simulate
    add_noise : bool
        Add photon and read noise
    seed : int, optional
        Random seed
    
    Returns
    -------
    SpeckleDataCube
        Data cube with all frames
    """
    if seed is not None:
        np.random.seed(seed)
    
    telescope = sim['telescope']
    atmosphere = sim['atmosphere']
    camera = sim['camera']
    photons = sim['photons_per_exposure']
    
    print(f"\nSimulating {num_frames} speckle frames...")
    print(f"  Exposure time: {EXPOSURE_TIME*1e3:.1f} ms")
    print(f"  Photons per frame: {photons:.2e}")
    
    # Time step (exposure time)
    dt = EXPOSURE_TIME
    
    images = []
    
    for i in range(num_frames):
        # Evolve atmosphere
        t = i * dt
        atmosphere.evolve_until(t)
        
        # Create wavefront
        wf = create_wavefront(telescope.pupil_grid, WAVELENGTH)
        
        # Apply atmospheric phase
        wf = atmosphere(wf)
        
        # Propagate through telescope
        focal_wf = telescope.propagate(wf)
        
        # Get PSF
        psf = focal_wf.intensity
        
        # Expose with camera
        image = camera.expose(psf, photons, add_noise=add_noise)
        images.append(image)
        
        if (i + 1) % 10 == 0:
            print(f"  Frame {i+1}/{num_frames}")
    
    # Reset atmosphere for future use
    atmosphere.reset()
    
    return SpeckleDataCube(images, telescope.focal_grid, EXPOSURE_TIME)


def compute_metrics(
    data_cube: SpeckleDataCube,
    telescope: HaleTelescope
) -> dict:
    """
    Compute speckle imaging metrics.
    
    Parameters
    ----------
    data_cube : SpeckleDataCube
        Simulated data cube
    telescope : HaleTelescope
        Telescope model
    
    Returns
    -------
    dict
        Metrics dictionary
    """
    # Mean image
    mean_img = data_cube.mean_image()
    
    # Power spectrum
    ps = data_cube.power_spectrum_average()
    
    # Autocorrelation
    ac = data_cube.autocorrelation_average()
    
    # Strehl of mean image
    mean_strehl = telescope.compute_strehl(mean_img)
    
    # Single frame Strehl (estimate from first frame)
    frame_strehl = telescope.compute_strehl(data_cube.get_frame(0))
    
    return {
        'mean_image': mean_img,
        'power_spectrum': ps,
        'autocorrelation': ac,
        'mean_strehl': mean_strehl,
        'single_frame_strehl': frame_strehl,
        'num_frames': data_cube.num_frames
    }


def plot_results(
    sim: dict,
    data_cube: SpeckleDataCube,
    metrics: dict,
    output_dir: Optional[Path] = None
):
    """
    Plot simulation results.
    
    Parameters
    ----------
    sim : dict
        Simulation parameters
    data_cube : SpeckleDataCube
        Simulated data
    metrics : dict
        Computed metrics
    output_dir : Path, optional
        Directory to save plots
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Single speckle frame
    ax = axes[0, 0]
    frame = data_cube.get_frame(0)
    im = ax.imshow(frame.shaped, origin='lower', cmap='hot')
    ax.set_title('Single Speckle Frame')
    plt.colorbar(im, ax=ax, label='Electrons')
    
    # Mean image (long exposure)
    ax = axes[0, 1]
    mean = metrics['mean_image']
    im = ax.imshow(mean.shaped, origin='lower', cmap='hot')
    ax.set_title(f'Mean Image ({data_cube.num_frames} frames)')
    plt.colorbar(im, ax=ax, label='Electrons')
    
    # Power spectrum (log scale)
    ax = axes[0, 2]
    ps = metrics['power_spectrum']
    im = ax.imshow(np.log10(ps + 1), origin='lower', cmap='viridis')
    ax.set_title('Average Power Spectrum (log)')
    plt.colorbar(im, ax=ax)
    
    # Autocorrelation
    ax = axes[1, 0]
    ac = metrics['autocorrelation']
    ny, nx = ac.shape
    # Zoom to central region
    cx, cy = nx//2, ny//2
    hw = min(nx, ny) // 4
    ac_crop = ac[cy-hw:cy+hw, cx-hw:cx+hw]
    im = ax.imshow(ac_crop, origin='lower', cmap='RdBu_r')
    ax.set_title('Autocorrelation (center)')
    plt.colorbar(im, ax=ax)
    
    # Reference PSF (diffraction limited)
    ax = axes[1, 1]
    ref_psf = sim['telescope'].reference_psf
    im = ax.imshow(ref_psf.shaped, origin='lower', cmap='hot')
    ax.set_title('Diffraction-Limited PSF')
    plt.colorbar(im, ax=ax, label='Normalized')
    
    # Info panel
    ax = axes[1, 2]
    ax.axis('off')
    info_text = f"""
    Simulation Parameters:
    ─────────────────────
    Magnitude: {sim['magnitude']}
    Zenith angle: {sim['zenith_angle']}°
    Wavelength: {WAVELENGTH*1e9:.0f} nm
    
    Atmospheric Conditions:
    ─────────────────────
    r₀: {sim['r0']*100:.1f} cm
    Seeing: {sim['seeing']:.2f}"
    τ₀: {sim['coherence_time']*1e3:.1f} ms
    
    Camera:
    ─────────────────────
    Exposure: {EXPOSURE_TIME*1e3:.1f} ms
    Frames: {data_cube.num_frames}
    Photons/frame: {sim['photons_per_exposure']:.2e}
    
    Results:
    ─────────────────────
    Mean Strehl: {metrics['mean_strehl']:.4f}
    Frame Strehl: {metrics['single_frame_strehl']:.4f}
    """
    ax.text(0.1, 0.9, info_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='top',
            fontfamily='monospace')
    
    plt.suptitle(f'Hale Telescope Speckle Imaging (M={sim["magnitude"]}, z={sim["zenith_angle"]}°)',
                 fontsize=14)
    plt.tight_layout()
    
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_dir / 'speckle_results.png', dpi=150, bbox_inches='tight')
        print(f"Saved plot to {output_dir / 'speckle_results.png'}")
    
    plt.show()


def main():
    """Main simulation entry point."""
    parser = argparse.ArgumentParser(
        description='Hale Telescope Speckle Imaging Simulation'
    )
    parser.add_argument(
        '-m', '--magnitude', type=float, default=5.0,
        help='Stellar magnitude (default: 5.0)'
    )
    parser.add_argument(
        '-z', '--zenith', type=float, default=0.0,
        help='Zenith angle in degrees (default: 0.0)'
    )
    parser.add_argument(
        '-n', '--nframes', type=int, default=100,
        help='Number of frames (default: 100)'
    )
    parser.add_argument(
        '--npix', type=int, default=256,
        help='Pupil grid pixels (default: 256)'
    )
    parser.add_argument(
        '--no-noise', action='store_true',
        help='Disable photon and read noise'
    )
    parser.add_argument(
        '--seed', type=int, default=None,
        help='Random seed for reproducibility'
    )
    parser.add_argument(
        '-o', '--output', type=str, default='./output',
        help='Output directory (default: ./output)'
    )
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("Hale Telescope Speckle Imaging Simulation")
    print("=" * 60)
    
    # Create simulation
    sim = create_simulation(
        magnitude=args.magnitude,
        zenith_angle_deg=args.zenith,
        num_pixels=args.npix,
        seed=args.seed
    )
    
    # Print SNR estimate
    snr = compute_snr(args.magnitude)
    print(f"\nSNR estimate:")
    print(f"  Signal: {snr['signal_electrons']:.0f} e-")
    print(f"  SNR: {snr['snr']:.1f}")
    
    # Run simulation
    data_cube = simulate_speckle_sequence(
        sim,
        num_frames=args.nframes,
        add_noise=not args.no_noise,
        seed=args.seed
    )
    
    # Compute metrics
    print("\nComputing metrics...")
    metrics = compute_metrics(data_cube, sim['telescope'])
    
    print(f"\nResults:")
    print(f"  Mean Strehl ratio: {metrics['mean_strehl']:.4f}")
    print(f"  Single frame Strehl: {metrics['single_frame_strehl']:.4f}")
    
    # Plot results
    print("\nGenerating plots...")
    plot_results(sim, data_cube, metrics, output_dir=args.output)
    
    print("\nSimulation complete!")


if __name__ == "__main__":
    main()
