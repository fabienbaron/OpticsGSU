"""
Atmospheric turbulence model using HCIpy.

Implements a three-layer atmosphere following Taylor's frozen flow hypothesis
with Hufnagel-Valley 5/7 Cn2 profile (scaled for good seeing conditions).
"""

import numpy as np
import hcipy as hp
from typing import Optional

from config import (
    LAYER_HEIGHTS,
    LAYER_R0_ZENITH,
    LAYER_WIND_SPEEDS,
    LAYER_WIND_DIRECTIONS,
    TOTAL_R0_ZENITH,
    WAVELENGTH,
    zenith_correction
)


def create_frozen_flow_atmosphere(
    pupil_grid: hp.Grid,
    zenith_angle_deg: float = 0.0,
    wavelength: float = WAVELENGTH,
    outer_scale: float = 25.0,
    scintillation: bool = False,
    seed: Optional[int] = None
) -> hp.MultiLayerAtmosphere:
    """
    Create a multi-layer atmospheric model with Taylor frozen flow.
    
    Uses HCIpy's InfiniteAtmosphericLayer for each layer derived from 
    HV-5/7 profile (scaled for good seeing).
    
    Parameters
    ----------
    pupil_grid : hp.Grid
        The pupil plane grid
    zenith_angle_deg : float
        Zenith angle in degrees (affects r0 via airmass)
    wavelength : float
        Observing wavelength in meters
    outer_scale : float
        Outer scale of turbulence in meters (von Karman)
    scintillation : bool
        Enable scintillation modeling (Fresnel propagation)
    seed : int, optional
        Random seed for reproducibility
    
    Returns
    -------
    hp.MultiLayerAtmosphere
        Multi-layer atmospheric model
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Airmass correction factor
    sec_z = 1.0 / np.cos(np.radians(zenith_angle_deg))
    
    layers = []
    
    for i, (height, r0_zenith, v_wind, theta_wind) in enumerate(zip(
            LAYER_HEIGHTS, LAYER_R0_ZENITH, LAYER_WIND_SPEEDS, LAYER_WIND_DIRECTIONS)):
        
        # Correct r0 for zenith angle: r0(z) = r0(0) * cos(z)^(3/5)
        r0_layer = zenith_correction(r0_zenith, zenith_angle_deg)
        
        # Convert r0 to Cn2 using HCIpy's utility
        cn2 = hp.Cn_squared_from_fried_parameter(r0_layer, wavelength)
        
        # Wind velocity vector
        v_x = v_wind * np.cos(theta_wind)
        v_y = v_wind * np.sin(theta_wind)
        velocity = np.array([v_x, v_y])
        
        # Effective height (for scintillation, height matters)
        effective_height = height * sec_z
        
        # Create layer
        layer = hp.InfiniteAtmosphericLayer(
            pupil_grid,
            Cn_squared=cn2,
            L0=outer_scale,
            velocity=velocity,
            height=effective_height
        )
        
        layers.append(layer)
    
    # Create multi-layer atmosphere
    atmosphere = hp.MultiLayerAtmosphere(layers, scintillation=scintillation)
    
    return atmosphere


def compute_seeing(r0: float, wavelength: float = WAVELENGTH) -> float:
    """
    Compute seeing FWHM from r0.
    
    FWHM = 0.98 * lambda / r0 (radians)
    
    Parameters
    ----------
    r0 : float
        Fried parameter in meters
    wavelength : float
        Wavelength in meters
    
    Returns
    -------
    float
        Seeing FWHM in arcseconds
    """
    fwhm_rad = 0.98 * wavelength / r0
    fwhm_arcsec = np.degrees(fwhm_rad) * 3600
    return fwhm_arcsec


def compute_coherence_time(r0: float, v_eff: float) -> float:
    """
    Compute atmospheric coherence time.
    
    tau0 = 0.314 * r0 / v_eff
    
    Parameters
    ----------
    r0 : float
        Fried parameter in meters
    v_eff : float
        Effective wind speed in m/s
    
    Returns
    -------
    float
        Coherence time in seconds
    """
    return 0.314 * r0 / v_eff


def get_effective_wind_speed() -> float:
    """
    Compute effective wind speed from layer parameters.
    
    v_eff = [sum(Cn2_i * v_i^(5/3)) / sum(Cn2_i)]^(3/5)
    
    For simplicity, we use a weighted average based on r0 contributions.
    
    Returns
    -------
    float
        Effective wind speed in m/s
    """
    # Weight by turbulence strength (proportional to r0^(-5/3))
    weights = np.array([r0**(-5/3) for r0 in LAYER_R0_ZENITH])
    weights /= weights.sum()
    
    v_53 = np.array(LAYER_WIND_SPEEDS)**(5/3)
    v_eff = (np.sum(weights * v_53))**(3/5)
    
    return v_eff


def get_atmospheric_parameters(zenith_angle_deg: float = 0.0) -> dict:
    """
    Get summary of atmospheric parameters.
    
    Parameters
    ----------
    zenith_angle_deg : float
        Zenith angle in degrees
    
    Returns
    -------
    dict
        Dictionary of atmospheric parameters
    """
    r0 = zenith_correction(TOTAL_R0_ZENITH, zenith_angle_deg)
    v_eff = get_effective_wind_speed()
    
    return {
        'r0': r0,
        'r0_cm': r0 * 100,
        'seeing_arcsec': compute_seeing(r0),
        'tau0': compute_coherence_time(r0, v_eff),
        'v_eff': v_eff,
        'zenith_angle_deg': zenith_angle_deg,
        'airmass': 1.0 / np.cos(np.radians(zenith_angle_deg)),
        'layer_r0_cm': [zenith_correction(r0z, zenith_angle_deg) * 100 
                       for r0z in LAYER_R0_ZENITH],
        'layer_heights': LAYER_HEIGHTS,
        'layer_wind_speeds': LAYER_WIND_SPEEDS
    }


if __name__ == "__main__":
    # Print atmospheric parameters
    print("Atmospheric Model Summary")
    print("=" * 50)
    
    for z in [0, 30, 45, 60]:
        params = get_atmospheric_parameters(z)
        print(f"\nZenith angle: {z}°")
        print(f"  r0: {params['r0_cm']:.1f} cm")
        print(f"  Seeing: {params['seeing_arcsec']:.2f} arcsec")
        print(f"  τ0: {params['tau0']*1000:.1f} ms")
        print(f"  Effective wind: {params['v_eff']:.1f} m/s")
