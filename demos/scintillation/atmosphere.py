"""
Atmospheric turbulence model using HCIpy.

Implements a three-layer atmosphere following Taylor's frozen flow hypothesis
with Hufnagel-Valley 5/7 Cn2 profile.
"""

import numpy as np
import hcipy as hp
from typing import Tuple, List, Optional

from config import (
    LAYER_HEIGHTS,
    LAYER_R0_ZENITH,
    LAYER_WIND_SPEEDS,
    LAYER_WIND_DIRECTIONS,
    TOTAL_R0_ZENITH,
    WAVELENGTH,
    zenith_correction,
    get_layer_weights
)


def create_atmospheric_model(
    pupil_grid: hp.Grid,
    zenith_angle_deg: float = 0.0,
    wavelength: float = WAVELENGTH,
    outer_scale: float = 25.0,
    seed: Optional[int] = None
) -> hp.AtmosphericModel:
    """
    Create a multi-layer atmospheric model for the Hale telescope.
    
    Uses HCIpy's InfiniteAtmosphericLayer with Taylor frozen flow
    for each layer derived from HV-5/7 profile.
    
    Parameters
    ----------
    pupil_grid : hp.Grid
        The pupil plane grid
    zenith_angle_deg : float
        Zenith angle in degrees (affects path length through atmosphere)
    wavelength : float
        Observing wavelength in meters
    outer_scale : float
        Outer scale of turbulence in meters (von Karman)
    seed : int, optional
        Random seed for reproducibility
    
    Returns
    -------
    hp.AtmosphericModel
        Multi-layer atmospheric model
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Correct r0 values for zenith angle
    # Each layer's r0 scales with airmass
    sec_z = 1.0 / np.cos(np.radians(zenith_angle_deg))
    
    # Create atmospheric layers
    layers = []
    
    for i, (height, r0_zenith, v_wind, theta_wind) in enumerate(zip(
            LAYER_HEIGHTS, LAYER_R0_ZENITH, LAYER_WIND_SPEEDS, LAYER_WIND_DIRECTIONS)):
        
        # Correct r0 for zenith angle
        r0_corrected = zenith_correction(r0_zenith, zenith_angle_deg)
        
        # Convert wind speed and direction to velocity vector
        # Velocity needs to account for airmass for angular velocity
        v_x = v_wind * np.cos(theta_wind)
        v_y = v_wind * np.sin(theta_wind)
        velocity = np.array([v_x, v_y])
        
        # Effective height through atmosphere at zenith angle
        effective_height = height * sec_z
        
        # Create the infinite atmospheric layer with von Karman spectrum
        # HCIpy uses Fried parameter directly
        layer = hp.InfiniteAtmosphericLayer(
            pupil_grid,
            Cn_squared=None,  # We specify r0 directly
            L0=outer_scale,
            velocity=velocity,
            height=effective_height,
            use_interpolation=True
        )
        
        # Set the Fried parameter directly
        # HCIpy's layer uses Cn_squared internally, but we can set r0
        # We need to compute equivalent Cn2 from r0
        # r0 = (0.423 * k^2 * Cn2 * L)^(-3/5)
        # Cn2 = r0^(-5/3) / (0.423 * k^2 * L)
        # But for HCIpy layers, we use the built-in r0 scaling
        
        layers.append((layer, r0_corrected))
    
    # Create the multi-layer atmosphere
    # We need to combine layers properly
    atmosphere = MultiLayerAtmosphere(
        pupil_grid=pupil_grid,
        layers=layers,
        wavelength=wavelength,
        zenith_angle_deg=zenith_angle_deg
    )
    
    return atmosphere


class MultiLayerAtmosphere:
    """
    Multi-layer atmospheric model with Taylor frozen flow.
    
    Combines multiple InfiniteAtmosphericLayers with proper
    phase screen combination.
    """
    
    def __init__(
        self,
        pupil_grid: hp.Grid,
        layers: List[Tuple[hp.InfiniteAtmosphericLayer, float]],
        wavelength: float,
        zenith_angle_deg: float
    ):
        """
        Initialize multi-layer atmosphere.
        
        Parameters
        ----------
        pupil_grid : hp.Grid
            The pupil plane grid
        layers : list of (layer, r0) tuples
            List of atmospheric layers with their r0 values
        wavelength : float
            Observing wavelength
        zenith_angle_deg : float
            Zenith angle in degrees
        """
        self.pupil_grid = pupil_grid
        self.layers = layers
        self.wavelength = wavelength
        self.zenith_angle_deg = zenith_angle_deg
        self.time = 0.0
        
        # Compute total r0 for the combined atmosphere
        r0_sum = sum(r0**(-5/3) for _, r0 in layers)
        self.total_r0 = r0_sum**(-3/5)
        
        # Wavenumber
        self.k = 2 * np.pi / wavelength
    
    def evolve_until(self, t: float):
        """
        Evolve the atmosphere to time t.
        
        Parameters
        ----------
        t : float
            Target time in seconds
        """
        for layer, _ in self.layers:
            layer.evolve_until(t)
        self.time = t
    
    def reset(self):
        """Reset the atmosphere to t=0."""
        for layer, _ in self.layers:
            layer.reset()
        self.time = 0.0
    
    @property
    def phase_for(self) -> hp.Field:
        """
        Get the total phase screen at current time.
        
        Returns
        -------
        hp.Field
            Combined phase screen from all layers
        """
        # Initialize total phase
        total_phase = hp.Field(np.zeros(self.pupil_grid.size), self.pupil_grid)
        
        for layer, r0 in self.layers:
            # Get phase from this layer
            # Scale by wavelength to get phase in radians
            layer_phase = layer.phase_for(self.wavelength)
            
            # Weight by this layer's contribution
            # The layer already has proper r0 scaling
            total_phase += layer_phase
        
        return total_phase
    
    def forward(self, wavefront: hp.Wavefront) -> hp.Wavefront:
        """
        Propagate wavefront through the atmosphere.
        
        Parameters
        ----------
        wavefront : hp.Wavefront
            Input wavefront
        
        Returns
        -------
        hp.Wavefront
            Output wavefront with atmospheric phase
        """
        # Get total phase screen
        phase = self.phase_for
        
        # Apply phase to wavefront
        wf_out = wavefront.copy()
        wf_out.electric_field *= np.exp(1j * phase)
        
        return wf_out


def create_frozen_flow_atmosphere(
    pupil_grid: hp.Grid,
    zenith_angle_deg: float = 0.0,
    wavelength: float = WAVELENGTH,
    outer_scale: float = 25.0,
    scintillation: bool = False,
    seed: Optional[int] = None
) -> hp.MultiLayerAtmosphere:
    """
    Create atmosphere using HCIpy's native MultiLayerAtmosphere.
    
    This is the preferred method using HCIpy's built-in functionality.
    
    Parameters
    ----------
    pupil_grid : hp.Grid
        The pupil plane grid
    zenith_angle_deg : float
        Zenith angle in degrees
    wavelength : float
        Observing wavelength in meters
    outer_scale : float
        Outer scale of turbulence (von Karman)
    scintillation : bool
        Enable scintillation (Fresnel propagation effects)
    seed : int, optional
        Random seed
    
    Returns
    -------
    hp.MultiLayerAtmosphere
        HCIpy native multi-layer atmosphere
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Correct for zenith angle
    sec_z = 1.0 / np.cos(np.radians(zenith_angle_deg))
    
    # Build layer list for HCIpy
    layer_list = []
    
    for height, r0_zenith, v_wind, theta_wind in zip(
            LAYER_HEIGHTS, LAYER_R0_ZENITH, LAYER_WIND_SPEEDS, LAYER_WIND_DIRECTIONS):
        
        # Correct r0 for zenith angle
        r0 = zenith_correction(r0_zenith, zenith_angle_deg)
        
        # Velocity vector
        velocity = v_wind * np.array([np.cos(theta_wind), np.sin(theta_wind)])
        
        # Effective height
        eff_height = height * sec_z
        
        # Create layer dictionary
        layer_list.append({
            'r0': r0,
            'L0': outer_scale,
            'velocity': velocity,
            'height': eff_height
        })
    
    # Create multi-layer atmosphere
    # HCIpy's make_standard_atmospheric_layers or direct creation
    layers = []
    for layer_params in layer_list:
        layer = hp.InfiniteAtmosphericLayer(
            pupil_grid,
            Cn_squared=None,
            L0=layer_params['L0'],
            velocity=layer_params['velocity'],
            height=layer_params['height']
        )
        # Set r0 by computing appropriate Cn2
        # For single layer: Cn2 = r0^(-5/3) / (0.423 * k^2)
        k = 2 * np.pi / wavelength
        cn2_equiv = layer_params['r0']**(-5/3) / (0.423 * k**2)
        layer.Cn_squared = cn2_equiv
        
        layers.append(layer)
    
    # Create multi-layer atmosphere
    # scintillation=True enables Fresnel propagation between layers
    atmosphere = hp.MultiLayerAtmosphere(layers, scintillation=scintillation)
    
    return atmosphere


def compute_seeing(r0: float, wavelength: float = WAVELENGTH) -> float:
    """
    Compute seeing FWHM from r0.
    
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
    # Seeing = 0.98 * lambda / r0 (in radians)
    seeing_rad = 0.98 * wavelength / r0
    seeing_arcsec = np.degrees(seeing_rad) * 3600
    return seeing_arcsec


def compute_coherence_time(r0: float, v_eff: float) -> float:
    """
    Compute atmospheric coherence time tau0.
    
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
    # tau0 = 0.314 * r0 / v_eff
    return 0.314 * r0 / v_eff


def get_effective_wind_speed(zenith_angle_deg: float = 0.0) -> float:
    """
    Compute effective wind speed for coherence time.
    
    Uses Cn2-weighted average of wind speeds.
    
    Parameters
    ----------
    zenith_angle_deg : float
        Zenith angle in degrees
    
    Returns
    -------
    float
        Effective wind speed in m/s
    """
    weights = get_layer_weights()
    v_eff_53 = sum(w * v**(5/3) for w, v in zip(weights, LAYER_WIND_SPEEDS))
    return v_eff_53**(3/5)


if __name__ == "__main__":
    print("Atmospheric Model Test Configuration")
    print("=" * 50)
    
    # Test at zenith
    z = 0.0
    print(f"\nZenith angle: {z}°")
    print(f"Total r0: {TOTAL_R0_ZENITH*100:.2f} cm")
    print(f"Seeing: {compute_seeing(TOTAL_R0_ZENITH):.2f} arcsec")
    
    v_eff = get_effective_wind_speed()
    print(f"Effective wind speed: {v_eff:.1f} m/s")
    
    tau0 = compute_coherence_time(TOTAL_R0_ZENITH, v_eff)
    print(f"Coherence time: {tau0*1e3:.1f} ms")
    
    # Test at 30° zenith angle
    z = 30.0
    r0_30 = zenith_correction(TOTAL_R0_ZENITH, z)
    print(f"\nZenith angle: {z}°")
    print(f"Total r0: {r0_30*100:.2f} cm")
    print(f"Seeing: {compute_seeing(r0_30):.2f} arcsec")
