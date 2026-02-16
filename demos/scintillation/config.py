"""
Configuration for Hale Telescope speckle imaging simulation.

Contains:
- Hale telescope parameters (5.08m primary with central obscuration)
- Hufnagel-Valley 5/7 atmospheric model for Cn2 profile
- Stellar magnitude to photon flux conversion
- Three-layer atmospheric model based on HV-57
"""

import numpy as np

# =============================================================================
# Hale Telescope Parameters
# =============================================================================

HALE_DIAMETER = 5.08  # meters - primary mirror diameter
HALE_OBSCURATION_RATIO = 0.36  # central obscuration ratio (secondary/primary)
# The Hale telescope has a 1.83m secondary, giving ~0.36 ratio

# =============================================================================
# Observational Parameters
# =============================================================================

WAVELENGTH = 500e-9  # meters (500 nm)
EXPOSURE_TIME = 5e-3  # seconds (5 ms speckle exposure)

# Reference for zero-point flux at 500nm (V-band approximate)
# Vega (mag 0) flux at 500nm ~ 3.64e10 photons/s/m^2/nm
ZERO_POINT_FLUX = 3.64e10  # photons/s/m^2/nm at mag=0

# Detector bandwidth (typical speckle filter)
BANDWIDTH = 20e-9  # 20 nm bandwidth in meters


def stellar_flux(magnitude: float) -> float:
    """
    Calculate photon flux from stellar magnitude.
    
    Parameters
    ----------
    magnitude : float
        Apparent magnitude of the star
    
    Returns
    -------
    float
        Photon flux in photons/s/m^2
    """
    # Convert bandwidth to nm for flux calculation
    bandwidth_nm = BANDWIDTH * 1e9
    flux = ZERO_POINT_FLUX * 10**(-0.4 * magnitude) * bandwidth_nm
    return flux


def photons_per_exposure(magnitude: float, 
                          diameter: float = HALE_DIAMETER,
                          obscuration: float = HALE_OBSCURATION_RATIO,
                          exposure: float = EXPOSURE_TIME) -> float:
    """
    Calculate total photons collected per exposure.
    
    Parameters
    ----------
    magnitude : float
        Apparent magnitude of the star
    diameter : float
        Primary mirror diameter in meters
    obscuration : float
        Central obscuration ratio
    exposure : float
        Exposure time in seconds
    
    Returns
    -------
    float
        Total photons collected
    """
    flux = stellar_flux(magnitude)
    collecting_area = np.pi * (diameter/2)**2 * (1 - obscuration**2)
    return flux * collecting_area * exposure


# =============================================================================
# Hufnagel-Valley 5/7 Atmospheric Model
# =============================================================================

def hufnagel_valley_cn2(h: float, v_rms: float = 21.0) -> float:
    """
    Hufnagel-Valley 5/7 Cn2 profile.
    
    The HV-5/7 model gives:
    - r0 = 5 cm at 500nm at zenith
    - isoplanatic angle = 7 microradians
    
    Parameters
    ----------
    h : float
        Altitude in meters
    v_rms : float
        RMS wind speed (typically 21 m/s for HV-5/7)
    
    Returns
    -------
    float
        Cn2 value in m^(-2/3)
    """
    # Convert to km for the standard formula
    h_km = h / 1000.0
    
    # HV-5/7 model (Hufnagel-Valley with A = 1.7e-14)
    # Cn2(h) = 5.94e-53 * (v/27)^2 * h^10 * exp(-h/1) 
    #        + 2.7e-16 * exp(-h/1.5) 
    #        + A * exp(-h/0.1)
    
    A = 1.7e-14  # Ground layer strength for HV-5/7
    
    term1 = 5.94e-53 * (v_rms / 27.0)**2 * (h_km * 1000)**10 * np.exp(-h_km)
    term2 = 2.7e-16 * np.exp(-h_km / 1.5)
    term3 = A * np.exp(-h_km / 0.1)
    
    return term1 + term2 + term3


def compute_r0_from_cn2(cn2: float, dh: float, wavelength: float = WAVELENGTH) -> float:
    """
    Compute r0 contribution from a single layer.
    
    Parameters
    ----------
    cn2 : float
        Cn2 value at the layer
    dh : float
        Effective thickness of the layer in meters
    wavelength : float
        Observing wavelength in meters
    
    Returns
    -------
    float
        r0 contribution in meters
    """
    # r0 = [0.423 * k^2 * integral(Cn2 dh)]^(-3/5)
    k = 2 * np.pi / wavelength
    integral = cn2 * dh
    r0 = (0.423 * k**2 * integral)**(-3/5)
    return r0


def zenith_correction(r0_zenith: float, zenith_angle_deg: float) -> float:
    """
    Correct r0 for zenith angle.
    
    Parameters
    ----------
    r0_zenith : float
        r0 at zenith in meters
    zenith_angle_deg : float
        Zenith angle in degrees
    
    Returns
    -------
    float
        Corrected r0 in meters
    """
    # r0(z) = r0(0) * cos(z)^(3/5)
    z_rad = np.radians(zenith_angle_deg)
    return r0_zenith * np.cos(z_rad)**(3/5)


# =============================================================================
# Three-Layer Atmospheric Model
# =============================================================================

# Layer heights (meters) and their vertical extents
# Ground layer: 0-500m
# Layer 1: 500-2000m (centered at 1km)
# Layer 2: 2000-5000m (centered at 3km)
LAYER_HEIGHTS = [0.0, 1000.0, 3000.0]  # Representative heights for each layer


def compute_layer_integrated_cn2(h_low: float, h_high: float) -> float:
    """
    Integrate Cn2 over a layer from h_low to h_high.
    
    Parameters
    ----------
    h_low : float
        Lower boundary in meters
    h_high : float
        Upper boundary in meters
    
    Returns
    -------
    float
        Integrated Cn2 in m^(1/3)
    """
    from scipy import integrate
    integral, _ = integrate.quad(hufnagel_valley_cn2, h_low, h_high)
    return integral


# Layer boundaries for integration
LAYER_BOUNDS = [(0, 500), (500, 2000), (2000, 5000)]

# Compute integrated Cn2 for each layer
_LAYER_INTEGRATED_CN2_RAW = [compute_layer_integrated_cn2(lo, hi) for lo, hi in LAYER_BOUNDS]

# Scale factor to simulate better seeing conditions
# To get r0 10x larger (good seeing ~50cm instead of ~5cm), scale Cn2 by 1/10^(5/3)
SEEING_SCALE_FACTOR = 10**(-5/3)  # ~0.0215
LAYER_INTEGRATED_CN2 = [cn2 * SEEING_SCALE_FACTOR for cn2 in _LAYER_INTEGRATED_CN2_RAW]

# Wind velocities based on Bufton wind model (typical for HV-5/7)
# Ground layer: typically slow
# Higher layers: jet stream influence
def bufton_wind_speed(h: float) -> float:
    """
    Bufton wind model for wind speed vs altitude.
    
    Parameters
    ----------
    h : float
        Altitude in meters
    
    Returns
    -------
    float
        Wind speed in m/s
    """
    h_km = h / 1000.0
    # Bufton model: V(h) = 5 + 30*exp(-((h-9.4)/4.8)^2)
    # Modified for lower altitudes
    v_ground = 5.0  # m/s
    v_tropo = 30.0 * np.exp(-((h_km - 9.4) / 4.8)**2)
    return v_ground + v_tropo


# Layer wind speeds (m/s)
LAYER_WIND_SPEEDS = [bufton_wind_speed(h) for h in LAYER_HEIGHTS]

# Wind directions (radians) - varied for realistic simulation
LAYER_WIND_DIRECTIONS = [0.0, np.pi/4, np.pi/2]  # 0°, 45°, 90°


def compute_layer_r0(layer_idx: int, wavelength: float = WAVELENGTH) -> float:
    """
    Compute r0 for a specific layer at zenith.
    
    Parameters
    ----------
    layer_idx : int
        Index of the layer (0, 1, or 2)
    wavelength : float
        Observing wavelength in meters
    
    Returns
    -------
    float
        r0 for that layer in meters
    """
    # LAYER_INTEGRATED_CN2 already contains integral(Cn2 dh) for each layer
    integrated_cn2 = LAYER_INTEGRATED_CN2[layer_idx]
    k = 2 * np.pi / wavelength
    r0 = (0.423 * k**2 * integrated_cn2)**(-3/5)
    return r0


def compute_total_r0(wavelength: float = WAVELENGTH) -> float:
    """
    Compute total r0 from all layers at zenith.
    
    The total r0 combines as: r0_total^(-5/3) = sum(r0_i^(-5/3))
    
    Parameters
    ----------
    wavelength : float
        Observing wavelength in meters
    
    Returns
    -------
    float
        Total r0 in meters
    """
    r0_sum = 0.0
    for i in range(len(LAYER_HEIGHTS)):
        r0_i = compute_layer_r0(i, wavelength)
        r0_sum += r0_i**(-5/3)
    return r0_sum**(-3/5)


def get_layer_weights() -> np.ndarray:
    """
    Compute relative Cn2 weights for each layer.
    
    Returns
    -------
    np.ndarray
        Normalized weights summing to 1
    """
    # LAYER_INTEGRATED_CN2 already contains integral(Cn2 dh) for each layer
    weights = np.array(LAYER_INTEGRATED_CN2)
    return weights / weights.sum()


# Precompute layer r0 values at 500nm zenith
LAYER_R0_ZENITH = [compute_layer_r0(i) for i in range(len(LAYER_HEIGHTS))]
TOTAL_R0_ZENITH = compute_total_r0()

# Print configuration summary
if __name__ == "__main__":
    print("=" * 60)
    print("Hale Telescope Speckle Imaging Configuration")
    print("=" * 60)
    print(f"\nTelescope:")
    print(f"  Diameter: {HALE_DIAMETER} m")
    print(f"  Obscuration ratio: {HALE_OBSCURATION_RATIO}")
    print(f"  Collecting area: {np.pi * (HALE_DIAMETER/2)**2 * (1 - HALE_OBSCURATION_RATIO**2):.2f} m^2")
    
    print(f"\nObservation:")
    print(f"  Wavelength: {WAVELENGTH*1e9:.0f} nm")
    print(f"  Exposure time: {EXPOSURE_TIME*1e3:.1f} ms")
    print(f"  Bandwidth: {BANDWIDTH*1e9:.0f} nm")
    
    print(f"\nAtmospheric Layers (HV-5/7):")
    print(f"  {'Layer':<10} {'Height':<10} {'Cn2':<15} {'r0':<10} {'Wind':<10} {'Dir':<10}")
    for i, (h, cn2, r0, v, d) in enumerate(zip(LAYER_HEIGHTS, LAYER_CN2, 
                                                 LAYER_R0_ZENITH, LAYER_WIND_SPEEDS,
                                                 LAYER_WIND_DIRECTIONS)):
        print(f"  {i:<10} {h:<10.0f} {cn2:<15.2e} {r0*100:<10.2f}cm {v:<10.1f} {np.degrees(d):<10.0f}°")
    
    print(f"\n  Total r0 at zenith: {TOTAL_R0_ZENITH*100:.2f} cm")
    
    mag = 5.0
    print(f"\nStellar flux (mag={mag}):")
    print(f"  Flux: {stellar_flux(mag):.2e} photons/s/m^2")
    print(f"  Photons per exposure: {photons_per_exposure(mag):.2e}")
