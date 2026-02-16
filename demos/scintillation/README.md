# Hale Telescope Speckle Imaging Simulation

HCIpy-based simulation of speckle imaging with the Palomar 200-inch Hale Telescope through a realistic three-layer atmospheric model.

## Overview

This package simulates short-exposure "speckle" imaging of a star of magnitude M at 500nm, observing at zenith angle z with the Hale telescope. The atmosphere follows the Hufnagel-Valley 5/7 model with Taylor frozen flow.

## Physical Model

### Hale Telescope (Palomar 200-inch)
- **Primary diameter**: 5.08 m
- **Central obscuration**: 1.83 m (36% ratio)
- **Spider vanes**: 4 vanes at 90° spacing
- **Diffraction limit at 500nm**: ~20 mas (λ/D)

### Atmospheric Model (HV-5/7)

Three layers following the Hufnagel-Valley 5/7 Cn² profile:

| Layer | Height | Cn² (m⁻²/³) | r₀ at zenith | Wind Speed | Direction |
|-------|--------|-------------|--------------|------------|-----------|
| Ground | 0 m | ~1.7×10⁻¹⁴ | ~5 cm | 5 m/s | 0° |
| Layer 1 | 1 km | ~2.7×10⁻¹⁶ | ~20 cm | 7 m/s | 45° |
| Layer 2 | 3 km | ~1.5×10⁻¹⁶ | ~30 cm | 12 m/s | 90° |

The HV-5/7 model gives:
- r₀ ≈ 5 cm at 500nm at zenith
- Isoplanatic angle θ₀ ≈ 7 μrad

### Taylor Frozen Flow
Each layer evolves according to Taylor's frozen flow hypothesis:
- Phase screens translate with the wind velocity
- No temporal decorrelation within the layer
- Atmospheric coherence time τ₀ ~ 0.314 × r₀/v_eff

### Zenith Angle Correction
For observation at zenith angle z:
- r₀(z) = r₀(0) × cos(z)^(3/5)
- Effective path length increases as sec(z)

### Stellar Photometry
Photon flux from magnitude M at 500nm:
- Zero-point: 3.64×10¹⁰ photons/s/m²/nm (Vega)
- Bandwidth: 20 nm (typical speckle filter)
- Photons collected per exposure scales with telescope area and exposure time

### Speckle Camera
- **Exposure time**: 5 ms (shorter than τ₀)
- **Quantum efficiency**: 90%
- **Read noise**: 3 e⁻ RMS
- **Dark current**: 0.01 e⁻/pixel/s

## Installation

```bash
pip install hcipy numpy matplotlib scipy
```

## Usage

### Command Line

```bash
# Basic simulation
python simulate.py --magnitude 5 --zenith 30 --nframes 100

# Full options
python simulate.py \
    --magnitude 5 \      # Stellar magnitude
    --zenith 30 \        # Zenith angle (degrees)
    --nframes 100 \      # Number of frames
    --npix 256 \         # Pupil resolution
    --seed 42 \          # Random seed
    --output ./results   # Output directory
```

### Python API

```python
from hale_speckle import create_simulation, simulate_speckle_sequence

# Create simulation with M=5 star at z=30°
sim = create_simulation(
    magnitude=5.0,
    zenith_angle_deg=30.0,
    num_pixels=256,
    seed=42
)

print(f"r0 = {sim['r0']*100:.1f} cm")
print(f"Seeing = {sim['seeing']:.2f} arcsec")
print(f"Coherence time = {sim['coherence_time']*1e3:.1f} ms")

# Generate speckle sequence
data_cube = simulate_speckle_sequence(
    sim,
    num_frames=100,
    add_noise=True
)

# Get average power spectrum (for speckle interferometry)
ps = data_cube.power_spectrum_average()

# Get mean (long-exposure) image
mean_img = data_cube.mean_image()
```

## Output

The simulation produces:
1. **Speckle data cube**: Stack of short-exposure images
2. **Mean image**: Long-exposure equivalent
3. **Power spectrum average**: For speckle imaging reconstruction
4. **Autocorrelation**: For Knox-Thompson analysis
5. **Strehl ratios**: Mean and per-frame

## Module Structure

```
hale_speckle/
├── __init__.py      # Package exports
├── config.py        # Constants and configuration
├── atmosphere.py    # HV-5/7 atmospheric model
├── telescope.py     # Hale telescope optics
├── camera.py        # Speckle camera simulation
├── simulate.py      # Main simulation script
└── README.md        # This file
```

## Key Functions

### config.py
- `stellar_flux(magnitude)`: Convert magnitude to photon flux
- `photons_per_exposure(magnitude)`: Total photons collected
- `hufnagel_valley_cn2(height)`: HV-5/7 Cn² profile
- `zenith_correction(r0, zenith)`: Airmass correction

### atmosphere.py
- `create_frozen_flow_atmosphere()`: Build multi-layer atmosphere
- `compute_seeing(r0)`: Seeing FWHM from r₀
- `compute_coherence_time(r0, v_eff)`: Atmospheric τ₀

### telescope.py
- `HaleTelescope`: Complete optical model class
- `create_hale_pupil()`: Pupil with obscuration and spider
- `compute_diffraction_limit()`: λ/D metrics

### camera.py
- `SpeckleCamera`: Camera simulator with noise model
- `SpeckleDataCube`: Container for image sequences
- `compute_snr()`: Signal-to-noise calculation

## References

1. Fried, D.L. (1966). "Optical Resolution Through a Randomly Inhomogeneous Medium"
2. Hufnagel, R.E. & Stanley, N.R. (1964). "Modulation Transfer Function Associated with Image Transmission through Turbulent Media"
3. Valley, G.C. (1980). "Isoplanatic degradation of tilt correction and short-term imaging systems"
4. Taylor, G.I. (1938). "The Spectrum of Turbulence"
5. Labeyrie, A. (1970). "Attainment of Diffraction Limited Resolution in Large Telescopes by Fourier Analysing Speckle Patterns"

## License

MIT License
