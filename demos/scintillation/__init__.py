"""
Hale Telescope Speckle Imaging Simulation Package

A HCIpy-based simulation of speckle imaging with the Palomar 200-inch
Hale telescope through a three-layer Hufnagel-Valley 5/7 atmosphere
following Taylor's frozen flow hypothesis.

Components:
-----------
config : Configuration and constants
    - Hale telescope parameters (5.08m, 36% obscuration)
    - HV-5/7 atmospheric Cn2 profile
    - Stellar magnitude to photon flux conversion

atmosphere : Atmospheric turbulence model
    - Three-layer model at 0m, 1km, 3km
    - r0 values from HV-5/7 profile
    - Taylor frozen flow with realistic wind profiles

telescope : Optical model
    - Hale pupil with central obscuration and spider vanes
    - Fraunhofer propagation to focal plane
    - Diffraction-limited PSF computation

camera : Speckle camera simulation
    - Short exposure (5ms) imaging
    - Photon noise (Poisson)
    - Read noise and dark current
    - Data cube handling for speckle processing

simulate : Main simulation script
    - End-to-end speckle sequence generation
    - Power spectrum and autocorrelation computation
    - Result visualization

Usage:
------
    python simulate.py --magnitude 5 --zenith 30 --nframes 100

Or programmatically:
    
    from hale_speckle import create_simulation, simulate_speckle_sequence
    
    sim = create_simulation(magnitude=5.0, zenith_angle_deg=30.0)
    data_cube = simulate_speckle_sequence(sim, num_frames=100)

Dependencies:
-------------
- hcipy >= 0.5
- numpy
- matplotlib
- scipy (optional, for some utilities)
"""

__version__ = '0.1.0'
__author__ = 'Fabien'

from .config import (
    HALE_DIAMETER,
    HALE_OBSCURATION_RATIO,
    WAVELENGTH,
    EXPOSURE_TIME,
    LAYER_HEIGHTS,
    TOTAL_R0_ZENITH,
    stellar_flux,
    photons_per_exposure,
    hufnagel_valley_cn2,
    zenith_correction
)

from .atmosphere import (
    create_frozen_flow_atmosphere,
    MultiLayerAtmosphere,
    compute_seeing,
    compute_coherence_time,
    get_effective_wind_speed
)

from .telescope import (
    HaleTelescope,
    create_hale_pupil,
    create_pupil_grid,
    compute_diffraction_limit
)

from .camera import (
    SpeckleCamera,
    SpeckleDataCube,
    compute_snr
)

from .simulate import (
    create_simulation,
    create_three_layer_atmosphere,
    simulate_speckle_sequence,
    compute_metrics,
    plot_results
)
