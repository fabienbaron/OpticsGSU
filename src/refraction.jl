#  function refract(λ; P=600,T=7.0,fwv=8.0)
#      # λ needs to be in microns
#      # Default dispersion values for observing site at ~2km altitude and +/- 30 degrees latitude
#      # P = 600 # air pressure in mm Hg
#      # T = 7.0 # air temperature in degrees C
#      # fwv = 8.0 # water vapor pressure in mm Hg
#      nl15760 = (64.328+(29498.1/(146-(1/λ)^2))+(255.4/(41-(1/λ)^2)))
#      nlTP = nl15760 * (P*(1+(1.049-0.0157*T)*1e-6*P)/(720.883*(1+0.003661*T)))
#      nlTP = nlTP-fwv*(0.06224-(0.000680/λ^2))/(1+0.003661*T)
#      return nlTP
#  end

function air_refractive_index_minus_one(wave; pressure=80.0, temperature=280.15, H2O_pressure=1.0665)
    """
    Adapted from GalSim to take meters 
    Return the refractive index of air as function of wavelength.

    Uses the formulae given in Filippenko (1982), which appear to come from Edlen (1953),
    and Coleman, Bozman, and Meggers (1960).  The units of the original formula are non-SI,
    being mmHg for pressure (and water vapor pressure), and degrees C for temperature.  This
    function accepts SI units, however, and transforms them when plugging into the formula.

    The default values for temperature, pressure and water vapor pressure are expected to be
    appropriate for LSST at Cerro Pachon, Chile, but they are broadly reasonable for most
    observatories.

    Parameters:
        wave:             Wavelength array in meters
        pressure:         Air pressure in kiloPascals.
        temperature:      Temperature in Kelvins.
        H2O_pressure:     Water vapor pressure in kiloPascals.

    Returns:
        the refractive index minus 1.
    """
    P = pressure * 7.50061683 # kPa -> mmHg
    T = temperature - 273.15 # K -> C
    W = H2O_pressure * 7.50061683 # kPa -> mmHg
    sigma_squared = 1.0 / (wave * 1e6)^2.0 # inverse wavenumber squared in micron^-2
    n_minus_one = (64.328 + (29498.1 / (146.0 - sigma_squared))+ (255.4 / (41.0 - sigma_squared))) * 1.e-6
    n_minus_one *= P * (1.0 + (1.049 - 0.0157 * T) * 1.e-6 * P) / (720.883 * (1.0 + 0.003661 * T))
    n_minus_one -= (0.0624 - 0.000680 * sigma_squared)/(1.0 + 0.003661 * T) * W * 1.e-6
    return n_minus_one
end

function get_refraction(wave, zenith_angle; pressure=69.328, temperature=293.15, H2O_pressure=1.067)
    """Compute the angle of refraction for a photon entering the atmosphere.

    Photons refract when transitioning from space, where the refractive index n = 1.0 exactly, to
    air, where the refractive index is slightly greater than 1.0.  This function computes the
    change in zenith angle for a photon with a given wavelength.  Output is a positive number of
    radians, even though the apparent zenith angle technically decreases due to this effect.

    Parameters:
        wave:            Wavelength array in nanometers
        zenith_angle:    `Angle` from object to zenith (radian)
        **kwargs:        Keyword arguments to pass to air_refractive_index() to override default
                         pressure, temperature, and/or H2O_pressure.

    Returns:
        the absolute value of change in zenith angle in radians.
    """
    nm1 = air_refractive_index_minus_one(wave, pressure = pressure, temperature = temperature, H2O_pressure=H2O_pressure)
    # The following line is equivalent to:
    # n_squared = (nm1 + 1)**2
    # r0 = (n_squared - 1.0) / (2.0 * n_squared)
    r0 = nm1 * (nm1+2) / 2.0 / (nm1^2 + 2*nm1 + 1) 
    #(α*(1-β)*tan(zenith_angle) - α*(β-0.5*α)*tan(zenith_angle)^3)
    return r0 * tan(zenith_angle)
end