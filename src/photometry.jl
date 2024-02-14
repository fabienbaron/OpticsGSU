function flux(mag; tel_surf=pi*(3.6/2)^2, airmass = 1.0, exptime=20e-3, filter="V")
  # Return number of incoming photons arriving at the camera level (= before detection)
  nphotons,extinct_coeff = mag_zeropoint(filter=filter);
  return exptime*tel_surf*10.0^(-0.4*mag)*nphotons*10.0^(-0.4*airmass*extinct_coeff);
end

function mag_zeropoint(;filter="none")
# Photons per square meter per second produced by a 0th mag star above the atmosphere.
# Assuming spectrum like Vega
    nphot = -1.0;
    coeff = 1e99; # extinction
    if (filter == "none") 
        nphot = 4.32e+10;
        coeff = 0.20;
    elseif (filter == "U")
        nphot = 5.50e+9;
        coeff = 0.60;
    elseif (filter == "B")
        nphot = 3.91e+9;
        coeff = 0.40;
    elseif (filter == "V")
        nphot = 8.66e+9;
        coeff = 0.20;
    elseif (filter == "R")
        nphot = 1.10e+10;
        coeff = 0.10;
    elseif (filter == "I")
        nphot = 6.75e+9;
        coeff = 0.08;
    end
  return nphot, coeff;
end



