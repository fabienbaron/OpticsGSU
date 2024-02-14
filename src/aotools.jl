using NumericalIntegration

function isoplanaticAngle(cn2, h; λ=500.E-9)
    # Calculates the isoplanatic angle from the Cn2 profile
    # Parameters:
    #     cn2 (array): Cn2 profile in m^2/3.
    #     h (Array): Altitude levels of cn2 profile in m
    #     λ : wavelength
    # Returns:
    #     isoplanatic angle in arcseconds
    Jh = sum(cn2.*(h.^(5.0/3.0)))
    iso = 0.057*λ^(6.0/5.0)*Jh^(-3.0/5.0)*180*3600/pi
    return iso
end

function CN2_Hufnagel(h, r0, λ, h0, h1, h2) # works for r0 , λ, h = floats or vectors
	c0power2 =  1.027e-3
	return c0power2*r0.^(-5/3).*(2pi./λ).^(-2).*( (h/h0).^10 .* exp.(-h/h1) + exp.(-h/h2))
end

function CN2_huffnagel_valley(h; V = 21.0, Cn2_0 = 1.7e-14)
    # h     in [m]
    # v     in [m/s]  :  rms upper (5-20km) atmospheric wind speed 
    # Cn2_0 in [m^-(2/3)] : CN2(h = 1m above ground) 
    # Default values yield a H-V 5/7 profile, with "typical" r0=5cm and θ0=7μrad (typical for λ = 500nm)
    z = h/1000; # [km]
    return 8.148e-26*V^2*(z.^10).*exp.(-h)+2.7e-16*exp.(-z/1.5)+Cn2_0*exp.(-10*z)
end

function CN2_huffnagel_valley_generalized(h; A=1.7e-14, hA=100, B=2.7e-16, hB=1500, C=3.59e-53, hC=1000, D=0, d=1, hD=0)
    # Everything in meters - default is HV 5-7
    # A is the coefficient for the surface (boundary layer) turbulence strength (m−2/3 )
    # hA is the height for its 1/e decay (meters)
    # B and hB similarly define the turbulence in the troposphere (up to about 10 km)
    # C and HC define the turbulence peak at the tropopause
    # D and HD define one or more isolated layers of turbulence, with d being the layer thickness (meters).
    return A*exp.(-h/hA) +  B*exp.(-h/hB) + C*(h.^10).*exp.(-h/hC)+ D*exp.(-(h.-hD).^2/(2*d^2))
end

function wind_profile_log(h; ustar = 1.5, κ = 0.41, d=0.5 , z0=0.5 )
 # h = altitude [m] -- Needs to be > d for this to work
 # ustar = friction velocity 
 # κ = the Von Kármán constant (~0.41)
 # d = zero displacement [m], approximated as 2/3 to 3/4 of the average height of the obstacles
 # z0 = roughness length [m] (open terrain: 0.01-0.05 m, bush: 0.1-0.25 m, forest: 0.5-1.0 m, urbain: 1-5m)
 return ustar/κ*log.((h .-d)/z0) 
end

function wind_profile_power(h, v_r, h_r; α=1/7)
    # given wind speed v_r at height h_r, returns a wind profile
    # better than log profile above 2000m
    return v_r*(h/h_r).^α
end


function plot_cn2(h, cn2) # altitude = x axis
    fig=figure("Cn2 profile")
    plot(h, cn2)
    yscale("log")
    xlabel("Altitude (m)")
    ylabel(latexstring("C_n^2~~(m^{-2/3})"))
end


function plot_cn2_alt(h, cn2)# altitude = y axis
    fig=figure("Cn2 profile")
    plot(cn2, h)
    xscale("log")
    ylabel("Altitude (m)")
    xlabel(latexstring("C_n^2~~(m^{-2/3})"))
    ylim([minimum(h), maximum(h)])
end

function cn2_profile_to_r0(cn2, h; ζ=0, λ=500e-9, unit="m")
    # "Top of atmosphere" result
    # Valid for plane wavefront (Schmidt pp 164)
    # Returns r0 in [m]
    # ζ  zenith angle of the observation direction from the observer to the target in degrees
    # HV57 should return r0~5cm
    # h=0:30000; cn2 = cn2_huffnagel_valley(h); cn2_profile_to_r0(cn2, h, unit="cm")
    k = 2pi/λ
    cn2_int = integrate(h, cn2, SimpsonEven())
    r0 = (0.423*k^2*sec(ζ/180*pi)*cn2_int)^(-3/5)
    if unit == "cm"
        r0 *= 100.0
    end
    return r0
end

function cn2_profile_to_θ0(cn2, h; ζ=0, λ=500e-9, unit = "rad")
    # Returns θ0 in [rad]
    # ζ in degrees
    k = 2pi/λ
    cn2_int = integrate(h, cn2.*(h).^(5/3), SimpsonEven())
    θ0 = (2.91*k^2*sec(ζ/180*pi)^(8/3)*cn2_int)^(-3/5)
    if unit == "arcsec"
        θ0 *= 180*3600/pi
    end
    return θ0
end


function cn2_profile_to_r0_range(cn2, h; ζ=0, λ=500e-9, unit="m")
    # Note that CN2 profile here is 
    # Returns r0 in [m] by default
    # ζ zenith angle of the observation direction from the observer to the target in degrees
    k = 2pi/λ
    Δh = maximum(h)
    cn2_int = integrate(h, cn2.*(1.0 .- h/Δh).^(5/3), SimpsonEven())
    r0 = (0.423*k^2*sec(ζ/180*pi)*cn2_int)^(-3/5)
    if unit == "cm"
        r0 *= 100.0
    end
    return r0
end

function cn2_profile_to_θ0_range(cn2, h; ζ=0, λ=500e-9, unit = "rad")
    # Returns θ0 in [rad]
    # ζ in degrees
    k = 2pi/λ
    Δh = maximum(h)
    cn2_int = integrate(h, cn2.*(h/Δh).^(5/3), SimpsonEven())
    θ0 = (2.91*k^2*(Δh)^(5/3)*sec(ζ/180*pi)^(8/3)*cn2_int)^(-3/5)
    if unit == "arcsec"
        θ0 *= 180*3600/pi
    end
    return θ0
end

function mtf_atm_fried(f, f0, r0, α )
# α = 0 : long exposure
# α = 1 : short exposure without scintillation
# α = 1/2 : short exposure with scintillation
return exp.(-3.44*(1.0 .- α*(f/(2f0)).^(1/3) ).*( f/(2f0)*(D/r0) ).^(5/3)    )
end

function mtf_diffraction(ρ, ρc)
    # ρc = D/(λl) where l=focal length of camera
    mtf = zeros(size(ρ))
    indx = findall(abs.(ρ) .<= 2*ρc)
    f = 0.5*abs.(ρ[indx]/ρc);
    mtf[indx] = (2/pi*(acos.(f) - f.*sqrt.(1.0 .-f.^2)))
    return mtf
end
    
    



# function cn2_to_seeing(cn2; λ=500.E-9)
#     #    Calculates the seeing angle from the integrated Cn2 value
#     #   Parameters:
#     #    cn2 (float) integrated Cn2 value in m^2/3
#     #    λ : wavelength
#     # Returns:
#     #    seeing angle in arcseconds
#     r0 = cn2_to_r0(cn2,λ)
#     seeing = r0_to_seeing(r0,λ)
#     return seeing
# end

# function seeing_to_cn2(seeing; λ=500.E-9)
    
#     # Calculates the integrated Cn2 value from the seeing
#     # Parameters:
#     #     seeing (float) seeing in arcseconds
#     #     λ : wavelength
#     # Returns:
#     #     integrated Cn2 value in m^2/3
    
#     r0 = seeing_to_r0(seeing,λ)
#     cn2 = r0_to_cn2(r0,λ)
#     return cn2
# end

# function cn2_to_r0(cn2; λ=500.E-9)
#     # Calculates r0 from the integrated Cn2 value
#     # Parameters:
#     #     cn2 (float) integrated Cn2 value in m^2/3
#     #     λ : wavelength
#     # Returns:
#     #     r0 in m
#     r0=(0.423*(2*pi/λ)^2*cn2)^(-3.0/5.0)
#     return r0
# end

# function r0_to_cn2(r0; λ=500.E-9)
    
#     # Calculates integrated Cn2 value from r0
#     # Parameters:
#     #     r0 (float) r0 in cm
#     #     λ : wavelength
#     # Returns:
#     #     cn2 (float) integrated Cn2 value in m^2/3
    
#     cn2 = r0^(-5.0/3.0)/(0.423*(2*pi/λ)^2)
#     return cn2
# end

# function r0_to_seeing(r0; λ=500.E-9)

#     # Calculates the seeing angle from r0
#     # Parameters:
#     #     r0 (float) Freid's parameter in cm
#     #     λ : wavelength
#     # Returns:
#     #     seeing angle in arcseconds
   
#     return (0.98*λ/r0)*180.0*3600.0/pi
# end

# function seeing_to_r0(seeing; λ=500.E-9)
    
#     # Calculates r0 from seeing
#     # Parameters:
#     #     seeing (float) seeing angle in arcseconds
#     #     λ : wavelength
#     # Returns:
#     #     r0 (float) Fried's parameter in cm
#     return 0.98*λ/(seeing*pi/(180.0*3600.))
# end
