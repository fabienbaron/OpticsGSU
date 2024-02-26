using FFTW, Statistics, Random, ProgressMeter, NumericalIntegration


mutable struct Atmosphere
    N::Int64 # phase/amp screen output size
    λ::Array{Float32,1}
    heights::Array{Float32,1}
    l0::Array{Float32,1}
    L0::Array{Float32,1}
    Cn2
    winds::Matrix{Float32}
    # Constructor should be able to fill the next ones:
    nλ::Int64
    nlayers::Int64
    r0::Array{Float32, 2}
    phase_screens::Array{Float32, 4}
    composite_amplitude::Array{Float32, 3}
    composite_phase::Array{Float32, 3}
    function Atmosphere(N::Int64, λ::Array{Float32,1}, heights::Array{Float32,1}, l0::Array{Float32,1}, L0::Array{Float32,1}, Cn2, wind_velocities::Matrix{Float32})
        nλ = length(λ)
        nlayers = length(heights)-1
        r0 = zeros(Float32, nlayers, nλ);
        phase_screens    = zeros(Float32, N, N, nlayers, nλ); # We will need 4x larger if we propagate
        composite_amplitude  = ones(Float32, N, N, nλ);
        composite_phase      = zeros(Float32, N, N, nλ);
        new(N, λ, heights, l0, L0,Cn2,wind_velocities,nλ,nlayers, r0, phase_screens, composite_amplitude,composite_phase);
    end
end

function str_fcn2_ft(ph, mask=1.0, δ=1.0)
    # Function for calculating the structure function of a randomly distributed field
    N  = size(ph, 1);
    ph = ph.*mask;
    P = ft2(ph, δ);
    S = ft2(ph.^2, δ);
    W = ft2(mask, δ);
    δ_f = 1/(N*δ);
    w2 = ift2(W.*conj(W),δ_f);
    D = 2*ift2(real(S.*conj(W)) - abs.(P).^2, δ_f)./ w2 .* mask;
    return real(D)
end


function phase_to_dr0(phase; mask=[])
    if mask == []
        dr0 = (var(phase)/1.0299)^(3/5) # original D/r0 of the *square* phase screen
        return dr0;
    else
        # Compute D/r0 over the relevant pixels (= within the aperture)
        indx= findall(mask.>1e-8);
        var_pupil = var(phase[indx]);
        dr0 = (var_pupil/1.0299)^(3/5);
        return dr0
    end
end



function Dr0_to_Cn2(Dr0, D, λ, Dz)
    k = 2pi/minimum(λ)
    Cn2 = (D/Dr0)^(-5/3) / (0.423 * (k.^2) * (cos(z)^-1)*Dz);
    return Cn2
end

@views function ang_spec_multi_prop(Uin, t, ρ2, sg, λ0, δ1, δn, heights; FTYPE=Float32)
	# Following Schmidt Listing 9.1, page 151
	# heights = heights of atmosphere phase screens, increasing in altitude, must include heights=0
	# t = complex transmission [:,:,nlayers]
	# δ1 = grid sampling for first layer 
	# δn = grid sampling for last layer
	# Function for simulating scintillation through the atmosphere
    # TEST Uin = U; t=cis.(atmosphere.phase_screens[:,:,:,l]); λ0 = atmosphere.λ[l]; heights = atmosphere.heights; ρ2=ρ.^2
    N = size(Uin, 1); # Number of side grid points
    n = length(heights);
    k = 2*π/λ0;    # Optical wavevector
    # Propagation distances
    Δheights = heights[2:n] - heights[1:n-1];    # Calculates distances between layers
    # Grid spacings
    α = heights / heights[n];
    δ = (FTYPE(1.0) .- α)*δ1 + α*δn;
    m = δ[2:n] ./ δ[1:n-1];
    Q1 = cis.(FTYPE(k/2)*(1-m[1])/Δheights[1]*ρ2*δ[1]^2);
    Uin = Uin.*Q1.*sg.*t[:,:,1];
    for idx in 1:n-1
        δf = 1/(N*δ[idx]); # Spatial frequencies of i^th plane
        Q2 = cis.(FTYPE(-2*π^2/k)*Δheights[idx]/m[idx]*δf^2*ρ2);   # Quadratic Phase Factor
        Uin = sg .* t[:,:,idx+1].*ift2(Q2.*ft2(Uin/m[idx], δ[idx]), δf); # Compute the propagated field
    end
    # Observation-plane coordinates
    Q3 = cis.(FTYPE(k/2)*(m[n-1]-1)/(m[n-1]*Δheights[n-1])*(ρ2*δ[n]^2));   # Δ_heights[n-1] is used instead of heights. Functionality should be the same.
    Uout = Q3 .* Uin;
    return Uout
end

function ft_scint_screen(r0, N, δ, L0, l0, h=5000.0, λ=0.5e-6)
     # Function for creating a random draw scintilation screen
	 # Creates a scint screen based on the FT method
     # Height of layer, assumed to be a high level [m] unless otherwise stated
	 # Setup the PSD
	 # Sedmak 2004  Implementation of fast-Fourier-transform-based
	 # simulations of extra-large atmospheric phase and scintillation screens
	 del_f = 1/(N*δ);		    # Frequency grid space [1/m]
	 fx = [-N/2:N/2-1;] * del_f;
	 # Frequency grid [1/m]
	 fxx = [j for i=fx, j=fx];
	 # Translate from cartesian to polar
	 #th = atan(fy/fx);
	 f = sqrt.(fxx.^2 + fxx'.^2);
	 fm = 5.92/l0/(2*π);		    # Inner scale frequency		# Use of 'π'
	 f0 = 1/L0;		       	    # Outer scaling frequency
     Cn2 = 7.53e12   # Calculated using values in Sedmak 2004 using h = 5000m
	 # Modified von Karman atmospheric scint PSD from Sedmak 2004 - Eqn. 10 (unfiltered).
     # Assuming single layer at h
     # This version is different to that in the course notes
     PSD_sci = 1.54 * λ^(-2) * f.^(-11/3) * Cn2 * (sin.(π*λ*h*f.^2)).^2;
     PSD_sci[div(N,2)+1, div(N,2)+1] = 0;
	 # Random draws of Fourier coefficients
	 cn = (randn(Float64, N,N)+im*randn(Float64, N,N)).*sqrt.(PSD_sci)*del_f;
	 # Synthesize the phase screen
	 scint = real(ift2(cn, 1.0)); # FB: why 1.0?
	 return scint
end

function ft_phase_screen(r0, N, δ, L0, l0, seed; FTYPE=Float32)
     # Function for creating a random draw phase screen
	 # Creates a phase screen based on the FT method
	 # Setup the PSD
	 del_f = 1/(N*δ);		    # Frequency grid space [1/m]
	 # Frequency grid [1/m]
	 fx, fy = meshgrid([-N/2:N/2-1;]*del_f);
	 f = sqrt.(fx.^2 + fy.^2);
	 fm = 5.92/l0/(2*π);		# Inner scale frequency		# Use of 'π'
	 f0 = 1/L0;		       	    # Outer scaling frequency
	 # Modified von Karman atmospheric phase PSD
	 PSD_phi = 0.0229*r0^(-5/3)*exp.(-(f/fm).^2).*(f.^2 .+ f0^2).^(-11/6);
	 PSD_phi[div(N,2)+1, div(N,2)+1] = 0;
	 # Random draws of Fourier coefficients
	 Random.seed!(seed)
	 cn = (randn(FTYPE, N,N)+im*randn(FTYPE, N,N)).*sqrt.(PSD_phi)*del_f;
	 # Synthesize the phase screen
	 return real(ift2(cn, 1));
end


function ft_sh_phase_screen(r0, N, δ, L0, l0)
	 # Following Schmidt Listing 9.3, page 170
	# Augments the above phase screen method to include subharmonics to improve the low spatial frquencies.
	 D = N*δ;
	 phz_hi = ft_phase_screen(r0, N, δ, L0, l0); 	 # High frequency screen from FFT method
	 x, y = meshgrid([-N/2:N/2-1;]*δ);	 # Spatial grid [m]
	 phz_lo = zeros(size(phz_hi)); 	 # Initialise low-freq screen,
	 # loop over frequency grids with spacing 1/(3^p*L)
	 for p = 1:3
	     del_f = 1 / (3^p*D);		# Frequency grid spacing [1/m]
		 fx, fy = meshgrid([-1.0:1.0;]*del_f)
		 f  = sqrt.(fx.^2+fy.^2)
	     fm = 5.92/l0/(2*π);		# Inner scale frequency [1/m]  # 'π' may not work - Matt
	     f0 = 1/L0; 			    # Outer scale frequency [1/m]
	     PSD_phi = 0.023*r0^(-5/3)* exp.(-(f/fm).^2)./ (f.^2 .+ f0^2).^(11/6); 	     # modified von Karmen atmospheric phase PSD
	     PSD_phi[2,2] = 0.0;
	     cn = (randn(3) + im*randn(3)).* sqrt.(PSD_phi)*del_f; 	     # Random draws of Fourier coefficients
	     SH = zeros(Float64, N,N);
	     # Loop over frequencies on this grid
	     for ii = 1:9
	     	 SH += real(cn[ii] * cis.(2*π*(fx[ii]*x+fy[ii]*y))); 
	     end
	 phz_lo += SH;			# Accumulate subharmonics
     end
	 phz_lo .-= mean(real(phz_lo));
	 return phz_lo, phz_hi
end


@views function instantiate_atmosphere(atmosphere::Atmosphere, D, z; verbose=false, propagate = false, r0avg = -1, FTYPE=Float32) # D will determine sampling
    #FTYPE=Float32; verbose = true; propagate = true;
    nλ=atmosphere.nλ;
    nlayers=atmosphere.nlayers
    k = (2*pi)./atmosphere.λ;
    seeds= shuffle(1:nlayers);
    if propagate == true 
        Nscreen_double = 2*atmosphere.N        # phase screen internal generation size
        Nscreen = atmosphere.N    # phase/amp screen output size
        atmosphere.phase_screens    = Array{FTYPE}(undef, Nscreen_double, Nscreen_double, nlayers, nλ); # We will need 4x larger if we propagate
        atmosphere.composite_amplitude  = Array{FTYPE}(undef, Nscreen, Nscreen, nλ);
        atmosphere.composite_phase      = Array{FTYPE}(undef, Nscreen, Nscreen, nλ);
        lo = div(Nscreen_double-Nscreen,2)+1
        hi = div(Nscreen_double+Nscreen,2)
        #########grid spacings#########
        del_z = atmosphere.heights[nlayers]-atmosphere.heights[nlayers-1]; # Fix this for nlayers = 1
        δ1= FTYPE(4*D/Nscreen_double); #3.52e-3; #source-plane grid spacing [m]
        δn= FTYPE(4*D/Nscreen_double); #3.52e-3; #observation-plane grid spacing [m]
        alpha = atmosphere.heights / atmosphere.heights[nlayers];
        δ = ( (FTYPE(1.0) .-alpha)*δ1+(alpha*δn) );
        # Big matrices for propagation
        ρ = meshrad([-Nscreen:Nscreen-1;] * δ1);
        sg = exp.(-(ρ/(FTYPE(0.47)*Nscreen)).^16);  # Super Gaussian absorbing boundary
        Utop = ones(FTYPE, Nscreen_double, Nscreen_double); # Amplitude at the top of the atmosphere
        for ilayer = 1:nlayers
            atmosphere.r0[ilayer,1] = cn2_profile_to_r0(atmosphere.Cn2(atmosphere.heights[ilayer]:atmosphere.heights[ilayer+1]), atmosphere.heights[ilayer]:atmosphere.heights[ilayer+1], ζ=z, λ=atmosphere.λ[1]) 
        end
        p = Progress(nλ, desc="Generating & propagating atmospheric layers")
        Threads.@threads for l = 1:nλ
            for ilayer = 1:nlayers
                atmosphere.r0[ilayer,l] =  atmosphere.r0[ilayer,1]*(atmosphere.λ[l]/atmosphere.λ[1])^(6/5)
                atmosphere.phase_screens[:,:,ilayer,l] = ft_phase_screen(atmosphere.r0[ilayer,l], Nscreen_double, δ[ilayer], atmosphere.L0[ilayer], atmosphere.l0[ilayer],seeds[ilayer]);
            end            
            U = ang_spec_multi_prop(Utop, cis.(@views atmosphere.phase_screens[:,:,:,l]), ρ.^2, sg, atmosphere.λ[l], δ1, δn, atmosphere.heights[2:end]);
            atmosphere.composite_amplitude[:,:,l] = abs.(U[lo:hi,lo:hi]);
            atmosphere.composite_phase[:,:,l]     = angle.(U[lo:hi,lo:hi]);
            next!(p)
            end
            finish!(p)
    else # No propagation
        Nscreen = atmosphere.N        # phase screen internal generation size
        atmosphere.phase_screens    = Array{FTYPE}(undef, Nscreen, Nscreen, nlayers, nλ); # We will need 4x larger if we propagate
        atmosphere.composite_amplitude  = ones(FTYPE, Nscreen, Nscreen, nλ);
        atmosphere.composite_phase      = Array{FTYPE}(undef, Nscreen, Nscreen, nλ);
        δ1= FTYPE(4*D/Nscreen); #3.52e-3; #source-plane grid spacing [m]
        δn= FTYPE(4*D/Nscreen); #3.52e-3; #observation-plane grid spacing [m]
        alpha = atmosphere.heights / atmosphere.heights[nlayers];
        δ = ( (FTYPE(1.0) .-alpha)*δ1+(alpha*δn) );
        p = Progress(nλ, desc="Generating & compositing atmospheric layers")
        for ilayer = 1:nlayers
            atmosphere.r0[ilayer,1] = cn2_profile_to_r0(atmosphere.Cn2(atmosphere.heights[ilayer]:atmosphere.heights[ilayer+1]), atmosphere.heights[ilayer]:atmosphere.heights[ilayer+1], ζ=z, λ=atmosphere.λ[1]) 
        end
        Threads.@threads for l = 1:nλ
            for ilayer = 1:nlayers
                atmosphere.r0[ilayer,l] =  atmosphere.r0[ilayer,1]*(atmosphere.λ[l]/atmosphere.λ[1])^(6/5)
                atmosphere.phase_screens[:,:,ilayer,l] = ft_phase_screen(atmosphere.r0[ilayer,l], Nscreen, δ[ilayer], atmosphere.L0[ilayer], atmosphere.l0[ilayer],seeds[ilayer]);
            end
            atmosphere.composite_phase[:,:,l] = dropdims(sum(atmosphere.phase_screens[:,:,:,l], dims=3),dims=3);
            next!(p) 
        end
        finish!(p)
    end
    return atmosphere
end


function phase_screen_size(D, t, winds, pixscale) # Computed the support size of the phase screens
    # D: diameter of the telescope
    # pixel_size: size of pixels on telescope pupil
    # Diameter of phase screen should be larger than diameter of telescope aperture
    # (ideally close to the outer scale size for the atmosphere)
    # Compute number of points to be compatible with 2^N points for the FFT
    Δt = maximum(t)-minimum(t)
    println("Maximum elapsed time = ", Δt, " seconds"); 
    Δpos= maximum((abs.([Δt*winds[:, 1].*sin.(winds[:, 2].*pi/180) Δt*winds[:, 1].*cos.(winds[:, 2].*pi/180)])))
    Δpos_pix = ceil(Int, Δpos/pixscale)
    println("Maximum length of swept phase = ", Δpos, " m, corresponding to ", Δpos_pix ," pixels");
    npix_pupil = 2*ceil(Int, D/pixscale) 
    npix_screen = nextprod((2, 3, 5, 7), Δpos_pix+npix_pupil)
    println("Number of pixels along one side of phase screen = ", npix_screen)
    return npix_screen
 end 
 
function get_non_overlapping_punch_outs(t, npix_pupil, atmosphere)
    nframes = length(t);
    # Maximimize the use of the phase screen
    npix_screens = size(atmosphere.composite_phase,1)
    nside = floor(Int, (npix_screens-1)/npix_pupil); 
    if (nside+1)^2<nframes
        error("Phase screen size is not large enough to get $nframes frames as $(nside) x $(nside) punch outs")
    end
    # Upper left corner of punch out
    x = (mod.(0:nframes-1, nside))*npix_pupil .+1;
    y = (div.(0:nframes-1, nside))*npix_pupil .+1;
    # toto = zeros(npix_screens,npix_screens)
    # for i=1:nframes 
    #     toto[x[i]:x[i]+npix_pupil-1, y[i]:y[i]+npix_pupil-1] .+=1.0
    # end
    # # using PyCall; patches = pyimport("matplotlib.patches")
    # for i=1:nframes
    #     rect = patches.Rectangle((x[i], y[i]), npix_pupil-1, npix_pupil-1, linewidth=1, edgecolor="r", facecolor="none")
    # end
    #imshow(atmosphere.composite_phase[:,:,1])
    #scatter(x,y)
    # We return crop operators that will apply to the phase
    M = [(z,k)->z[x[n]:x[n]+npix_pupil-1,y[n]:y[n]+npix_pupil-1,k] for n=1:nframes]; # punch out function
    return M
end



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
    # h=0:30000; cn2 = CN2_huffnagel_valley(h); cn2_profile_to_r0(cn2, h, unit="cm")
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
    θ0 = (2.914*k^2*sec(ζ/180*pi)^(8/3)*cn2_int)^(-3/5)
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
