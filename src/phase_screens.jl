using FFTW, Random 

mutable struct Atmosphere
    N::Int64 # phase/amp screen output size
    λ::Array{Float64,1}
    heights::Array{Float64,1}
    l0::Array{Float64,1}
    L0::Array{Float64,1}
    Cn2::Array{Float64,1}
    wind_velocities::Matrix{Float64}
    # Constructor should be able to fill the next ones:
    nλ::Int64
    nlayers::Int64
    r0::Array{Float64, 2}
    phase_screens::Array{Float64, 4}
    composite_amplitude::Array{Float64, 3}
    composite_phase::Array{Float64, 3}
    function Atmosphere(N::Int64, λ::Array{Float64,1}, heights::Array{Float64,1}, l0::Array{Float64,1}, L0::Array{Float64,1}, Cn2::Array{Float64,1}, wind_velocities::Matrix{Float64})
        nλ = length(λ)
        nlayers = length(heights)
        r0 = zeros(Float64, nlayers, nλ);
        phase_screens    = zeros(Float64, 2*N, 2*N, nlayers, nλ);
        composite_amplitude  = zeros(Float64, N, N, nλ);
        composite_phase      = zeros(Float64, N, N, nλ);
        new( N, λ, heights, l0, L0,Cn2,wind_velocities,nλ,nlayers, r0, phase_screens, composite_amplitude,composite_phase);
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


function ang_spec_multi_prop(Uin, t, λ, δ1, δn, z)
	# Following Schmidt Listing 9.1, page 151
	# z = heights of atmosphere phase screens, increasing in altitude, must include z=0
	# t = complex transmission [:,:,nlayers]
	# δ1 = grid sampling for first layer 
	# δn = grid sampling for last layer
	# Function for simulating scintillation through the atmosphere

# TEST Uin = U; t=sg.*cis.(phz[:,:,:,l]); λ = lam[l];  δ1=delta1; δn = deltan; z = dz

    N = size(Uin, 1); # Number of side grid points
    ρ = meshrad([-N/2:N/2-1;])
	ρ2 = ρ.^2
    k = 2*π/λ;    # Optical wavevector
    sg = exp.(-ρ2.^8/(0.47*N)^16);  # Super Gaussian absorbing boundary
	n = length(z);
    # Propagation distances
    Δz = z[2:n] - z[1:n-1];    # Calculates distances between layers
    # Grid spacings
    α = z / z[n];
    δ = (1 .- α)*δ1 + α*δn;
    m = δ[2:n] ./ δ[1:n-1];
    Q1 = cis.(k/2*(1-m[1])/Δz[1]*ρ2*δ[1]^2);
    Uin = Uin.*Q1.*t[:,:,1];
    for idx in 1:n-1
        δf = 1/(N*δ[idx]); # Spatial frequencies of i^th plane
        Q2 = cis.(-2*π^2*Δz[idx]/m[idx]/k*δf^2*ρ2);   # Quadratic Phase Factor
        Uin = sg .* t[:,:,idx+1].*ift2(Q2.*ft2(Uin/m[idx], δ[idx]), δf); # Compute the propagated field
    end
    # Observation-plane coordinates
    Q3 = cis.(k/2*(m[n-1]-1)/(m[n-1]*Δz[n-1])*(ρ2*δ[n]^2));   # Δ_z[n-1] is used instead of Z. Functionality should be the same.
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

function ft_phase_screen(r0, N, δ, L0, l0, seed)
     # Function for creating a random draw phase screen
	 # Creates a phase screen based on the FT method
	 # Setup the PSD
	 del_f = 1/(N*δ);		    # Frequency grid space [1/m]
	 # Fequency grid [1/m]
	 fx, fy = meshgrid([-N/2:N/2-1;]*del_f);
	 f = sqrt.(fx.^2 + fy.^2);
	 fm = 5.92/l0/(2*π);		# Inner scale frequency		# Use of 'π'
	 f0 = 1/L0;		       	    # Outer scaling frequency
	 # Modified von Karman atmospheric phase PSD
	 PSD_phi = 0.0229*r0^(-5/3)*exp.(-(f/fm).^2).*(f.^2 .+ f0^2).^(-11/6);
	 PSD_phi[div(N,2)+1, div(N,2)+1] = 0;
	 # Random draws of Fourier coefficients
	 Random.seed!(seed)
	 cn = (randn(Float64, N,N)+im*randn(Float64, N,N)).*sqrt.(PSD_phi)*del_f;
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


function instantiate_polychromatic_atmosphere(atmosphere::Atmosphere, D, z) # D will determine sampling
    nλ=atmosphere.nλ;
    nlayers=atmosphere.nlayers
    k = (2*pi)./atmosphere.λ;
    N = 2*atmosphere.N        # phase screen internal generation size
    N_final = atmosphere.N    # phase/amp screen output size
    #########grid spacings#########
    delta1= 4*D/N; #3.52e-3; #source-plane grid spacing [m]
    deltan= 4*D/N; #3.52e-3; #observation-plane grid spacing [m]
    alpha = atmosphere.heights / atmosphere.heights[nlayers];
    delta = (1.0 .-alpha)*delta1+(alpha*deltan);
    del_z = atmosphere.heights[nlayers]-atmosphere.heights[nlayers-1];
    #R = Dz; # wavefront radius of curvature [m]
    ρ = meshrad([-N/2:N/2-1;] * delta1);
    sg = exp.(-(ρ/(0.47*N)).^16);  # Super Gaussian absorbing boundary
    seeds= shuffle(1:nlayers);
    Threads.@threads for l = 1:nλ
        println("Propagating $l/$nλ")
        atmosphere.r0[:,l] = (0.423 * (k[l].^2) * (cos(z)^-1) * atmosphere.Cn2 * del_z).^(-3/5); 
        for ilayer = 1:nlayers
           atmosphere.phase_screens[:,:,ilayer,l] = ft_phase_screen(atmosphere.r0[ilayer,l], N, delta[ilayer], atmosphere.L0[ilayer], atmosphere.l0[ilayer],seeds[ilayer]);
        end
        U = ones(N, N); # Amplitude at the top of the atmosphere
        U = ang_spec_multi_prop(U, sg.*cis.(atmosphere.phase_screens[:,:,:,l]), atmosphere.λ[l], delta1, deltan, atmosphere.heights[2:end]);
        lo = div(N-N_final,2)+1
        hi = div(N+N_final,2)
        atmosphere.composite_amplitude[:,:,l] = abs.(U[lo:hi,lo:hi]);
        atmosphere.composite_phase[:,:,l]     = angle.(U[lo:hi,lo:hi]);
    end
    return atmosphere
end


function phase_screen_size(D, pixel_size, wind, Δt)
    # D: diameter of the telescope
    # pixel_size: size of pixels on telescope pupil

    # Diameter of phase screen should be larger than diameter of telescope aperture
    # (ideally close to the outer scale size for the atmosphere)
    Δpos_meters = [Δt.*wind[:, 1].*sin.(wind[:, 2].*pi/180) Δt.*wind[:, 1].*cos.(wind[:, 2].*pi/180)]'
    Δpos_pix = Δpos_meters ./ pixel_size
    Δpos_pix_total = Δpos_pix .* nepochs
    npix_screen = Int(2^ceil(log(2, maximum(ceil.(abs.(Δpos_pix_total)))+pdim)))

    D1=2.0*D
    npts=round(Int, D1/pixel_size) 
    # Compute number of points to be compatible with 2^N points for the FFT
    dim=ceil(Int, 2^ceil(log10(npts)/log10(2.0)+1));
    println("Number of pixels along one side of phase screen = ",dim);
    Dscreen=dim*D1/float(npts);
    println("Physical size of phase screen = ",Dscreen, " meters"); 
    println("Physical size of pixels = ",Dscreen/dim, " meters"); 
    return dim
 end 
    