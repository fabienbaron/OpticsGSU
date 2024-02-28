using LinearAlgebra, Statistics, JLD2, Random, OpticsGSU
#Random.seed!(1667)
# Diameter of telescope aperture [m]#
D=3.6 
pixscale_wanted = 4e-2 # [m/pix] pupil min sampling
N = pupil_support_size(D, pixscale_wanted) # pick size ensuring minimum sampling
pixscale = D/(N/2) 

# min and max observing wavelengths [m]
λmin=400.0e-9
λmax=500.0e-9
resolution = 0.25*λmin/D*1e6/2.0 # assuming Nyquist sampling of pupil at 400 nm
println("Nyquist pixel scale = ", resolution," [arcsec]")

#
# DETECTOR SETUP
#
nframes = 1000 # Number of frames
exptime = 5e-3; # exposure time for each frame
timestamps = (0:nframes-1)*exptime 
detector = Detector(true, true, UInt16, 1.0, 1.0, 2^16-1, 2.0, exptime)
#poisson::Bool, adu::Bool, aduTYPE::DataType, qe::Array{Float32,1}, gain::Float32, saturation::Int32,σ_ron::Float32, exptime::Float32 

#
# APERTURE SETUP
#
aperture_mask, λ = generate_aperture_chromatic_steps(N, λmin, λmax, delta_slice=2)
nwavs = length(λ);

# reverse so λ increasing?

#
# OBJECT SETUP
#

#object, abundances, spectra = generate_hyperspectral_object(N, λ, template= "./data/sat_template2.fits");
object, ~, ~ = generate_hyperspectral_object(N, λ, template= "./data/sat_template2.fits");
z = 17/360.0*2*pi;            # observation: angular distance from zenith [radians]
mag1 = 4
f1 = flux(mag1, tel_surf=pi*(D/2)^2, airmass = 1.0, exptime=detector.exptime);
object *= f1 / sum(object)
FTobject = object_to_fto(object);
object_distance = 35786e3 # [m] - closest distance to object
object_sampling = copy(pixscale)

#
# ATMOSPHERE 
#
elevation = 2400; # observation: elevation
nlayers = 3; # number of atmospheric layers
winds = Float32.([  4.6 45.0 ; 5.0 -30; 3.1 25 ])     # (m/s, deg) 0deg = East, increases clockwise
l0 = collect(range(3e-3,3e-2,length=nlayers));
L0 = collect(range(10,2000,length=nlayers));
layer_heights = elevation .+ [0, 10000, 30000, 40000]; 
Cn2 = h->CN2_huffnagel_valley_generalized(Float64.(h), A=9.9e-15, C = 5.94e-53);
Nscreens = phase_screen_size(D, timestamps, winds, pixscale); # atmosphere screen size so that pixels are ~2cm
atmosphere = Atmosphere(Nscreens, λ, Float32.(layer_heights), Float32.(l0), Float32.(L0), Cn2, Float32.(winds)); # This will construct the atmosphere with blank screens
atmosphere = instantiate_atmosphere(atmosphere, D, z, verbose = true, propagate=false);
I = generate_isoplanatic_frozen_flow_phase_extractors(atmosphere, timestamps, N, object_distance, pixscale);

#Note: by default generate_isoplanatic_data creates broadband data
data, psfs, otfs, pupil_amps, pupil_phases = generate_isoplanatic_data(timestamps, λ, I, FTobject, aperture_mask, atmosphere, detector, disperse = false, λref_disp = 500e-9, zenith = z, pixscale = pixscale);

savefile= "speckle_sat_n"*string(N)*"_nframes"*string(nframes)*"_nλ"*string(length(λ))*"_mag"*string(mag1)*"_great.jld2"
println("Saving file as: $savefile")
jldsave(savefile; λ,  timestamps, I, aperture_mask, data, psfs, otfs, pupil_amps, pupil_phases, detector, object, FTobject, atmosphere)

writefits(psf_broad, "psfs.fits")
writefits(data, "data.fits")



data=readfits("data.fits")

shift_img = shift_and_add(Float64.(data))
imview(shift_img)






