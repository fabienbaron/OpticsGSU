#
# WARNING: Noiseless
#

using OpticsGSU
using FFTW
obj=rotl90(readfits("./data/jupiter.fits")*1.0)
npix=size(obj,1)
# Exercice 1 : create Gaussian PSF with sigma=5 pixels
psf_gaussian=gaussian2d(npix,npix,5)
image = conv_psf_obj(psf_gaussian, obj);
imview(image, title="Image");

# Exercice 2: Read phase screen
phase=readfits("./data/Dr03_phases.fits.fz");  # turbulence
α = 2.72219 # scaling factor for phase
phase *= α
aperture=circular_aperture(npix=npix,diameter=npix/2.,centered=true);
npad = div(max(npix-size(phase,1),0),2); # we may need to pad the phase screen



# Compute speckle data
nframes = size(phase)[3]
image_data = zeros(Float64, npix, npix, nframes )
for i=1:nframes
   psf = pupil_to_psf(aperture, pad(phase[:,:,i],npad));
   image_data[:,:,i] = conv_psf_obj(psf, obj);
end
writefits(image_data, "./data/speckle_fake_data.fits") # One can use Qfitsview to examine this



# Direct inversion on single frame
psf = pupil_to_psf(aperture, pad(phase[:,:,1],npad));


# no noise
image_spec = conv_psf_obj(psf, obj) + 0*randn(npix,npix);
direct_obj=real(fftshift(ifft(fft(image_spec)./(fft(psf).+1.0e-15))))      # adding small number because psf have 0, and we couldn't divided by 0.
imview2(image_spec, direct_obj, caption1="Speckle image", caption2="Direct Inversion", figtitle="Direct Inversion")

# A slight amount of noise
image_spec = conv_psf_obj(psf, obj) + randn(npix,npix)/10;
direct_obj=real(fftshift(ifft(fft(image_spec)./(fft(psf).+1.0e-15))))      # adding small number because psf have 0, and we couldn't divided by 0.
imview2(image_spec, direct_obj, caption1="Speckle image", caption2="Direct Inversion", figtitle="Direct Inversion")


# The Wiener filter is more resilient to noise and uses the entire dataset
wiener_numer = zeros(Complex{Float64},npix,npix);
wiener_denom = zeros(Complex{Float64},npix,npix);
image_long_exposure = zeros(Float64,npix,npix);
# Wiener filter over several frames
for i=1:500
   psf = pupil_to_psf(aperture, pad(phase[:,:,i],npad));
   image_spec = conv_psf_obj(psf, obj) + randn(npix,npix)/10;
   image_long_exposure += image_spec;
   wiener_numer += fft(image_spec).*conj(fft(psf))
   wiener_denom += abs2.(fft(psf)) .+ 1.0e-15
end
wiener_obj =real(fftshift(ifft(wiener_numer./wiener_denom)))
imview2(image_long_exposure, wiener_obj, caption1="Long exposure image", caption2="Wiener Inversion", figtitle="Wiener Inversion")
