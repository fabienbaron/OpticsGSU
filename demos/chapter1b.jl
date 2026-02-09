#
# Simple speckle deconvolution via Wiener filter
#

using OpticsGSU, FFTW, LinearAlgebra, Noise
obj=rotl90(readfits("./data/jupiter.fits")*1.0)
npix=size(obj,1)
# Exercice 1 : create Gaussian PSF with sigma=5 pixels
psf_gaussian=gaussian2d(npix,npix,5);
image = convolve(psf_gaussian, obj);
imview(image, title="Image");

# Read phase screen
phases=readfits("./data/Dr03_phases.fits.fz");  # turbulence
α = 2.72219 # scaling factor for phase
phases *= α
aperture=circular_aperture(npix=npix,diameter=npix/2.,centered=true);
npad = div(max(npix-size(phases,1),0),2); # we may need to pad the phase screen


# Compute speckle data  - Gaussian Noise
nframes = size(phases)[3]
image_data_noiseless = zeros(Float64, npix, npix, nframes )
image_data_noisy = zeros(Float64, npix, npix, nframes )
psfs = zeros(Float64, npix, npix, nframes )
otfs = zeros(ComplexF64, npix, npix, nframes )
σ = 1/10
for i=1:nframes
   psfs[:,:,i] = pupil_to_psf(aperture, pad(phases[:,:,i],npad));
   otfs[:,:,i] = ft2(psfs[:,:,i]);
   image_data_noiseless[:,:,i] = convolve(psfs[:,:,i], obj);
   image_data_noisy[:,:,i] = image_data_noiseless[:,:,i] + σ*randn(npix,npix);
end
writefits(image_data_noiseless, "./data/speckle_fake_data_noiseless.fits")
writefits(image_data_noisy, "./data/speckle_fake_data_noisy.fits") # One can use Qfitsview to examine this
writefits(psfs, "./data/speckle_fake_data_psfs.fits")
#writefits(otfs, "./data/speckle_fake_data_otfs.fits")


# Direct inversion on single frame
psf = pupil_to_psf(aperture, pad(phases[:,:,1],npad));

# no noise
image_spec = convolve(psf, obj) + 0*randn(npix,npix);
direct_obj=real(fftshift(ifft(fft(image_spec)./(fft(psf).+1.0e-5))))      # adding small number because psf have 0, and we couldn't divided by 0.
imview2(image_spec, direct_obj, caption1="Speckle image", caption2="Direct Inversion (no noise on data)", figtitle="Direct Inversion")

# A slight amount of noise
image_spec = convolve(psf, obj) + randn(npix,npix)/1000;
direct_obj=real(fftshift(ifft(fft(image_spec)./(fft(psf).+1.0e-5))))      # adding small number because psf have 0, and we couldn't divided by 0.
imview2(image_spec, direct_obj, caption1="Speckle image", caption2="Direct Inversion", figtitle="Direct Inversion")


# The Wiener filter is more resilient to noise and uses the entire dataset
# Note that it is forward deconvolution (otf are known), setting the gradient of ∑_t || H_t x - y||^2  to be = 0
# Note σ does not appear since it simplifies out as long as the same noise is used
# Of course this won't work for Poisson noise
image_long_exposure = dropdims(sum(image_data_noisy, dims=3), dims=3);

wiener_numer = zeros(Complex{Float64},npix,npix);
wiener_denom = zeros(Float64,npix,npix);
for i=1:500
   wiener_numer += ft2(image_data_noisy[:,:,i]).*conj(otfs[:,:,i]);
   wiener_denom += abs2.(otfs[:,:,i]) .+ 1e-15;
end
wiener_obj =real(ift2(wiener_numer./wiener_denom));
imview2(image_long_exposure, wiener_obj, caption1="Long exposure image", caption2="Wiener Inversion", figtitle="Wiener Inversion")
norm(wiener_obj-obj)
