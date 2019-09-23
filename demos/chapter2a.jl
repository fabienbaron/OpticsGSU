using OpticsGSU
using FFTW
obj=readfits("./data/jupiter.fits")*1.0
npix=size(obj,1)
# Exercice 1 : create Gaussian PSF with sigma=5 pixels
psf_gaussian=gaussian2d(npix,npix,5)
image = conv_psf_obj(psf_gaussian, obj);
image = image + maximum(image)/50*randn(Float64, size(image))
imview(image, title="Image");

# Exercice 2: Read phase screen
phase=readfits("./data/Dr03_phases.fits.fz");  # turbulence
aperture=circular_aperture(npix=npix,diameter=npix/2.,centered=true);
npad = div(max(npix-size(phase,1),0),2); # we may need to pad the phase screen


# Direct inversion on single frame
psf = pupil_to_psf(aperture, pad(phase[:,:,1]*2.72219,npad));
image_spec = conv_psf_obj(psf, obj);
image_spec = image_spec + maximum(image_spec)/50*randn(Float64, size(image))
direct_obj=real(fftshift(ifft(fft(image_spec)./(fft(psf).+1.0e-15))))      # adding small number because psf have 0, and we couldn't divided by 0.
imview2(image_spec, direct_obj, caption1="Detector image", caption2="Direct Inversion", figtitle="Direct Inversion")

# Wiener filter over several frames
wiener_numer = zeros(Complex{Float64},npix,npix);
wiener_denom = zeros(Complex{Float64},npix,npix);
image_blurred = zeros(Float64, size(image_spec));
for i=1:500
   psf = pupil_to_psf(aperture, pad(phase[:,:,i]*2.72219,npad));
   image_spec_instant = conv_psf_obj(psf, obj);
   image_spec_instant = image_spec_instant + maximum(image_spec_instant)/50*randn(Float64, size(image));
   global image_blurred += image_spec_instant;
   global wiener_numer += fft(image_spec).*conj(fft(psf))
   global wiener_denom += abs2.(fft(psf)) .+ 1.0e-5
end
wiener_obj =real(fftshift(ifft(wiener_numer./wiener_denom)))
imview2(image_blurred, wiener_obj, caption1="Detector image", caption2="Wiener Inversion", figtitle="Wiener Inversion")
