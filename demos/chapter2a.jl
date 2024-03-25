using OpticsGSU
using FFTW
obj=readfits("./data/jupiter.fits")*1.0
npix=size(obj,1)
# Exercice 1 : create Gaussian PSF with sigma=5 pixels
psf_gaussian=gaussian2d(npix,npix,5)
image = convolve(psf_gaussian, obj);
image = image + maximum(image)/50*randn(Float64, size(image))
imview(image, title="Image");

# Exercice 2: Read phase screen
phase=readfits("./data/Dr03_phases.fits.fz");  # turbulence
aperture=circular_aperture(npix=npix,diameter=npix/2.,centered=true);
npad = div(max(npix-size(phase,1),0),2); # we may need to pad the phase screen

psf = pupil_to_psf(aperture, pad(phase[:,:,1]*5,npad));
image = convolve(psf, obj);
σ=maximum(obj)/50
image = image + σ*randn(Float64, size(image))

# Direct inversion on single frame
direct_obj=real(fftshift(ifft(fft(image)./(fft(psf).+1.0e-10))))      # adding small number because psf have 0, and we couldn't divided by 0.
imview2(image, direct_obj, caption1="Detector image", caption2="Direct Inversion", figtitle="Direct Inversion")

# Wiener filter over several frames
nframes=500
wiener_numer = zeros(Complex{Float64},npix,npix);
wiener_denom = zeros(Complex{Float64},npix,npix);
images = zeros(Float64, npix, npix, nframes);
psfs = zeros(Float64, npix, npix, nframes);
σ=200000
μ=1e6
for i=1:nframes
   psfs[:,:,i] = pupil_to_psf(aperture, pad(phase[:,:,i]*5,npad));
   images[:,:,i] = convolve(psf, obj) + σ*randn(Float64, size(obj));
   otf_i = ft2(psfs[:,:,i])
   wiener_numer += ft2(images[:,:,i]).*conj(otf_i)
   wiener_denom += abs2.(otf_i) .+ μ
end
wiener_obj =real(ift2(wiener_numer./wiener_denom))
imview2(images[:,:,123], wiener_obj, caption1="Detector image", caption2="Wiener Inversion", figtitle="Wiener Inversion")
