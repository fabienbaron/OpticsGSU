using OpticsGSU
obj=read(FITS("jupiter.fits")[1])*1.0
npix=size(obj,1)
# Exercice 1 : create Gaussian PSF with sigma=5 pixels
psf_gaussian=gaussian2d(npix,npix,5)
image = conv_psf_obj(psf_gaussian, obj);

# Exercice 2: Read phase screen
phase=read((FITS("Dr03_phases.fits")[1]));  # turbulence
aperture=circular_aperture(npix=npix,diameter=npix/2,centered=true);

wiener_numer = zeros(Complex64,npix,npix);
wiener_denom = zeros(Complex64,npix,npix);
npad = div(max(npix-size(phase,1),0),2); # we may need to pad the phase screen


# Direct inversion on single frame
psf = pupil_to_psf(aperture, pad(phase[:,:,1]*2.72219,npad));
image_spec = conv_psf_obj(psf, obj);
direct_obj=real(fftshift(ifft(fft(image_spec)./(fft(psf).+1.0e-15))))      # adding small number because psf have 0, and we couldn't divided by 0.

# Wiener filter over several frames
for i=1:30
   psf = pupil_to_psf(aperture, pad(phase[:,:,i]*2.72219,npad));
   image_spec = conv_psf_obj(psf, obj);
   wiener_numer += fft(image_spec).*conj(fft(psf))
   wiener_denom += abs2.(fft(psf)) +1.0e-15
end
wiener_obj =real(fftshift(ifft(wiener_numer./wiener_denom)))
