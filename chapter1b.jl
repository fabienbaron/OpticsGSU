include("optics.jl");
obj=read(FITS("jupiter.fits")[1])*1.0
npix=size(obj,1)

# Exercice 1 : create Gaussian PSF with sigma=5 pixels
psf_gaussian=gaussian2d(npix,npix,5)
image = conv_psf_obj(psf_gaussian, obj);




# Exercice 2: Read phase screen
phase=read((FITS("Dr03_phases.fits")[1]));  # turbulence
aperture=circular_aperture(npix,npix/2,centered=true);

wiener_numer = zeros(Complex64,npix,npix);
wiener_denom = zeros(Complex64,npix,npix);

npad = div(max(npix-size(phase,1),0),2); # we may need to pad the phase screen

for i=1:30
   psf = pupil_to_psf(aperture, pad(phase[:,:,i]*2.72219,npad));
   image_spec = conv_psf_obj(psf, obj);
   wiener_numer += fft(image_spec).*conj(fft(psf))
   wiener_denom += fft(psf).*conj(fft(psf)) +1.0e-15
end



# Direct inversion
direct_obj=real(fftshift(ifft(fft(image_spec)./(fft(psf_sum).+1.0e-15))))      # adding small number because psf have 0, and we couldn't divided by 0.

# Wiener filter
wiener_obj =real(fftshift(ifft(wiener_numer./wiener_denom)))
