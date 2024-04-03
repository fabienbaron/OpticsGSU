using OpticsGSU, FFTW, LinearAlgebra, Noise
phase = 10*readfits("./data/Dr03_phases.fits.fz")[:,:,1];  # turbulence
npix = 240
aperture = circular_aperture(npix=npix,diameter=npix/2.,centered=true);


# Grid 5x5

aperture_small = circular_aperture(npix=npix,diameter=20,centered=true);

apertures_shack = zeros(5,5,npix,npix)
psfs_shack = zeros(5,5,npix,npix)
for i=1:5
    for j=1:5
        apertures_shack[i,j,:,:] = circshift(aperture_small, (20*(i-3),20*(j-3)))
        psfs_shack[i,j,:,:] = pupil_to_psf(apertures_shack[i,j,:,:].*exp.(im*phase));
    end
end

