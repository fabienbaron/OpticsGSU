using OpticsGSU, PyPlot, LinearAlgebra

D = 0.7 #m 
lambda = collect(range(400.0, 1000.0, step=50))*1e-9 #m

npix = 256
nyquist_pixels = 128
lambda_min = minimum(lambda)

s = nyquist_pixels*lambda_min./lambda

obj_ref=rotl90(readfits("./data/jupiter.fits")*1.0) # assumed observed at lambda_ref=600nm
lambda_ref = 600e-9
obj = zeros(npix, npix, length(lambda))
pupil = zeros(npix, npix, length(lambda))
psf = zeros(npix, npix, length(lambda))

hyperspectral_images = zeros(npix, npix, length(lambda))

for i=1:length(lambda)
    pupil[:,:,i] = circular_aperture(npix=npix, diameter=s[i], centered=true)
    psf[:,:,i] = pupil_to_psf(pupil[:,:,i], normalize=true)
    obj[:,:,i] = obj_ref*(lambda_ref/lambda[i])^-4
    hyperspectral_images[:,:,i]= convolve(obj[:,:,i], psf[:,:,i])
end
broadband_image = sum(hyperspectral_images, dims=3)


