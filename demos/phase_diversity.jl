using OpticsGSU
using Statistics
using LinearAlgebra
#
# Simulate
#

# Object
o =rotl90(readfits("./data/jupiter.fits")*1.0)
N = size(o,1)
# PSF1 & PSF2 (defocus)
amplitude = circular_aperture(npix=N, diameter=N/4, centered=true)#-circular_aperture(npix=N, diameter=N/4, centered=true);
phase = -10*zernike(5, npix=N, diameter=N/4)+17*zernike(7, npix=N, diameter=N/4) 
pupil1 = amplitude.*exp.(im*phase)
PSF1 = abs2.(ift2(pupil1))
image1 = convolve(o, PSF1)

diversity = 40*zernike(4, npix=N, diameter=N/4)
pupil2 = amplitude.*exp.(im*(phase+diversity))
PSF2 = abs2.(ift2(pupil2))
image2 = convolve(o, PSF2)

#
# Reconstruct
#

# Initial values of o -> focal plane image1
#           and phase on PSF1... ? zero


# Residual

# Gradients

#psf1 = abs2.(ift2(amplitude.*exp.(im*phase)))
#psf2 = abs2.(ift2(amplitude.*exp.(im*(phase+diversity)))

ϵ(o, phase) = (o,phase) -> norm(image1 - convolve(o, abs2.(ift2(amplitude.*exp.(im*phase)))))^2 + norm(image2 - convolve(o, abs2.(ift2(amplitude.*exp.(im*(phase+diversity))))))^2

using Zygote #automatic differentiation

gradient(ϵ, o)


#using OptimPackNextGen

#
# Comparison
#