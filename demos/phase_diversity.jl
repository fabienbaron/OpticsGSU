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
amplitude = circular_aperture(npix=N, diameter=N/4, centered=true)
amplitude = amplitude/norm(amplitude)

phase = -10*zernike(5, npix=N, diameter=N/4)+17*zernike(7, npix=N, diameter=N/4) 
pupil1 = amplitude.*exp.(im*phase)
PSF1 = abs2.(ift2(pupil1)*N)
σ = 5.0
image1 = convolve(o, PSF1) + σ*randn(N,N)

diversity = 40*zernike(4, npix=N, diameter=N/4)
pupil2 = amplitude.*exp.(im*(phase+diversity))
PSF2 = abs2.(ift2(pupil2)*N)
image2 = convolve(o, PSF2) + σ*randn(N,N)

#
# Reconstruct
#

# Initial values of o -> focal plane image1
#           and phase on PSF1... ? zero


f = (x,ϕ) -> 0.5*norm((image1 - convolve(x, abs2.(ift2(N*amplitude.*exp.(im*ϕ)))))/σ)^2/(N*N) + 0.5*norm((image2 - convolve(x, abs2.(ift2(N*amplitude.*exp.(im*(ϕ+diversity))))))/σ)^2/(N*N)

f(o,phase)
# Residual

# Gradients

#psf1 = abs2.(ift2(amplitude.*exp.(im*phase)))
#psf2 = abs2.(ift2(amplitude.*exp.(im*(phase+diversity)))
o_start = (image1+image2)/2
o_start = o_start.*(o_start.>0)
phase_start = zeros(N,N)
f(o_start, phase_start)

g_o = zeros(N,N)


function fg(x, ϕ, g_o)
    f1 = 0.5*norm((image1 - convolve(x, abs2.(ift2(N*amplitude.*exp.(im*ϕ)))))/σ)^2/(N*N) 
    f2 = 0.5*norm((image2 - convolve(x, abs2.(ift2(N*amplitude.*exp.(im*(ϕ+diversity))))))/σ)^2/(N*N)
    pupil1 = amplitude.*exp.(im*ϕ)
    PSF1 = abs2.(ift2(pupil1)*N)
    pupil2 = amplitude.*exp.(im*(ϕ+diversity))
    PSF2 = abs2.(ift2(pupil2)*N)
    R1 = (image1 - convolve(x, PSF1))/σ 
    R2 = (image2 - convolve(x, PSF2))/σ 
    g_o[:] = -2*correlate(R1,PSF1) -2*correlate(R2,PSF2)
    return f1+f2
end

fg(o_start, phase_start, g_o)

using OptimPackNextGen

x = vmlmb(fg, o_start, lower=0)


#
# Comparison
#