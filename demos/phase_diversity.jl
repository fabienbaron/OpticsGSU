using OpticsGSU
using Statistics
using LinearAlgebra

# Object
o = rotl90(readfits("./data/jupiter.fits")*1.0)
# o = zeros(256, 256)
# o[110, 129] = 4.0
# o[140, 89] = 1.0
# o[220, 159] = .1


N = size(o,1)
# PSF1 & PSF2 (defocus)
amplitude = circular_aperture(npix=N, diameter=N/4, centered=true)
amplitude = amplitude/norm(amplitude)
phase = -10*zernike(5, npix=N, diameter=N/4)+17*zernike(7, npix=N, diameter=N/4) 
pupil1 = amplitude.*exp.(im*phase)
PSF1 = abs2.(ift2(pupil1)*N)

#σ = 2e-3
σ = 1.0

image1 = convolve(o, PSF1) + σ*randn(N,N)

diversity = 40*zernike(4, npix=N, diameter=N/4)
pupil2 = amplitude.*exp.(im*(phase+diversity))
PSF2 = abs2.(ift2(pupil2)*N)
image2 = convolve(o, PSF2) + σ*randn(N,N)

f = (x,ϕ) -> 0.5*norm((image1 - convolve(x, abs2.(ift2(N*amplitude.*exp.(im*ϕ)))))/σ)^2/(N*N) + 0.5*norm((image2 - convolve(x, abs2.(ift2(N*amplitude.*exp.(im*(ϕ+diversity))))))/σ)^2/(N*N)
f(o,phase)

o_start = (image1+image2)/2
o_start = o_start.*(o_start.>0)
phase_start = zeros(N,N)
f(o_start, phase_start)
g_o = zeros(N,N)
ϕ = copy(phase)


# function fg(x, g_o)
#     pupil1 = amplitude.*exp.(im*ϕ)
#     PSF1 = abs2.(ift2(pupil1)*N)
#     pupil2 = amplitude.*exp.(im*(ϕ+diversity))
#     PSF2 = abs2.(ift2(pupil2)*N)
#     R1 = (image1 - convolve(x, PSF1))/σ 
#     R2 = (image2 - convolve(x, PSF2))/σ 
#     f1 = 0.5*norm(R1)^2/(N*N) 
#     f2 = 0.5*norm(R2)^2/(N*N)
#     g_o[:] = -1/(N*N)*(correlate(R1,PSF1) + correlate(R2,PSF2))
#     return f1+f2
# end

# fg(o_start, g_o)

# x = vmlmb(fg, o_start, lower=0, maxiter=1000, verb=true)

using Zygote, OptimPackNextGen

function f_o(x)
    pupil1 = amplitude.*exp.(im*ϕ)
    PSF1 = abs2.(ift2(pupil1)*N)
    pupil2 = amplitude.*exp.(im*(ϕ+diversity))
    PSF2 = abs2.(ift2(pupil2)*N)
    R1 = (image1 - convolve(x, PSF1))/σ 
    R2 = (image2 - convolve(x, PSF2))/σ 
    f1 = 0.5*norm(R1)^2/(N*N) 
    f2 = 0.5*norm(R2)^2/(N*N)
    return f1+f2
end
x = vmlmb(f_o, o_start, lower=0, maxiter=300, autodiff=true, verb=true)



function tv(x; ϵ=1e-8)
    xijplus1  = circshift(x, (0,-1))
    xiplus1j  = circshift(x, (-1,0))
    d1 = sqrt.((xiplus1j-x).^2  + (xijplus1-x).^2  .+ ϵ^2)
    return sum(d1.-ϵ)
end


function likelihood_and_prior(x)
  return f_o(x) + 5e-7*tv(x)
end


x = vmlmb(likelihood_and_prior, rand(N,N), lower=0, maxiter=300, autodiff=true, verb=true)


function f_ϕ(ϕ)
    pupil1 = amplitude.*exp.(im*ϕ)
    PSF1 = abs2.(ift2(pupil1)*N)
    pupil2 = amplitude.*exp.(im*(ϕ+diversity))
    PSF2 = abs2.(ift2(pupil2)*N)
    R1 = (image1 - convolve(o, PSF1))/σ 
    R2 = (image2 - convolve(o, PSF2))/σ 
    f1 = 0.5*norm(R1)^2/(N*N) 
    f2 = 0.5*norm(R2)^2/(N*N)
    return f1+f2
end

x = vmlmb(f_ϕ, rand(N,N), maxiter=300, autodiff=true, verb=true)
