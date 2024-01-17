using FFTW, LinearAlgebra

function ftn2(x)
    return fftshift(fft(fftshift(x)))/sqrt(length(x))
end
   
function iftn2(x)
   return ifftshift(ifft(ifftshift(x)))/sqrt(length(x))
end

function ft2(x)
    return fftshift(fft(fftshift(x)))
end
   
function ift2(x)
   return ifftshift(ifft(ifftshift(x)))
end
   
function rft2(x)
    return fftshift(rfft(fftshift(x)))
end
   
function irft2(x,n)
    return ifftshift(irfft(ifftshift(x),n))
end

function irft2(x, n, dims::Tuple)
    return ifftshift(irfft(ifftshift(x,dims),n, dims),dims)
end

function rft2(x, dims::Tuple)
    return fftshift(rfft(fftshift(x,dims),dims),dims)
end

function ft2(g, δ)
    # Function for calculating a two dimensional fourier transfer
    return ft2(g) * δ^2;
end

function ift2(G, δ_f)
    # Function for calculating a two dimensional inverse fourier transfer
    N = size(G, 1);
	return ift2(G) * (N*δ_f)^2;
end

function corr2_ft(u1,u2, mask, δ)
    # Function for calculating the correlation between two arrays
    N = size(u1, 1);
    c = zeros(N,N) + im*zeros(N,N);
    δ_f = 1/(N*δ); #frequency grid spacing [m]
    U1 = ft2(u1 .* mask, δ);
    U2 = ft2(u2 .* mask, δ);
    U12corr = ift2(conj(U1).*U2, δ_f);
    maskcorr = ift2(abs.(ft2(mask, δ)).^2, δ_f) * δ^2;
    idx = find(maskcorr.!=0);
    c[idx] = U12corr[idx] ./ maskcorr[idx] .* mask[idx];
    return c
end

function convolve_isoplanatic_all(object::Array{Float64,2}, psfs::Array{Float64,3})
    images= Array{Float64, 3}(undef, size(object,1), size(object,2), size(psfs,3))
    Threads.@threads for i=1:size(psfs,3)
        images[:,:,i]= convolve(object, (@view psfs[:,:,i]))
    end
    return images
end

# function convolve_classic(a::Array{Float64,2}, b)
#     return real(ift2(ft2(a).*ft2(b)))
# end
  
function convolve(object, psf)
    return irft2(rft2(object).*rft2(psf),size(object,1))
end

function correlate(object, psf)
    return irft2(rft2(object).*conj(rft2(psf)),size(object,1))
end


function set_fourier_plans(direct, inverse)
    F = plan_rfft(direct, flags=FFTW.ESTIMATE, timelimit=Inf)
    iF = plan_irfft(inverse, size(inverse,2), flags=FFTW.ESTIMATE, timelimit=Inf)
    return F, iF
 end
 
# function convolve_planned(object::Array{Float32,2}, otf::Array{ComplexF32,2}, F::FFTW.rFFTWPlan, iF)
#     return ifftshift(iF*(ifftshift(fftshift(F*fftshift(object)).*otf)))
# end

function convolve_planned(object, otf, F::FFTW.rFFTWPlan, iF)
   return ifftshift(iF*(ifftshift(fftshift(F*fftshift(object)).*otf)))
end

function correlate_planned(object, otf, F::FFTW.rFFTWPlan, iF)
    return ifftshift(iF*(ifftshift(fftshift(F*fftshift(object)).*conj(otf))))
end
 


function convolve(object::Array{Float64,2}, psf::SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true})
    return irft2(rft2(object).*rft2(psf),size(object,1))
end
      
function convolve_pow(object::Array{Float64,2}, psf::Array{Float64,2})
    return irft2(abs2.(rft2(object)).*rft2(psf),size(object,1))
end

function convolve_tr(object::Array{Float64,2}, psf::Array{Float64,2})
    return irft2(conj(rft2(object)).*rft2(psf),size(object,1))  
end

function convolve_otf(object::Array{Float64,2}, otf::Array{ComplexF64,2})
    return real.(ift2(ft2(object).*otf))
end
   
function convolve_otf_tr(object::Array{Float64,2}, otf::Array{ComplexF64,2})
    return real.(ift2(ft2(object).*conj.(otf)))  
end

 
function convolve_otf(object::Array{Float64,3}, otf::Array{ComplexF64,3})
    return irft2(rft2(object, (1,2)).*otf,size(object,2))
end
   
function convolve_otf_tr(object::Array{Float64,3}, otf::Array{ComplexF64,3})
    return irft2(rft2(object, (1,2)).*conj.(otf),size(object,2))  
end

      
function convolve_fto(ft_object::Matrix{ComplexF64}, psf::Matrix{Float64})
    return irft2(ft_object.*rft2(psf),size(ft_object,2))
end
   
      
function convolve_fto(ft_object::Matrix{ComplexF32}, psf::Matrix{Float64})
    return irft2(ft_object.*rft2(psf),size(ft_object,2))
end
   
function convolve_fto_tr(ft_object::Array{ComplexF64,2}, psf::Array{Float64,2})
    return irft2(conj(ft_object).*rft2(psf),size(ft_object,2))  
end

function convolve_gaussblur(x::Array{Float64,2}, fwhm )
    nx = size(x,1)
    ny = size(x,2)
    psf = gaussian2d(nx,ny,fwhm)
    return convolve(x, psf)
end

function pupil_to_psf(pup_amplitude, pup_phase; normalize = false) 
    pupil= pup_amplitude.*cis.(pup_phase)
    return pupil_to_psf(pupil, normalize = normalize);
end

function pupil_to_psf(pupil; normalize=false)
    if normalize == true
        pupil /= norm(pupil)
    end
    return abs2.(ift2(pupil)*size(pupil,1))
end

function pad(x::Array{Float64,2}, nborder)
 npad = size(x,1)+nborder*2;
 xpadded = zeros(Float64, npad,npad);
 xpadded[nborder+1:npad-nborder,nborder+1:npad-nborder] .= x;
 return xpadded
end

function pad(x::Array{Complex{Float64},2}, nborder)
 npad = size(x,1)+nborder*2;
 xpadded = zeros(Complex{Float64}, npad,npad);
 xpadded[nborder+1:npad-nborder,nborder+1:npad-nborder] .= x;
 return xpadded
end

function pad_array(array0, size_desired)
   # Function to pad out an array to the desired size.
   size0 = size(array0,1);
   size_diff = size_desired - size0
   size_pad = div(size_diff,2);
   size1 = size0 + 2*size_pad;
   if size_diff%2 == 0
       array_return = hcat(zeros(size1,size_pad), vcat(zeros(size_pad,size0), array0, zeros(size_pad,size0)), zeros(size1,size_pad));
   else
       array_return = hcat(zeros(size1-1,size_pad), vcat(zeros(size_pad,size0), array0, zeros(size_pad-1,size0)), zeros(size1-1,size_pad-1));
   end
   return array_return
end

function pupil_to_psf(pup_amplitude, pup_phase)
    pupil= pup_amplitude.*cis.(pup_phase)
    pupil ./= norm(pupil);
    psf = abs2.(ifft(pupil)*size(pupil,1));
    return ifftshift(psf);
end