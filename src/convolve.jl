using FFTW, LinearAlgebra

function conv_psf_obj(psf::Array{Float64,2}, object::Array{Float64,2})
return real.(fftshift(ifft(fft(object).*fft(psf))));
end

function convolve(a::Array{Float64,2}, b::Array{Float64,2})
return real.(fftshift(ifft(fft(a).*fft(b))));
end

function correlate(a::Array{Float64,2}, b::Array{Float64,2})
return real.(fftshift(ifft(fft(a).*conj(fft(b)))));
end

function pupil_to_psf(pup_amplitude, pup_phase)
    pupil= pup_amplitude.*cis.(pup_phase)
    pupil ./= norm(pupil);
    psf = abs2.(ifft(pupil)*size(pupil,1));
    return fftshift(psf);
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

using FITSIO
function readfits(fitsfile; normalize = false, vectorize=false)
if fitsfile[end-1:end]=="fz"
x = (read((FITS(fitsfile))[2]))
else
    x = (read((FITS(fitsfile))[1]))
end
if normalize == true
 x ./= sum(x)
end
if vectorize == true
    x = vec(x)
end
return x;
end

function writefits(data, fitsfile)
f = FITS(fitsfile, "w");
write(f, data);
close(f);
end
