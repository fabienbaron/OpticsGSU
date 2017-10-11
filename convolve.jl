function conv_psf_obj(psf::Array{Float64,2}, object::Array{Float64,2})
return real(fftshift(ifft(fft(object).*fft(psf))));
end


function pupil_to_psf(pup_amplitude, pup_phase)
    pupil= pup_amplitude.*cis.(pup_phase)
    pupil /= norm(pupil);
    psf = abs2.(ifft(pupil)*size(pupil,1));
    return fftshift(psf);
end


function pad(x::Array{Float64,2}, nborder)
 npad = size(x,1)+nborder*2;
 xpadded = zeros(npad,npad);
 xpadded[nborder+1:npad-nborder,nborder+1:npad-nborder] = x;
 return xpadded
end

function pad(x::Array{Complex64,2}, nborder)
 npad = size(x,1)+nborder*2;
 xpadded = zeros(Complex64,npad,npad);
 xpadded[nborder+1:npad-nborder,nborder+1:npad-nborder] = x;
 return xpadded
end
