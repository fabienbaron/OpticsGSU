using ProgressMeter
import Distributions

function poisson_likelihood_image(model, data)
return -sum(loglikelihood.(Distributions.Poisson.(model), data))
end

function add_poisson_noise(x)
  return rand.(Distributions.Poisson.(x))  
end

mutable struct Detector
   poisson::Bool
   adu::Bool
   aduTYPE::DataType
   qe::Float32 #Array{Float32,1}
   gain::Float32
   saturation::Int32
   σ_ron::Float32
   exptime::Float32
end

function convert_to_adu(x, detector)
   x = floor.(detector.aduTYPE, min.(max.(detector.gain*x, 0), detector.saturation));
   return x; 
end

function is_saturated(x, detector)
   return sum(x.>=detector.saturation).>0
end

function where_saturated(x, detector)
   return findall(x.>=detector.saturation)
end

@views function generate_isoplanatic_data(t, λ, I, FTobject, aperture_mask, atmosphere, detector::Detector; zenith = 0,  pixscale = 1e99, disperse = false, λref_disp = 500e-9, FTYPE=Float32, CTYPE=ComplexF32)
   # Note: λref_disp is the reference wavelength for zero dispersion
   # DEBUG: t=timestamps; zenith = z;  disperse= true;  FTYPE=Float32; CTYPE=ComplexF32;λref_disp = 500e-9; noise=true
   dim       = size(aperture_mask, 1)
   nwavs     = atmosphere.nλ
   nframes   = length(t)
   psf_broad  = zeros(FTYPE, dim, dim, nframes);
   otfs       = zeros(CTYPE, div(dim,2)+1, dim, nwavs, nframes);
   data       = zeros(FTYPE, dim, dim, nframes); 
   pupil_phases = extract_composite_phases(atmosphere.phase_screens, I); # Frozen flow model, composite phase = add phases
   pupil_amps = zeros(FTYPE, dim, dim, nwavs, nframes);
     
   θref = 0;
   if disperse == true
     θref = get_refraction(λref_disp, zenith);
   end 
   p = Progress(nframes, desc="Computing images")
   Threads.@threads for n=1:nframes
        for k=1:nwavs
           pupil_amps[:,:,k,n] = aperture_mask[:,:,k]; #.*M[n](pupil_amp,k);
           psftmp = abs2.(ift2(pupil_amps[:,:,k,n].*cis.(pupil_phases[:,:,k,n]))*dim)
           if disperse == true
               delta = zeros(FTYPE, dim,dim)
               θk  = get_refraction(λ[k], zenith);
               dpix = 206265.0*(θk-θref)/pixscale #Atmospheric differential refraction, convert from rad to arcsec
               delta[round(Int, div(dim,2)+1-dpix),div(dim,2)+1]= 1.0 
               psftmp = convolve(delta, psftmp)
           end
           psf_broad[:,:,n] +=  psftmp*real(FTobject[div(div(dim,2)+1,2)+1,div(dim,2)+1,k]) # psf_λ x sum(obj_λ)
           otfs[:,:,k,n]    = rft2(psftmp)
           data[:,:,n] += irft2(FTobject[:,:,k].*otfs[:,:,k,n], dim) 
       end
       if detector.poisson == true
         data[:,:,n] = detector.qe*add_poisson_noise.(data[:,:,n]) + detector.σ_ron*randn(FTYPE, dim, dim)
       end
       next!(p)
   end       
   finish!(p)

   # Conversion to ADU if needed
   if detector.adu == true 
      data = convert_to_adu(data, detector)
   else
      data = max.(data, FTYPE(0))
   end
   return data, psf_broad, otfs, pupil_amps, pupil_phases
end


@views function generate_liom_data(offsets, t, λ, I, aperture_mask, atmosphere, detector::Detector; zenith = 0,  pixscale = 1e99, disperse = false, λref_disp = 500e-9, FTYPE=Float32, CTYPE=ComplexF32)
   # Note: λref_disp is the reference wavelength for zero dispersion
   # DEBUG: t=timestamps; zenith = z;  disperse= true;  FTYPE=Float32; CTYPE=ComplexF32;λref_disp = 500e-9; noise=true
   dim       = size(aperture_mask, 1)
   nwavs     = atmosphere.nλ
   nframes   = length(t)
   psf       = zeros(FTYPE, dim, dim, nframes);
   data       = zeros(FTYPE, dim, dim, nframes); 
   pupil_phases = extract_composite_phases(atmosphere.phase_screens, I); # Frozen flow model, composite phase = add phases
   pupil_amps = zeros(FTYPE, dim, dim, nwavs, nframes);
   θref = 0;
   if disperse == true
     θref = get_refraction(λref_disp, zenith);
   end 
   p = Progress(nframes, desc="Computing images")
   Threads.@threads for n=1:nframes
      for k=1:nwavs
           data[:,:,n] += abs2.(ift2(aperture_mask[:,:,k].*cis.(offsets[:,:,k,n]+pupil_phases[:,:,k,n]))*dim)
      end
       next!(p)
   end       
   finish!(p)
   return data, psfs
end


