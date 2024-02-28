using LinearAlgebra, LinearInterpolators, TwoDimensional, ProgressMeter

function pupil_support_size(D, pixsize)
    return nextprod([2], ceil(2*D/pixsize))
end

function generate_aperture(diameter1, diameter2, pupil_width_pixels, pupil_sampling, normalize=true)
    # Function for creating a mask on the data representing the the field of
    # view of a telescopic aperture.
	x = [j for i = -pupil_width_pixels/2:pupil_width_pixels/2-1, j = -pupil_width_pixels/2:pupil_width_pixels/2-1] * pupil_sampling;
    aperture = Float32.(((diameter2/2)^2).<=(x.^2 + x'.^2).<=(diameter1/2)^2)
    if normalize == true
    	aperture /= norm(aperture)
    end
    return aperture
end

function generate_aperture_chromatic_steps(N, λmin, λmax; delta_slice = 1) 
    # generates a telescope aperture for multiwav use
    rad_max = div(N, 4)
    rad_min = ceil(Int, rad_max*(λmin/λmax))
    nwavs = round(Int, (rad_max - rad_min)/delta_slice)+1
    rad_min = rad_max - nwavs*delta_slice+delta_slice
    rad = zeros(Float32, nwavs);
    λ = zeros(Float32, nwavs);
    aperture_mask=zeros(Float32, N, N,nwavs);
    #tel_mask=zeros(Bool, dim,dim,nwavs); 
    println("rad_max = ",rad_max)
    println("rad_min = ",rad_min)
    println("nwavs = ",nwavs)
    # increment in pixels for change in radius of pupil
    Threads.@threads for k in 1:nwavs
        rad[k] = rad_min + (k-1)*delta_slice
        λ[k]=λmin*(float(N)/4.0)/rad[k]
        aperture_mask[:,:,k]= generate_aperture(2*rad[k], 0, N, 1.0)
        println("rad = ",rad[k], " λ = ",λ[k])
    end
    return reverse(aperture_mask, dims=3), reverse(λ);
end

@views function generate_anisoplanatic_frozen_flow_phase_extractors(atmosphere, timestamps, patches, N, source_height, pupil_sampling, source_sampling; T=Float32)
    # Frozen flow, isoplanatic, monochromatic, no propagation
    npatches = patches.npatches;
    nlayers = atmosphere.nlayers;
    nλ = atmosphere.nλ
    λmin = minimum(atmosphere.λ)
    nframes = length(timestamps);
    I = Array{TwoDimensionalTransformInterpolator{T, LinearSpline{T, Flat}, LinearSpline{T, Flat}},4}(undef, npatches, nlayers, nλ, nframes)
    ker1 = LinearSpline(T)
    Φ_size = size(atmosphere.phase_screens,1), size(atmosphere.phase_screens,2) # size of a 2D phase screen
    centers_x = patches.isopatch_positions[:,1].+ (N÷2+1)
    centers_y = patches.isopatch_positions[:,2] .+ (N÷2+1)
    for n=1:nframes
            for l=1:nλ
                for j=1:npatches
                    for i=1:nlayers
                        dx = (atmosphere.heights[i]/source_height*source_sampling*centers_x[j] + timestamps[n] * atmosphere.winds[i, 1]* cos(atmosphere.winds[i, 2]/180*pi) ) / pupil_sampling;
                        dy = (atmosphere.heights[i]/source_height*source_sampling*centers_y[j] + timestamps[n] * atmosphere.winds[i, 1]* sin(atmosphere.winds[i, 2]/180*pi) ) / pupil_sampling;
                        # println("iframe = ", n, "\t layer=" , i, "\tdx= ", dx, "\tdy= ", dy)
                        if ( ((dx + N) > Φ_size[1]) || ((dy + N)> Φ_size[2])) || (dx < 0 ) || (dy<0)
                            println("Out of phase screen bounds at frame $n")
                        end
                        scaling = (1.0 - atmosphere.heights[i]/source_height)*atmosphere.λ[l]/λmin
                        R = (1/scaling)*AffineTransform2D{T}() + scaling .* (dx,dy)
                        # Q: isopatches.npix_isopatch_width????
                        I[j,i,l,n] = TwoDimensionalTransformInterpolator((isopatches.npix_isopatch_width, isopatches.npix_isopatch_width), Φ_size, ker1, R)
                    end
                end
            end
    end
   return I
end



@views function generate_isoplanatic_frozen_flow_phase_extractors(atmosphere, timestamps, N, source_height, pupil_sampling; T=Float32)
    # Frozen flow, isoplanatic, monochromatic, no propagation
    nlayers = atmosphere.nlayers;
    nλ = atmosphere.nλ
    λmin = minimum(atmosphere.λ)
    nframes = length(timestamps);
    I = Array{TwoDimensionalTransformInterpolator{T, LinearSpline{T, Flat}, LinearSpline{T, Flat}},3}(undef, nlayers, nλ, nframes)
    ker1 = LinearSpline(T)
    Φ_size = size(atmosphere.phase_screens,1), size(atmosphere.phase_screens,2) # size of a 2D phase screen
    for n=1:nframes
        for l=1:nλ
            for i=1:nlayers
                dx = (timestamps[n] * atmosphere.winds[i, 1]* cos(atmosphere.winds[i, 2]/180*pi)) / pupil_sampling;
                dy = (timestamps[n] * atmosphere.winds[i, 1]* sin(atmosphere.winds[i, 2]/180*pi)) / pupil_sampling;
             #   println("iframe = ", n, "\t layer=" , i, "\tdx= ", dx, "\tdy= ", dy)
                # if ( ((dx + N) > Φ_size[1]) || ((dy + N)> Φ_size[2]))
                #         println("Out of phase screen bounds at frame $n")
                # end
                scaling = (1.0 - atmosphere.heights[i]/source_height)*(atmosphere.λ[l]/λmin)
                R = (1/scaling)*AffineTransform2D{T}() + scaling .* (dx,dy)
                I[i,l,n] = TwoDimensionalTransformInterpolator((N, N), Φ_size, ker1, R)
            end
        end
    end
   return I
end


@views function extract_composite_phases(Φ_layers, I; FTYPE=Float32, verbose=true)
    # Monochromatic, isoplanatic, no propagation
    nlayers = size(I,1)
    nλ = size(I,2) 
    nframes = size(I,3)    
    pupil_width_pixels = I[1].rows[1]
    phases = zeros(FTYPE, pupil_width_pixels, pupil_width_pixels, nλ, nframes);
    if verbose==true
        p = Progress(nframes, desc="Computing composite phases")
    end
    Threads.@threads for n=1:nframes
        for l=1:nλ
            phases[:, :, l, n] = reduce(+,[I[i,l,n]*Φ_layers[:,:,i,l] for i=1:nlayers])
        end 
       if verbose==true
         next!(p)
       end
    end
    if verbose == true
        finish!(p)
    end
    return phases
end


function pupil_to_psf(A_pupil, Φ_pupil, psf_size, delta)
    Np         = 2*size(Φ_pupil,1);   # Pads out the array so there are no aliasing problems
    npatches   = size(Φ_pupil, 3);
    nframes    = size(Φ_pupil, 4);
    psfs = zeros(Float64, psf_size, psf_size, npatches, nframes);
    Threads.@threads for iframe in 1:nframes
                        for ipos in 1:npatches
                            psfs[:,:, ipos, iframe] = punch_out(abs2.(ft2(pad_array(A_pupil.*cis.(Φ_pupil[:,:,ipos,iframe]),Np), delta)), psf_size);
                        end
                    end
    return psfs
end
