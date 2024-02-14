function object_to_fto(object; CTYPE=ComplexF32)
    N = size(object,1)
    nwavs = size(object,3)
    FTobject = zeros(CTYPE, div(N, 2) + 1, N, nwavs)
    for k = 1:nwavs
        FTobject[:, :, k] = CTYPE.(rft2(object[:, :, k]))
    end
    return FTobject
end


function support_box(object)
    if ndims(object)==3
        object = dropdims(sum(abs.(object), dims=3),dims=3)
    else
        object = abs.(object)
    end
    x = vec(sum(object, dims=1)).>0
    y = vec(sum(object, dims=2)).>0
    return findfirst(x):findlast(x), findfirst(y):findlast(y)
end



function generate_binary(N, λ; f1 = 1.0, f2 = 1.0, T1=6000::Float64, T2=4000::Float64)
    # Blackbody data for binary star object
    c = 2.997e8 # speed of light in m/s
    hp = 6.6260693e-34 # Planck"s constant
    kb = 1.380658e-23 # Boltzmann"s constant
    nwavs = length(λ)
    object = zeros(Float64, N, N, nwavs)
    x1 = div(N, 2) + 1
    if (T1 == -1) || (T2 == -1)
        error("Please specific binary temperatures")
    end
    # Generate binary star object with Blackbody behavior
    y1 = div(N, 2) + 20 + 1
    y0 = div(N, 2) - 20 + 1
    for k = 1:nwavs
        object[y1, x1, k] = floor(8 * pi * hp * c / (λ[k])^5 * (1 / (exp(hp * c / (λ[k] * kb * T1)) - 1)))
        object[y0, x1, k] = floor(8 * pi * hp * c / (λ[k])^5 * (1 / (exp(hp * c / (λ[k] * kb * T2)) - 1)))
    end
    return object
end


@views function generate_hyperspectral_object(N, λ; template="./data/sat_template1.fits")
    # Blackbody data for binary star object
    nwavs = length(λ)
    object = zeros(Float64, N, N, nwavs)
    FTobject = zeros(ComplexF64, div(N, 2) + 1, N, nwavs)
    if template == ""
        error("Please specify object material template")
    end
    nmaterials = 6
    materials = Array{Array{Float64}}(undef, nmaterials)
    # Solar Panel (AR coated)
    materials[1] = [9.667504e2, -4.583580e0, 8.008073e-3, -6.110959e-6, 1.73911e-9]
    # Kapton
    materials[2] = [-1.196472E5, 1.460390E3, -7.648047E0, 2.246897E-2, -4.056309E-5, 4.615807E-8, -3.238676E-11, 1.283035E-14, -2.200156E-18]
    # Aluminized Mylar
    materials[3] = [-7.964498E3, 7.512566E1, -2.883463E-1, 5.812354E-4, -6.488131E-7, 3.801200E-10, -9.129042E-14]
    # Kapton - aged:4
    materials[4] = [-7.668973E4, 9.501756E2, -5.055507E0, 1.509819E-2, -2.771163E-5, 3.204811E-8, -2.283223E-11, 9.171536E-15, -1.591887E-18]
    # Aluminized Mylar - aged:4
    materials[5] = [-2.223456E4, 2.305905E2, -1.007697E0, 2.404597E-3, -3.379637E-6, 2.797231E-9, -1.262900E-12, 2.401227E-16]
    # Faint solar panels 
    materials[6] = 1e-2*materials[1];
    # Read in data image (FITS format)
    obj_coeffs = Int.(crop_to(readfits(template), (N, N)))
    spectra =Array(hcat([[sum([materials[i][ll+1] * (λ[k]*1e9).^ll for ll in 0:length(materials[i])-1]) for k = 1:nwavs] for i=1:nmaterials]...)')
    # for i=1:nmaterials
    # plot(λ*1e9, specs[:,i])
    # end
    # ylim(0, 100)
    abundances = zeros(Float64, N, N, nmaterials)
    for i = 1:nmaterials
        indx = findall(obj_coeffs .== i)
        abundances[indx, i] .= 1.0
        for j in indx
            object[j, :] .= spectra[i,:]
        end
    end
    return object, abundances, spectra
end

function crop_to(x, (wantx, wanty))
    nx = size(x,1)
    ny = size(x,2)
    lo_x = div(nx-wantx,2)+1
    hi_x = div(nx+wantx,2)
    lo_y = div(ny-wanty,2)+1
    hi_y = div(ny+wanty,2)
    return x[lo_x:hi_x,lo_y:hi_y]
end
