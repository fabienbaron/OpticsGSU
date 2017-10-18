using OptimPack
include("optics.jl")
obj = read(FITS("../opticscourse/jupiter.fits")[1])*1.;
# Generate data
npix = size(obj,1);
psf = zeros(npix,npix); psf[125:131,125:131]=1;
x_true = convolve(psf,obj);
σ = maximum(x_true)/10.;
x_noisy = x_true + σ.*randn(size(obj));

# Function to minimize + gradient w/r object pixels
function epsilon(model, g, data, σ, psf)
    res = (data - convolve(model,psf))/σ;
    f= sum(res.^2);
    g[:] = -2/σ*correlate(res,psf);
    return f;
end

initial_guess = rand(size(x_noisy));#copy(image_noisy);
crit = (x,g)->epsilon(x, g, x_noisy, σ, psf);
x_sol = OptimPack.vmlmb(crit, initial_guess, verb=true, lower=0, maxiter=500);
imview4(x_true, x_noisy,x_sol, convolve(x_sol,psf))
