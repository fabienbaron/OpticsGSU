using OptimPack
include("optics.jl")
obj = read(FITS("../opticscourse/jupiter.fits")[1])*1.;
# Generate data
npix = size(obj,1);
psf = gaussian2d(npix,npix,5);
image_true = convolve(psf,obj);
σ = maximum(image_true)/10.;
image_noisy = image_true + σ.*randn(size(obj));

# Function to minimize + gradient w/r object pixels
function epsilon(model, g, data, σ, psf)
    res = (data - convolve(model,psf))/σ;
    f= sum(res.^2);
    g[:] = -2/σ*correlate(res,psf);
    return f;
end

initial_guess = rand(size(image_noisy));#copy(image_noisy);
crit = (x,g)->epsilon(x, g, image_noisy, σ, psf);
x_sol = OptimPack.vmlmb(crit, initial_guess, verb=true, lower=0, maxiter=500);



# numerical gradient
# model = copy(initial_guess)
# g= -2/σ*correlate((image_noisy - convolve(psf,model))/σ,psf)
# res0 = sum( ((image_noisy - convolve(psf,model))/σ).^2);
# model[23,56] += 1e-6
# res1 = sum( ((image_noisy - convolve(psf,model))/σ).^2);
# println("Numerical: ", (res1-res0)/1e-6, " Analytic: " ,  g[23,56])
