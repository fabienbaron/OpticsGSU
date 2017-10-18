using OptimPack
include("optics.jl")
x_true = read(FITS("jupiter.fits")[1])*1.;
# Generate data
npix = size(x_true,1);
psf = gaussian2d(npix,npix,3);
#psf= zeros(npix,npix); psf[125:131,128]=1.0; psf[128,125:131]=1.0;
x_conv = convolve(psf,x_true);
σ = 1e-5#maximum(x_true)/10.;
x_noisy = x_conv + σ.*randn(size(x_true));

# Function to minimize + gradient w/r object pixels
function epsilon(model, g, data, σ, psf)
    res = (data - convolve(model,psf))/σ;
    f= sum(res.^2);
    g[:] = -2/σ*correlate(res,psf);
    return f;
end

initial_guess = rand(size(x_noisy));
crit = (x,g)->epsilon(x, g, x_noisy, σ, psf);
x_sol = OptimPack.vmlmb(crit, initial_guess, verb=true, lower=0, maxiter=500);
imview4(x_true, x_noisy, x_sol, convolve(x_sol,psf))
