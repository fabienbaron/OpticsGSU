using OpticsGSU
using OptimPackNextGen
using PyPlot

x_true = readfits("./data/jupiter.fits")*1.0;
# Generate data
npix = size(x_true,1);
psf = gaussian2d(npix,npix,3);
#psf= zeros(npix,npix); psf[125:131,128]=1.0; psf[128,125:131]=1.0;
x_conv = convolve(psf,x_true);
σ = maximum(x_conv)/10.;
x_noisy = x_conv + σ.*randn(Float64, size(x_true));

# Function to minimize + gradient w/r object pixels
function epsilon(model, g, data, σ, psf, μ)
    res = (data - convolve(model,psf))/σ;
    f= sum(res.^2) + μ*sum(model.^2); # chi2 + regularization
    g[:] = -2/σ*correlate(res,psf)+2*μ*model;
    return f;
end

x_init = rand(Float64, size(x_noisy));

for μ=[1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1,1e2,1e3,1e4]
crit = (x,g)->epsilon(x, g, x_noisy, σ, psf, μ);
x_sol = OptimPackNextGen.vmlmb(crit, x_init, verb=true, lower=0, maxiter=2000);
imview4(x_true, x_noisy, x_sol, convolve(x_sol,psf));
readline();clf();
end
