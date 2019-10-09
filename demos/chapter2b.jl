using OpticsGSU, OptimPackNextGen, PyPlot, LinearAlgebra, SparseArrays
# Note: to install OptimPackNextGen, do the following with Julia's package manager:
# add https://github.com/emmt/LazyAlgebra.jl.git https://github.com/emmt/OptimPackNextGen.jl.git

x_true = readfits("./data/jupiter.fits")*1.0;
# Generate data
npix = size(x_true,1);
psf = gaussian2d(npix,npix,5);
#psf= zeros(npix,npix); psf[125:131,128]=1.0; psf[128,125:131]=1.0;
x_conv = convolve(psf,x_true);
σ = maximum(x_conv)/100.;
x_noisy = x_conv + σ.*randn(Float64, size(x_true));

# Function to minimize + gradient w/r object pixels
function epsilon_tik(object, g, data, σ, psf, μ)
    res = (data - convolve(object,psf) ) / σ;
    f = sum(res.^2) + μ*sum(object.^2); # chi2 + Tikhonov regularization
    g[:] = -2/σ*correlate(res,psf) + 2*μ*object;
    return f;
end

x_init = rand(Float64, size(x_noisy));
dummy = similar(x_init);
μ=0.0005
#for μ=[0, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 1]
crit = (x,g)->epsilon_tik(x, g, x_noisy, σ, psf, μ);
x_sol = OptimPackNextGen.vmlmb(crit, x_init, verb=true, lower=0, maxiter=2000); # positivity is imposed here
imview4(x_true, x_conv, x_noisy, x_sol);
println("μ = $(μ)\t ϵ = $(crit(x_sol, dummy ))\t dist = $(norm(x_sol - x_true,1))\n");
#readline();
#end

#
# With another regularization
#
o = ones(npix);
D_1D = spdiagm(-1=>-o[1:npix-1],0=>o)
D = [kron(spdiagm(0=>ones(npix)), D_1D) ;  kron(D_1D, spdiagm(0=>ones(npix)))];
DtD= D'D

function epsilon_tv(object, g, data, σ, psf, μ)
    res = (data - convolve(object,psf) ) / σ;
    x= vec(object);
    f= sum(res.^2) + μ*norm(D*x,2)^2; # chi2 + Tikhonov regularization
    g[:] = -2/σ*correlate(res,psf)+2*μ*reshape(DtD*x, size(object));
    return f;
end
#for μ=[1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1]
μ=0.005
crit = (x,g)->epsilon_tv(x, g, x_noisy, σ, psf, μ);
x_sol = OptimPackNextGen.vmlmb(crit, x_init, verb=true, lower=0, maxiter=2000); # positivity is imposed here
imview4(x_true, x_conv, x_noisy, x_sol);
println("μ = $(μ)\t ϵ = $(crit(x_sol, dummy ))\t dist = $(norm(x_sol - x_true,1))\n");
readline();
#end
