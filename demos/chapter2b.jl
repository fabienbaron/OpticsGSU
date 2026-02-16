using OpticsGSU, OptimPackNextGen, PyPlot, LinearAlgebra, SparseArrays

x_true = readfits("./data/jupiter.fits")*1.0;
# Generate data
npix = size(x_true,1);
psf = gaussian2d(npix,npix,4);
#psf= zeros(npix,npix); psf[125:131,128].=1.0; psf[128,125:131].=1.0;
x_conv = convolve(psf,x_true);
σ = maximum(x_conv)/10.;
x_noisy = x_conv + σ.*randn(Float64, size(x_true));
ψ=1e2
# Function to minimize + gradient w/r object pixels
function epsilon(x, g, data, σ, psf)
    res = (data - convolve(x,psf) ) / σ;
    f = sum(res.^2) + 0.5*ψ*norm(x)^2
    g[:] = -2/σ*correlate(res,psf)+ψ*x
    return f;
end

x_init = deepcopy(x_noisy);
crit = (x,g)->epsilon(x, g, x_noisy, σ, psf);
x_sol = OptimPackNextGen.vmlmb(crit, x_init, verb=true, lower=0.0, maxiter=2000); # positivity is imposed here
imview4(x_true, x_conv, x_noisy, x_sol);
println("ϵ = $(crit(x_sol, dummy ))\t dist = $(norm(x_sol - x_true,1))\n");

#
# With regularization
#



o = ones(npix);
D_1D = spdiagm(-1=>-o[1:npix-1],0=>o)
D = [kron(spdiagm(0=>ones(npix)), D_1D) ;  kron(D_1D, spdiagm(0=>ones(npix)))];
DtD= D'D

function epsilon_tv(object, g, data, σ, psf, μ)
    res = (data - convolve(object,psf) ) / σ;
    x= vec(object);
    f= sum(res.^2) + μ*norm(D*x)^2; # chi2 + Tikhonov regularization
    g[:] = -2/σ*correlate(res,psf)+2*μ*reshape(DtD*x, size(object));
    return f;
end
for μ=[1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1]
#μ=0.005
crit = (x,g)->epsilon_tv(x, g, x_noisy, σ, psf, μ);
x_sol = OptimPackNextGen.vmlmb(crit, x_init, verb=true, lower=0, maxiter=2000); # positivity is imposed here
imview4(x_true, x_conv, x_noisy, x_sol);
dummy=similar(x_sol)
println("μ = $(μ)\t ϵ = $(crit(x_sol, dummy ))\t dist = $(norm(x_sol - x_true,1))\n");
readline();
end


# Sobel gradient + autodiff
using Zygote
N = npix
c = N÷2+1;
Gx = zeros(Float64, N,N); 
Gy = zeros(Float64, N,N); 
Gd1 = zeros(Float64, N,N); 
Gd2 = zeros(Float64, N,N); 
#Gtv_x[c-1:c+1, c-1:c+1] = [0 0 0; -1 1 0; 0 0 0]
#Gtv_y[c-1:c+1, c-1:c+1] = [0 -1 0; 0 1 0; 0 0 0]
Gx[c-1:c+1, c-1:c+1] = [ -1 0 1; -2 0 2; -1 0 1]
Gy[c-1:c+1, c-1:c+1] = [-1 -2 -1; 0 0 0; 1 2 1]
Gd1[c-1:c+1, c-1:c+1] = [0 1 2 ;-1 0 1; -2 -1 0]
Gd2[c-1:c+1, c-1:c+1] = [-2 -1 0; -1 0 1; 0 1 2]
∇ = X -> hcat(convolve(X, Gx), convolve(X, Gy), convolve(X, Gd1), convolve(X,Gd2));
function epsilon_sobel(x, data, σ, psf, μ)
    res = (data - convolve(x,psf) ) / σ;
    f= sum(res.^2) + μ*norm(∇(x),2)^2; 
    return f;
end
μ=5e-4
crit = x->epsilon_sobel(x, x_noisy, σ, psf, μ);
x_sol = OptimPackNextGen.vmlmb(crit, x_init, verb=true, autodiff=true, lower=0, maxiter=2000); 
imview4(x_true, x_conv, x_noisy, x_sol);
