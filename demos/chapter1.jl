using OpticsGSU, PyPlot, LinearAlgebra
#Create a disc pupil
# Because of the Fourier transform later on, we would like it to be in a double-sized support
pupil_disc = circular_aperture(npix=256, diameter=128, centered=true);

# Create an annulus pupil
pupil_annulus = circular_aperture(npix=256, diameter=128, centered=true)-circular_aperture(npix=256, diameter=64, centered=true);

# Orthogonality check of disc zernikes
nz=20;
npix = 256
zprod = zeros(nz,nz);
for i=1:nz
  for j=1:nz
    Zi = zernike(i, npix=npix, diameter=npix, centered=true)
    Zj = zernike(j, npix=npix, diameter=npix, centered=true)
    zprod[i,j]=dot(Zi,Zj)
  end
end
imview(zprod, title="Zernike Orthogonality on Disc")
#If we want to check deeper, we can go in log plot, but we need to remove tiny <0 values
zprod=abs.(zprod)
imview(log.(zprod), title="Zernike Orthogonality on Disc -- Deep check")

# Orthogonality check of annuli zernikes
nz=20;
npix=256
Z = zeros(npix, npix, nz)
for i=1:nz
  Z[:,:,i] = zernike(i, npix=npix, diameter=npix, centered=true)
end
pupil_annulus_2 = circular_aperture(npix=npix, diameter=npix, centered=true)-circular_aperture(npix=npix, diameter=npix/4, centered=true);
zprod_ann = zeros(nz,nz);
for i=1:nz
  for j=1:nz
    zprod_ann[i,j]=dot(Z[:,:,i],Z[:,:,j].*pupil_annulus_2)
  end
end
imview(zprod_ann, title="Zernike Orthogonality on Annulus")

#Decomposition of a phase into Zernikes
phase=readfits("./data/atmosphere_d_r0_10.fits");
imview(phase,title="Original phase");
npix_phase = (size(phase))[1]
nz = 200; #let's decompose into 50 modes
a = zeros(nz); # here is the array to store the decomposition factors
for i=1:nz
    Zi = zernike(i, npix=npix_phase, diameter=npix_phase, centered=true)
    a[i]=dot(Zi, phase)/dot(Zi, Zi)
end
println("Decomposition coefficients: ", a);
recomposed_phase = zeros(size(phase)) # array of zeros of the same size as the original phase
for i=1:nz
  recomposed_phase += a[i]*zernike(i, npix=npix_phase, diameter=npix_phase, centered=true)
end
imview(recomposed_phase,title="Recomposed phase");


# GOLAY pupils

npix=1024
centers_x=[-0.5,0.5,0]*npix/4
centers_y = [0, 0, 0]*npix/4

diam = 64 #sub-aperture diameter
aperture = zeros(npix,npix)
for i=1:length(centers_x)
 aperture += circular_aperture(npix=npix, diameter=diam, cent_x =(npix+1)/2+centers_x[i], cent_y = (npix+1)/2+centers_y[i])
end
imview(aperture, title="Linear")
aperture = pad(aperture,npix÷2)
psf = abs2.(ift2(aperture));
otf = ft2(psf);
mtf = abs.(otf);
imview(mtf.^.05)

#Golay-3 sub-aperture positions
npix=1024
centers_x=[-0.5,0.5,0]*npix/4
centers_y = [-sqrt(3)/6,-sqrt(3)/6,sqrt(3)/3]*npix/4

diam = 64 #sub-aperture diameter
aperture = zeros(npix,npix)
for i=1:length(centers_x)
 aperture += circular_aperture(npix=npix, diameter=diam, cent_x =(npix+1)/2+centers_x[i], cent_y = (npix+1)/2+centers_y[i])
end
imview(aperture, title="Golay-3")

#Golay-6 sub-aperture positions
centers_x=[1,3/2,0,-1,-1,-1/2]*npix/4
centers_y=[2,-1,-4,-4,2,5]*sqrt(3)/6*npix/4
diam = 128 #sub-aperture diameter
aperture = zeros(npix,npix)
for i=1:length(centers_x)
 aperture += circular_aperture(npix=npix, diameter=diam, cent_x = (npix+1)/2+centers_x[i], cent_y = (npix+1)/2+centers_y[i])
end
imview(aperture, title="Golay-6")
aperture = pad(aperture,npix÷2)
psf = abs2.(ift2(aperture));
otf = ft2(psf);
mtf = abs.(otf);
imview(mtf.^.05)



ntels=5
small_D = 1.0 #m
big_D   = 3.5 #m
npup = 1024
#pixscale_pupil = 0.0069 #m/pix
pixscale_pupil = 100 #pix/m
θ = (0:ntels-1).*360/ntels
centers_x = npup÷2+1 .+ 0.5*big_D*pixscale_pupil*cos.(θ*pi/180)
centers_y = npup÷2+1 .+ 0.5*big_D*pixscale_pupil*sin.(θ*pi/180)
# Precompute Zernikes
Z = zeros(npup, npup, ntels, 3);
for i=1:ntels
    for j=1:3
        Z[:,:,i,j] = zernike(j, npix=npup, diameter=small_D*pixscale_pupil, cent_x = centers_x[i], cent_y = centers_y[i]);
    end
end
aperture = dropdims(sum(Z[:,:,:,1],dims=3), dims=3); aperture /= norm(aperture);
psf = abs2.(ift2(aperture));
otf = ft2(psf);
mtf = abs.(otf);
imview(mtf)









# Making a PSF

npix=1024;
aperture = circular_aperture(npix=npix, diameter=32, centered=true); # max npix/2 because of Nyquist sampling
aperture = aperture./norm(aperture);  # pupil normalization
#phase= zernike(4, npix, npix/2, centered=true);
phase = zeros(npix, npix)
pupil = aperture.*cis.(phase);
psf = abs2.(ift2(pupil)*npix); #the npix factor is for the normalization of the fft
sum(psf) # should be == 1  !
imview(psf, zoom=8, color="Greys") #view psf from the top
clf(); plot(collect(1:npix), psf[div(npix,2)+1,:]); #plot a slice

# Optical Transfer Function
otf = ft2(psf);
# Modulus Transfer Function 
mtf = abs.(otf);

clf(); plot(collect(1:npix), mtf[div(npix,2)+1,:]); #plot a slice
imsurf(mtf) #3d view of the mtf


## Aberrations
##

# Visualize the first 25 Zernikes one by one
npix=512;
aperture = circular_aperture(npix=npix, diameter=32, centered=true); # npix/2 because FFT needs padded pupil by a factor 2
aperture = aperture/norm(aperture);  # pupil normalization
fig = figure("PSF affected by single Zernike mode",figsize=(12,12))
ax = fig.subplots(5,5)
maxpsfzero=1
strehl = 1
for i=1:25
  phase = 100*zernike(i, npix=npix, diameter=64, centered=true);
  pupil = aperture.*cis.(phase);
  psf=abs2.(ift2(pupil)*npix); #the npix factor is for the normalization of the fft
  if i==1
    strehl = 1
    maxpsfzero = maximum(psf)
  else
    strehl = maximum(psf)/maxpsfzero
  end
  subplot(5,5,i)
 # ax.axes.get_xaxis().set_ticks([]);
 # ax.axes.get_yaxis().set_ticks([]);
  imview_add(psf.^.5, zoom=4, color="gist_yarg")
  title("Zernike $i")
  println("Noll: ", i, " Flux : ", sum(psf), " Strehl: ", strehl)
end
