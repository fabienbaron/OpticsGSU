using OpticsGSU;
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
    zprod[i,j]=sum(Zi.*Zj)
  end
end
imview(zprod, title="Zernike Orthogonality on Disc")
#If we want to check deeper, we can go in log plot, but we need to remove tiny <0 values
zprod=abs.(zprod)
imview(log.(zprod), title="Zernike Orthogonality on Disc -- Deep check")

# Orthogonality check of annuli zernikes
nz=20;
npix=256
pupil_annulus_2 = circular_aperture(npix=npix, diameter=npix, centered=true)-circular_aperture(npix=npix, diameter=npix/4, centered=true);
zprod_ann = zeros(nz,nz);
for i=1:nz
  for j=1:nz
    Zi = zernike(i, npix=npix, diameter=npix, centered=true)
    Zj = zernike(j, npix=npix, diameter=npix, centered=true)
    zprod_ann[i,j]=sum(Zi.*Zj.*pupil_annulus_2)
  end
end
imview(zprod_ann, title="Zernike Orthogonality on Annulus")

#Decomposition of a phase into Zernikes
using FITSIO
phase=read((FITS("./data/atmosphere_d_r0_10.fits"))[1]);
imview(phase,title="Original phase");
npix_phase = (size(phase))[1]
nz = 50; #let's decompose into 20 modes
a = zeros(nz); # here is the array to store the decomposition factors
for i=1:nz
    Zi = zernike(i, npix_phase, npix_phase, centered=true)
    a[i]=sum(Zi.*phase)
end
println("Decomposition coefficients: ", a);
recomposed_phase = zeros(size(phase)) # array of zeros of the same size as the original phase
for i=1:nz
   recomposed_phase += a[i]*zernike(i, npix_phase, npix_phase, centered=true)
end
imview(recomposed_phase,title="Recomposed phase");


# GOLAY pupils

#Golay-3 sub-aperture positions
npix=256
centers_x=[-0.5,0.5,0]*npix/4
centers_y = [-sqrt(3)/6,-sqrt(3)/6,sqrt(3)/3]*npix/4

diam = 64 #sub-aperture diameter
aperture = zeros(npix,npix)
for i=1:length(centers_x)
 aperture += circular_aperture(npix, diam, (npix+1)/2+centers_x[i], (npix+1)/2+centers_y[i])
end
imview(aperture, title="Golay-3")

#Golay-6 sub-aperture positions
centers_x=[1,3/2,0,-1,-1,-1/2]*npix/4
centers_y=[2,-1,-4,-4,2,5]*sqrt(3)/6*npix/4

diam = 64 #sub-aperture diameter
aperture = zeros(npix,npix)
for i=1:length(centers_x)
 aperture += circular_aperture(npix, diam, (npix+1)/2+centers_x[i], (npix+1)/2+centers_y[i])
end
imview(aperture, title="Golay-6")

# Making a PSF
npix=1024;
aperture = circular_aperture(npix=npix, diameter=npix/16, centered=true); # npix/2 because FFT needs padded pupil by a factor 2
aperture = aperture/norm(aperture);  # pupil normalization
#phase= zernike(4, npix, npix/2, centered=true);
phase = 0
pupil=aperture.*cis.(phase);
psf=abs2.(ifft(pupil)*npix); #the npix factor is for the normalization of the fft
psf = fftshift(psf); # fft is centered on [1,1], but we want it on npix/2,npix/2
sum(psf) # should be == 1  !
imview(psf, zoom=8, color="Greys") #view psf from the top
clf()
plot(collect(1:npix), psf[div(npix,2),:]); #plot a slice

## Aberrations
##

# Visualize the first 16 Zernikes one by one
npix=512;
aperture = circular_aperture(npix=npix, diameter=npix/16, centered=true); # npix/2 because FFT needs padded pupil by a factor 2
aperture = aperture/norm(aperture);  # pupil normalization
fig = figure("PSF affected by single Zernike mode",figsize=(12,12))
maxpsfzero=1
strehl = 1
for i=1:16
  phase= 20.*zernike(i, npix, npix/16, centered=true);
  pupil=aperture.*cis.(phase);
  psf=abs2.(ifft(pupil)*npix); #the npix factor is for the normalization of the fft
  psf = fftshift(psf); # fft is centered on [1,1], but we want it on npix/2,npix/2
  if i==1
    strehl = 1
    maxpsfzero = maximum(psf)
  else
    strehl = maximum(psf)/maxpsfzero
  end

  ax=fig[:add_subplot](4,4,i)
  ax[:axes][:get_xaxis]()[:set_ticks]([]);
  ax[:axes][:get_yaxis]()[:set_ticks]([]);
  imview_add(psf, zoom=4, color="Greys")
  title("Zernike $i")
  println("Noll: ", i, " Flux : ", sum(psf), " Strehl: ", strehl)
end

# OTF
otf = fftshift(fft(psf)); #fft result always need to be shifted
mtf = abs.(otf);
imsurf(mtf) #3d view of the mtf
