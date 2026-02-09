using OpticsGSU
using Statistics
#Read in Dr03 phase screens  (N.B. these screens give D/r0=3 over the entire data cube)
phase_orig=readfits("./data/Dr03_phases.fits.fz");

# Number of pixels across D/r0 ~ 3
npix = size(phase_orig,1)

# Number of frames
nframes = size(phase_orig,3)

# Compute original D/r0 using variance formula
var_screen = var(phase_orig);
dr0_orig = (var_screen/1.0299)^(3/5) # original D/r0 of the *square* phase screen

# The D/r0 over a circular aperture will be smaller - let's check
aperture = circular_aperture(npix=npix,diameter=npix,centered=true);
phase = phase_orig .* aperture;
indx= findall(repeat(aperture,outer=[1,1,nframes]).>0.5);

# Compute D/r0 over the relevant pixels (= within the aperture)
var_pupil = var(phase[indx]);
dr0_circ = (var_pupil/1.0299)^(3/5);

# Now rescale the phase to be consistent with required D/r0
phase = phase  * (dr0_orig/dr0_circ)^(5/6);

# Check we get the desired D/r0 result now
var_pupil = var(phase[indx]);
dr0 = (var_pupil/1.0299)^(3/5)  # <--- we want it to be close to dr0_orig


# Now we will show how the D/r0 is affected by tip/tilt removal
# Perform Zernike decomposition for tip and tilt and remove from phases
# N.B. dot product gives pi x coeff, therefore need to normalize coeff by pi

# Compute second and third Zernike maps (for tip and tilt): center is at [npix/2+1,npix/2+1]
zern1 = zernike(1,npix=npix,diameter=npix,centered=true);
zern2 = zernike(2,npix=npix,diameter=npix,centered=true);
zern3 = zernike(3,npix=npix,diameter=npix,centered=true);

for j = 1:nframes
  a1 = 0. #sum(zern1.*phase[:,:,j])/3.1415926;
  a2 = dot(zern2,phase[:,:,j])/pi;
  a3 = dot(zern3,phase[:,:,j])/pi;
  phase[:,:,j] =  phase[:,:,j] - (a1*zern1 + a2*zern2 + a3*zern3);
end

#Output D/r0 value in section
println("D/r0 over pupil: ", dr0);
#Output standard deviation for sub-section of phase screen in pupil before t/t correction
println("std. dev. of phase: ",sqrt(var_pupil),"  [radians]");
#Output standard deviation for sub-section of phase screen in pupil after t/t correction
var_pupil_tt = var(phase[indx]);
println("std. dev of phase after tip/tilt subtraction: ",sqrt(var_pupil_tt),"  [radians]");

# Output coefficient for tip/tilt corrected phases
 println("Coeff. for relation between variance and D/r0 with tip/tilt removal: ", var_pupil_tt/dr0^(5/3.));
