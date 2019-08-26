include("optics.jl");

#Read in Dr03 phase screens  (N.B. these screens give D/r0=3 over the entire data cube)
f=FITS("Dr03_phases.fits"); phase_orig = read(f[1]); close(f);

# Number of pixels across Dr03 phase array
npix = size(phase_orig,1)

# Number of phase screens to use (max=500)
jlim=500

#Initialize array
phase=zeros(size(phase_orig));

# Compute second and third Zernike maps (for tip and tilt): center is at [npix/2+1,npix/2+1]
zern1 = zernike(1,npix,npix,centered=true);
zern2 = zernike(2,npix,npix,centered=true);
zern3 = zernike(3,npix,npix,centered=true);

#Only look at phase inside aperture
aperture = circular_aperture(npix,npix,centered=true);
phase = phase_orig .* aperture;

indx= find(repeat(aperture,outer=[1,1,500]).>0.5);
# Compute d/r0 over subsection all phase pixels
std_dev = std(phase[indx]);
dr0 = (std_dev^2/1.03)^(3/5);

#scale phase to be consistent with required d/r0
#phase = phase  * (newdr0/dr0)^(5/6);
#std_dev1 = std(phase[50:200,50:200,1:jlim]);
#dr0=(std_dev1^2/1.03)^(3/5);

# Perform Zernike decomposition for tip and tilt and remove from phases
# N.B. dot product gives pi x coeff, therefore need to normalize coeff by pi
for j = 1:jlim
  a1 = 0. #sum(zern1.*phase[:,:,j])/3.1415926;
  a2 = sum(zern2.*phase[:,:,j])/3.1415926;
  a3 = sum(zern3.*phase[:,:,j])/3.1415926;
  phase[:,:,j] =  phase[:,:,j] - (a1*zern1 + a2*zern2 + a3*zern3);
end

#Output fits image of corrected phases
f=FITS("corrected.fits","w"); write(f,phase); close(f);

#Output D/r0 value in section
println("D/r0 over pupil: ",dr0);
#Output standard deviation for sub-section of phase screen in pupil before t/t correction
println("std. dev. of phase: ",std_dev,"  [radians]");
#Output standard deviation for sub-section of phase screen in pupil after t/t correction
std_dev2 = std(phase[indx]);
println("std. dev of phase after tip/tilt subtraction: ",std_dev2,"  [radians]");

# Output coefficient for tip/tilt corrected phases
 println("Coeff. for relation between variance and D/r0 with tip/tilt removal: ",std_dev2^2/dr0^(5/3.));
