module OpticsGSU
include("fits.jl")
include("utils_functions.jl")
include("zernikes.jl")
include("convolve.jl")
include("atmosphere.jl")
include("apertures.jl")
include("generate_data.jl")
include("photometry.jl")
include("objects.jl")
include("refraction.jl")
#include("plots.jl");
include("view.jl")

export imview, imsurf, imview_add, imview4, imview2
# General
export pad, ft2, ift2, rft2, irft2, convolve, correlate, convolve_isoplanatic_all, convolve_tr, convolve_fto, convolve_fto_tr, convolve_otf, convolve_otf_tr, pupil_to_psf, pad_array, gaussian2d, cog, recenter
export zernike, circular_aperture
export generate_aperture_chromatic_steps
export ft_phase_screen, ang_spec_multi_prop, ft_sh_phase_screen, generate_aperture, generate_aperture_with_obscuration, generate_pupil, generate_pupil_all_frames
export pupil_to_psf, create_images, generate_phase_screens
export CN2_huffnagel_valley_generalized, cn2_profile_to_r0, cn2_profile_to_θ0
export ftn2, iftn2
export generate_isoplanatic_data, generate_hyperspectral_object, generate_binary
export generate_tel_aperture_chromatic_steps
export phase_screen_size
export get_refraction, air_refractive_index_minus_one
export Atmosphere
export instantiate_atmosphere, get_non_overlapping_punch_outs
export phase_screen_size, pupil_support_size, phase_to_dr0
export loglikelihood_object, loglikelihood_phases, loglikelihood_object_f, loglikelihood_phases_f, loglikelihood_phase_screens_f, loglikelihood_phase_screens_alt_f, loglikelihood_phase_screens, loglikelihood_phases_annealed

export flux, mag_zeropoint
export object_to_fto, convolve_classic, crop_to, phase_to_dr0
# Util
export cog, posmax, entropy, shift_and_add, recenter, bartlett_hann1d, bartlett_hann1d_centered, bartlett_hann2d, gaussian2d, meshgrid, meshrad, meshpol, cart2pol
export set_fourier_plans, convolve_planned, correlate_planned
# Poisson
export add_poisson_noise, poisson_likelihood_image
# Proximal operators
export prox_l0, prox_l0_plus, prox_l1, prox_l2, prox_l2dist, prox_l2sq, prox_l2dist_simplex3d, prox_quad, prox_quad_diag, prox_poisson, prox_mask, prox_mask_cplx, prox_canonical_simplex, prox_z2_minus_o
export prox_structured_sparsity, prox_structured_sparsity_plus, prox_l1_plus, prox_positivity

# Blind deconv
export solve_object_normal_tv, solve_object_poisson_tv, solve_object_poisson_l1, solve_object_poisson_l0, solve_object_poisson_psd, solve_psfs

# Regularizations
export structured_sparsity, TV_functions, TV_functions_rfft, fit_psd_o, create_covar_psd

# Anisoplanatism
export patch_struct, generate_anisopatches, create_point_grid
export convolve_Nagy_single_frames, convolve_Nagy_all_frames, convolve_Nagy_all_frames_alt, convolve_eff_single_frame, convolve_eff_all_frames

# FITS
export readfits, writefits

#Plot materials 
export plot_materials_avg, plot_materials_pix

#criteria_and_gradients
export loglikelihood_gaussian_object_f, loglikelihood_gaussian_object_fg, regularized_loglikelihood_object

export Detector, convert_to_adu, is_saturated, where_saturated
export generate_isoplanatic_frozen_flow_phase_extractors, generate_anisoplanatic_frozen_flow_phase_extractors, extract_composite_phases
end
