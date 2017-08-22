#
# Based on poppy.zernike
#
function zernike_rad(n, m, rho)
"""Compute R[n, m], the Zernike radial polynomial

    Parameters
    ----------
    n, m : int
        Zernike function degree
    rho : array
        Image plane radial coordinates. `rho` should be 1 at the desired pixel radius of the unit circle
"""
    m = Int(abs(m))
    n = Int(abs(n))
    output = zeros(size(rho))
    if isodd(n - m)
        return 0
    else
        for k=0:div(n - m,2)
            output += ((-1)^k * factorial(n - k) / (factorial(k) * factorial(div(n + m,2) - k) * factorial(div(n - m, 2) - k))) * rho.^(n - 2 * k)
        end
        return output
    end
end

function zernike_nm(n, m, npix=128, diameter=128, cent_x=64, cent_y=64, outside=0, noll_normalize=true)
"""Return the Zernike polynomial Z[m,n] for a given pupil.

    For this function the desired Zernike is specified by 2 indices m and n.
    See zernike_nm for an equivalent function in which the polynomials are
    ordered by n,m

    You may specify the pupil in one of two ways:
     zernike(n, m, npix)       where npix specifies a pupil diameter in pixels.
                               The returned pupil will be a circular aperture
                               with this diameter, embedded in a square array
                               of size npix*npix.

    The expressions for the Zernike terms follow the normalization convention
    of Noll et al. JOSA 1976 unless the `noll_normalize` argument is False.

    Parameters
    ----------
    n, m : int
        Zernike function degree
    npix: int
        Number of pixels of the support.
    diameter:
        Desired pixel diameter for the circular pupil
    noll_normalize : bool
        As defined in Noll et al. JOSA 1976, the Zernike definition is
        modified such that the integral of Z[n, m] * Z[n, m] over the
        unit disk is pi exactly. To omit the normalization constant,
        set this to False. Default is True.

    Returns
    -------
    zern : 2D numpy array
        Z(m,n) evaluated at each (rho, theta)
"""

if n<m
    error("Zernike index m must be >= index n")
else
    x = (collect(linspace(1,npix,npix)) - cent_x) / (diameter / 2.)
    y = (collect(linspace(1,npix,npix)) - cent_y) / (diameter / 2.)
    xx = repmat(x,1,npix)
    yy = repmat(y,1,npix)'
    rho = sqrt.(xx.^2 + yy.^2)
    theta = atan2.(yy, xx)
    aperture = ones(size(rho))
    aperture[find(rho.>1)] = outside  # this is the aperture mask
    norm_coeff=1.0
    if (noll_normalize==true)
        norm_coeff = 1./(2*(n+1)/(1+(m==0)))^0.5
    end
    if (m > 0)
        return norm_coeff*zernike_rad(n, m, rho).*cos.(m * theta).*aperture
    end
    if (m < 0)
        return norm_coeff*zernike_rad(n, m, rho).*sin.(-m * theta).* aperture
    end
    return norm_coeff*zernike_rad(n, 0, rho).*aperture
end
end

function noll_indices(j)
    """Convert from 1-D to 2-D indexing for Zernikes or Hexikes.

    Parameters
    ----------
    j : int
        Zernike function ordinate, following the convention of Noll et al. JOSA 1976.
        Starts at 1.

    """

    if j < 1
        error("Zernike index j must be a positive integer")
    end
    n = 0
	j1 = j-1
	while (j1 > n)
		n += 1
		j1 -= n
    end
	m = (-1)^j * (mod(n,2) + 2 * div(j1+mod(n+1,2),2))
    return (n, m)
end



function zernike(j, x...)
"""Return the Zernike polynomial Z[j] for a given pupil.
For this function the desired Zernike is specified by the Noll index j.
See zernike_nm for an equivalent function in which the polynomials are ordered by n,m
The expressions for the Zernike terms follow the normalization convention
of Noll et al. JOSA 1976 unless the `noll_normalize` argument is False.

        Parameters
        ----------
        j : int
            Noll order
        npix: int
            Number of pixels of the support.
        diameter:
            Desired pixel diameter for the circular pupil
        noll_normalize : bool
            As defined in Noll et al. JOSA 1976, the Zernike definition is
            modified such that the integral of Z[n, m] * Z[n, m] over the
            unit disk is pi exactly. To omit the normalization constant,
            set this to False. Default is True.

        Returns
        -------
        zern : 2D array Z[j] evaluated at each (rho, theta)
"""
  n, m = noll_indices(j);
  return zernike_nm(n,m,x...);
end

using PyPlot

function imview(x)
    fig = figure("pyplot_surfaceplot",figsize=(10,10))
    imshow(x, cmap=ColorMap("coolwarm"), interpolation="none");
    tight_layout()
end
function imsurf(z)
    fig = figure("pyplot_surfaceplot",figsize=(10,10))
    x = collect(1:size(z,1))
    y = collect(1:size(z,2))
    zz = copy(z)
    zz[isnan.(z)]=0
    surf(x,y,z, facecolors=get_cmap("coolwarm")(zz/maximum(zz)), linewidth=0.25, rstride=1, edgecolors="k", cstride=1,antialiased=true, shade=false)
    tight_layout()
end
