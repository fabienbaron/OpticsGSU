
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


function stirling(n) # improved stirling formula
if(n==0)
     return 1
 else
     return (sqrt(2*pi*n)*(n/e)^n*exp(1/(12*n)-1/(360*n^3)))
end
end

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
            if((n>20)|(m>20))
                output += (-1)^k * stirling(n - k) / (stirling(k) * stirling(div(n + m,2) - k) * stirling(div(n - m, 2) - k)) * rho.^(n - 2 * k)
            else
                output += (-1)^k * factorial(n - k) / (factorial(k) * factorial(div(n + m,2) - k) * factorial(div(n - m, 2) - k)) * rho.^(n - 2 * k)
            end
        end
        return output
    end
end

function zernike(j, npix=128, diameter=128, cent_x=64.5, cent_y=64.5; outside=0, noll_normalize=true, centered=false)
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
if centered == true
  cent_x = (npix+1)/2
  cent_y = (npix+1)/2
end

n, m = noll_indices(j);
x = (collect(linspace(1,npix,npix)) - cent_x) / (diameter / 2.)
y = (collect(linspace(1,npix,npix)) - cent_y) / (diameter / 2.)
xx = repmat(x,1,npix)
yy = repmat(y,1,npix)'
rho = sqrt.(xx.^2 + yy.^2)
theta = atan2.(yy, xx)
aperture = ones(size(rho))
aperture[find(rho.>1)] = outside  # this is the aperture mask
norm_coeff=1.0/(diameter/2)
if (noll_normalize==true)
    norm_coeff = (2*(n+1)/(1+(m==0)))^0.5/(diameter/2)
end
if (m > 0)
    return norm_coeff*zernike_rad(n, m, rho).*cos.(m * theta).*aperture
end
if (m < 0)
    return norm_coeff*zernike_rad(n, m, rho).*sin.(-m * theta).*aperture
end

return norm_coeff*zernike_rad(n, 0, rho).*aperture
end

function circular_aperture(npix=128, diameter=128, cent_x=64.5, cent_y=64.5; outside=0, centered = false)
  """
  Returns a 2D aperture of the desired diameter pixels, centered on (cent_x,cent_y) and on support npix X npix
  """
if centered == true
  cent_x = (npix+1)/2
  cent_y = (npix+1)/2
end

x = (collect(linspace(1,npix,npix)) - cent_x) / (diameter / 2.)
y = (collect(linspace(1,npix,npix)) - cent_y) / (diameter / 2.)
xx = repmat(x,1,npix)
yy = repmat(y,1,npix)'
rho = sqrt.(xx.^2 + yy.^2)
aperture = ones(size(rho))
aperture[find(rho.>1)] = outside  # this is the aperture mask
return aperture
end
