# Generating a PSF

To use the provided package:
```julia
include("optics.jl")
```

Array pointwise product:
```julia
z=x.*y
```

Vector of values from 1 to 512
```julia
z=collect(1:512)
```

Maximum of a vector x
```julia
z=maximum(x)
```

Sum function:
```julia
z=sum(x)
```

Complex exponentiation:
```julia
z=cis.(x)=exp.(im*x)
```


Squared modulii of vector values:
```julia
z=abs2.(x)
```

Shifting arrays:
```julia
shifted=circshift(nonshifted, (deltax, deltay))
```

Doing the fft and inverse fft:
```julia
fft()
ifft()
```
