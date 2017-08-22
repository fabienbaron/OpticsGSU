
# Generating a PSF

To use the provided package:
```julia
include("zernikes.jl")
```

Array pointwise product:
```julia
z=x.*y
```

Sum function:
```julia
z=sum(x)
```

Squared modulii of vector values:
```julia
z=abs2.(x)
```

Shifting arrays:
```julia
circshift()
```

Doing the fft and inverse fft:
```julia
fft()
ifft()
```
