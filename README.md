# High Angular Resolution Imaging ASTR8800

## Prerequisites/Installation

OpticsGSU uses OptimPackNextGen, which requires Eric Thiebault's repository as well as the default repo.
At the moment we are using PyPlot with the Qt5 backend but will move to Makie soon.

```julia
ENV["PYTHON"]=""
ENV["MPLBACKEND"]="qt5agg"
using Pkg
pkg"registry add General"
pkg"registry add https://github.com/emmt/EmmtRegistry"
pkg"add https://github.com/fabienbaron/OpticsGSU.git"
```

Only proceed with this installation once you've got all the prerequisites installed.
From Julia's command line, press ```]``` to access the package manager.
Add the package for this course: ```add https://github.com/fabienbaron/OpticsGSU```

