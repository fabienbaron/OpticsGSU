# High Angular Resolution Imaging ASTR8800

## Prerequisites/Installation

OpticsGSU uses OptimPackNextGen, which requires Eric Thiebault's repository as well as the default repo.
At the moment we are using PythonPlot with the Qt5 backend.

```julia
ENV["PYTHON"]=""
using Pkg
pkg"registry add General"
pkg"registry add https://github.com/emmt/EmmtRegistry"
pkg"add https://github.com/fabienbaron/OpticsGSU.git"
```

Once installed, you can update the package by pressing ```]``` to access the package manager, then update with ```up```.


