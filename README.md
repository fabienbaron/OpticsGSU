# High Angular Resolution Imaging ASTR8800

## Julia setup

Install Julia from https://julialang.org/downloads/
Install Pyplot, see here: https://github.com/JuliaPy/PyPlot.jl
Install the Atom editor https://atom.io/.
Install the language-julia Atom package from within Atom.

## Check that you can use the REPL
Launch julia
Check that it can do basic operations, e.g. input ```2+2``` then press enter.
Check that it can plot basic functions using pyplot

```julia
using PyPlot
x=linspace(-pi, pi, 100)
plot(x,sin.(x))
```

## Learning Julia

These two links should get you started
https://learnxinyminutes.com/docs/julia/
https://math.mit.edu/~stevenj/Julia-cheatsheet.pdf

If you already master Matlab or Python also check https://cheatsheets.quantecon.org/
For astronomy students, also see http://juliaastro.github.io for useful functions.
