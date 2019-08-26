# High Angular Resolution Imaging ASTR8800

## Installation

From Julia's command line, press ```]``` to access the package manager.
Add the package for this course: ```add https://github.com/fabienbaron/OpticsGSU```
Exit the package manager, and do: ```using OpticsGSU```.

## Updates
This is a reminder that you can use the Julia's package manager to update this package: ```update OpticsGSU```.

## Prerequisites

## Julia setup

Install the latest Julia from https://julialang.org/downloads/

Install Pyplot library from the Julia command line. Mac users will probably need to have XQuartz installed before. See here for more details: https://github.com/JuliaPy/PyPlot.jl

Install the Atom editor https://atom.io/.

Install the language-julia Atom package from within Atom: this enables syntax highlighting of Julia code.

Install the latex-completions Atom package from within Atom: this enables the autocompletion of LaTeX symbols (α, β, ...)

## Check that you can use the REPL

1. Launch julia from the command line.

2. Check that it can do basic operations, e.g. input ```2+2``` then press enter.

3. Check that it can plot basic functions using the PyPlot library.

```julia
using PyPlot
x=range(-pi, pi, length=100)
plot(x,sin.(x))
grid()
```

## Learning Julia

These two links should get you started:
1. https://learnxinyminutes.com/docs/julia/
2. https://www.sas.upenn.edu/~jesusfv/Chapter_HPC_8_Julia.pdf

If you already master Matlab or Python also check https://cheatsheets.quantecon.org/

Also see http://juliaastro.github.io for useful functions.
