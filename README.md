# High Angular Resolution Imaging ASTR8800

## How to use this repository

First clone the repository
```
git clone https://github.com/fabienbaron/opticscourse.git
```

To update it during the semester: ```git pull```

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
