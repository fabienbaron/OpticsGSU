# High Angular Resolution Imaging ASTR8800

## How to use this repository

First clone the repository
```
git clone https://github.com/fabienbaron/opticscourse.git
```

To update it during the semester: ```git pull```

## Julia setup

Install Julia from https://julialang.org/downloads/

Install Pyplot by doing ```Pkg.add("PyPlot")``` on the Julia command line. Mac users will probably need to have XQuartz installed before. See here for more details: https://github.com/JuliaPy/PyPlot.jl

Install the Atom editor https://atom.io/.

Install the language-julia Atom package from within Atom (Linux: CTRL+comma then Install tab, Mac:COMMAND+comma then Install tab).

More advanced users can directly use Juno (a Julia + Atom remix) http://junolab.org/

## Check that you can use the REPL

1. Launch julia from the command line.

2. Check that it can do basic operations, e.g. input ```2+2``` then press enter.

3. Check that it can plot basic functions using the PyPlot library.

```julia
using PyPlot
x=linspace(-pi, pi, 100)
plot(x,sin.(x))
```

## Learning Julia

These two links should get you started:
1. https://learnxinyminutes.com/docs/julia/
2. https://math.mit.edu/~stevenj/Julia-cheatsheet.pdf (ignore MIT specifics)

If you already master Matlab or Python also check https://cheatsheets.quantecon.org/

Also see http://juliaastro.github.io for useful functions.
