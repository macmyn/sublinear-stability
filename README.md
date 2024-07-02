# Sublinear Stability

This repo was originally a clone of the [repo](https://zenodo.org/records/10476102) from "Diversity begets stability: sublinear growth and competitive coexistence across ecosystems", Ian A. Hatton, Onofrio Mazzarisi, Ada Altieri, Matteo Smerlak

### How to use (from scratch)
1. [Install Julia](https://julialang.org/downloads/)
    
    To use Julia, type `julia` in the command line. To run a particular file, type `include("path/to/file.jl")`. The first time round it may be very slow as it has to precompile things. The second time (even if you change params) it will be quick 🚀
    
    To add packages, type `]` in the Julia command line -- this takes you to the `pkg` package manager. Backspace to get out. To actually add a package, type `]`, then `add DrWatson` (for example)

2. [Activate project](https://juliadynamics.github.io/DrWatson.jl/dev/project/#Reproducibility) 
    
    `DrWatson` is a scientific project management package for Julia. The idea is to make things maximally reproducible. The two or three steps here should completely install everything you need (note the `julia>` and `pkg>` prompts on the left; type `]` to get into `pkg` mode). You should do this in the top-level directory of this project (e.g. `user/Documents/sublinear-stability/`)

3. (Optional) Set up Julia in VS Code for maximal efficiency:
    - Install Julia extension for VS code
    - This previews each plot you produce in the Julia Plots interface inside a VS Code window


### Files
- `DrWatson` splits files into `src` and `scripts`. Source code (functions etc) is in `src`; you don't run these directly. Scripts are for running things. 

Here I'm listing the important files that I made. Compare the git history or compare to the paper's original repo to see which files are used in the original paper.
- `src/NonlinearStability.jl`: all the dirty functions etc for simulating, finding Jacobians, etc. Use the line `include("NonlinearStability.jl")` at the top of each script you want to use
- `scripts/simulate_single_community.jl`: run a single community with params (including numbers of species) defined in the `p` dictionary. The functions are to be found in `NonlinearStability.jl`
- `scripts/sweep.jl`: simulates a number of communities for a range of alpha and beta; saves these to `datadir`
- `scripts/plot_d_s.jl`: plots the data above (used for the heatmap in the write-up)
- `scripts/paper`: various plots that are (or aren't) used in the write-up according to which folder they're in. Should all be self-contained (i.e. don't need data)
- `scripts/apples_pears.jl`: simulating with constants in both inter and intra terms
- `scripts/apples_pears(_lambda)_sweep`: create the butterfly/fish-looking plot and the eigenvalue plot for the model above