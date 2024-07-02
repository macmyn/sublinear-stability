using Distributions, LinearAlgebra, LaTeXStrings, PyPlot, PyCall, ColorSchemes, DrWatson
import Plots.palette

## USED IN PAPER for May div/stab scaling (randomly samples matrix elements)

plt.rc("text", usetex=true)
PyCall.PyDict(plt."rcParams")["text.latex.preamble"] = "\\usepackage{amsmath}"# \\usepackage{underset}"

means = []
stddevs = []
const SIG = 3.0
dist = Normal(0,SIG^2)

Ns = unique(trunc.(Int,10 .^ range(0, stop=2.4, length=10)))
# cur_colors = palette(:auto)
# cmap = get_cmap(:auto)
# Ns = trunc.(Ns)
repeats = 5
for N in Ns
    lambas_repeat = []
    for i in 1:repeats
    r = rand(dist,(N,N))

    A = r - I

    eigs = eigen(A).values

    l_max = maximum(real.(eigs))
    push!(lambas_repeat,l_max)
    
end

    lambda_avg = mean(lambas_repeat)
    lambda_std = std(lambas_repeat)
    push!(means, lambda_avg)
    push!(stddevs, lambda_std)

end

fig = figure(figsize=(3,1.5))
ax = PyPlot.axes()
ax.spines["top"].set_visible(false)
ax.spines["right"].set_visible(false)
# plot(figsize=(3,3))
ax.plot(0:maximum(sqrt.(Ns)),(0:maximum(sqrt.(Ns)))*SIG^2,linestyle="--",color="grey")
ax.errorbar(sqrt.(Ns), means, yerr=stddevs, linestyle="", marker=".", ecolor="k",capsize=0)
ax.set_xlabel(L"\sqrt{N}")
ax.set_ylabel(L"\max\limits_i \Re (\lambda_i)")


savefig(plotsdir("div_stab_may.png"),bbox_inches="tight",dpi=500)
display(fig)