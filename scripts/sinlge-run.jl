#= This script allows to evolve a community specified by P. It is possible to plot results and compare with predictions from the cavity solution. =#

using DrWatson, Glob, Revise
@quickactivate
foreach(includet, glob("*.jl", srcdir()))

using ProgressMeter, Suppressor, DataFrames
using Plots, LaTeXStrings

#= Define the system =#

(P = Dict{Symbol,Any}(
    :scaled => false,
    :S => 50,
    :μ => 0.1,
    :C => 1.0,
    :σ => 0.01,
    :k => 0.75,
    :b0 => 1.0,
    :K => 1e5,
    :λ => 0,
    :z => 0,
    :r => 1,
    :N => 1,
    :threshold => false,
    :dist => "normal",
    :symm => false,
    :seed => 17,
);

#= Evolve =#
evolve!(P; trajectory=true);

#= Trajectories =#
plot(P[:trajectory],
    ylabel=L"B_i(t)",
    xlabel=L"t",
    #yscale = :log,
    linewidth=3,
    legend=false,
    alpha=0.5,
    grid=false,
    palette=:YlOrBr_9,
    c=9
))

#= Density distribution =#
histogram([(P[:equilibrium])],
    normalize=true,
    alpha=0.5,
    color=:gray,
    ylabel=L"P(B^*)",
    xlabel=L"B^*",
    label=false,
    grid=false,
)

#= cavity solution =#
ϕ, e1, e2 = Cavity(P, n_max=100)
X = [n for n in 0.95*minimum(P[:equilibrium]):0.0001:1.05*maximum(P[:equilibrium])]
plot!(X, [P_n(n, e1, e2, P) for n in X],
    labels="cavity",
    alpha=1,
    linewidth=2,
    linecolor=:black,
)

#= cavity solution with gaussian approximation =#
X = [n for n in 0.5*minimum(P[:equilibrium]):0.0001:1.5*maximum(P[:equilibrium])]
plot!(X, [pdf(last(P_n_gauss(P)), n) for n in X],
    labels="cavity - gaussian",
    linewidth=2,
    alpha=1,
    linecolor=:black,
    linestyle=:dash,
    legend=:bottomleft
)

#= eigenvalues spectrum =#
boundary(P)