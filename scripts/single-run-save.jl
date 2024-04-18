using DrWatson, Glob, Revise, LaTeXStrings
@quickactivate
foreach(includet, glob("*.jl", srcdir()))
# includet()

using ProgressMeter, Suppressor, DataFrames
using Plots
using ArgCheck

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false)

DrWatson._wsave(s::String, plot::Plots.Plot) = savefig(plot, s)
DrWatson.default_allowed(::Dict) = (Real, String, Vector, Dict)

allparams = Dict{Symbol,Any}(
    :scaled => false,
    :S => [8,20,50],
    # :μ => Derived(:S, x -> 1 / x),
    :μ => 0.1,
    :C => 1.0,
    # :σ => Derived(:S, x -> 0.1 / sqrt(x)),
    :σ => 0.01,
    :k => 0.75,
    :b0 => 1.0,
    :K => 1e6,
    :λ => 0,
    :z => 0,
    :r => 1,
    :N => 1,
    :threshold => false,
    :dist => "normal",
    :symm => false,
    :seed => 17,
)

# allparams = Dict{Symbol,Any}(
#     :scaled => false,
#     :S => [8,20,50],
#     :μ => 0.01,
#     # :C => 1.0,
#     :σ => 0.005,
#     :k => 1.0,
#     :n0 => 1e-8,
#     :b0 => 1,
#     :K => 8,
#     :λ => 0,
#     :z => 0,
#     :r => 1,
#     :N => 1,
#     :threshold => false,
#     :dist => "normal",
#     :symm => false,
#     :seed => 17,
# )

function logistic_stability_threshold(r, K, μ, σ)
    return ((r / K) - μ) / σ
end


dicts = dict_list(allparams)
pl = Plots.plot()
for (i, d) in enumerate(dicts)
    # @argcheck d[:S] < logistic_stability_threshold(d[:r], d[:K], d[:μ], d[:σ])
    println("Running:\n$d")
    evolve!(d; trajectory=true)
    println(d[:richness])
    global pl = boundary(d, overprint=true)
end

plot!(xlims=(-2,0))
display(pl)

# safesave(plotsdir(savename(allparams, "png")), pl)


