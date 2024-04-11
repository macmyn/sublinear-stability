using DrWatson, Glob, Revise
@quickactivate
foreach(includet, glob("*.jl", srcdir()))
# includet()

using ProgressMeter, Suppressor, DataFrames
using Plots, LaTeXStrings
using ArgCheck

DrWatson._wsave(s::String, plot::Plots.Plot) = savefig(plot, s)
DrWatson.default_allowed(::Dict) = (Real, String, Vector, Dict)

# allparams = Dict{Symbol,Any}(
#     :scaled => false,
#     :S => 1000,
#     :μ => Derived(:S, x -> 1 / x),
#     :C => 1.0,
#     :σ => Derived(:S, x -> 0.1 / sqrt(x)),
#     :k => 0.75,
#     :b0 => 1.0,
#     :K => 1e5,
#     :λ => 0,
#     :z => 0,
#     :r => 1,
#     :N => 10,
#     :threshold => false,
#     :dist => "normal",
#     :symm => false,
#     :seed => 17,
# )

allparams = Dict{Symbol,Any}(
    :scaled => false,
    :S => [6],
    :μ => 0.01,
    # :C => 1.0,
    :σ => 0.005,
    :k => 1.0,
    :n0 => 1e-8,
    :b0 => 1,
    :K => 20,
    :λ => 0,
    :z => 0,
    :r => 1,
    :N => 1,
    :threshold => false,
    :dist => "normal",
    :symm => false,
    :seed => 19,
)

function logistic_stability_threshold(r, K, μ, σ)
    return ((r / K) - μ) / σ
end


dicts = dict_list(allparams)
pl = Plots.plot()
for (i, d) in enumerate(dicts)
    @argcheck d[:S] < logistic_stability_threshold(d[:r], d[:K], d[:μ], d[:σ])
    println("Running:\n$d")
    evolve!(d; trajectory=true)
    println(d[:richness])
    global pl = boundary(d, overprint=true)
end



println(typeof(pl))
display(pl)

safesave(plotsdir(savename(allparams, "png")), pl)


