using DrWatson, Glob, Revise, LaTeXStrings
@quickactivate
foreach(includet, glob("*.jl",srcdir()))

using Plots

plot_font = "Computer Modern"

default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false)

DrWatson._wsave(s::String, plot::Plots.Plot) = savefig(plot, s)
DrWatson.default_allowed(::Dict) = (Real, String, Vector, Dict)

gr(fmt=:png)

## Beta ##
allparams = Dict{Symbol,Any}(
    :scaled => false,
    :S => [10,20,30],
    :μ => 0.1,
    # :C => 1.0,
    :σ => 0.01,
    :k => 1.0,
    :n0 => 1e-8,
    :b0 => 1,
    :K => 1e6,
    :λ => 0,
    :z => 0,
    :r => 1,
    :N => 1,
    # :betaa => 2,
    # :betab => 2.5,
    :threshold => false,
    :dist => "normal",
    :symm => false,
    :seed => 17,
)

betasum = 4
# betaa = 2
betabs = 2.05:0.025:3.5
anim = @animate for i in eachindex(betabs)
# for i in eachindex(betabs)
    plot()
    # plot!(size=(800,300))
    bb = betabs[i]
    aa = betasum - betabs[i]
    allparams[:betaa] = aa
    allparams[:betab] = bb

    dicts = dict_list(allparams)
    # plot()
    # scatter(randn(100))
    for (i,d) in enumerate(dicts)
        evolve!(d; trajectory=true)
        global pl = boundary(d,overprint=true)
    end
    inset(allparams)
    plot!(xlim=(-5,-1),ylim=(-0.05,0.1),subplot=1)
    # frame(anim)
    # display(pl)
end
# display(pl)

# allparams[:betaa] = 2


# plot(pl,i)

gif(anim, "testanim_constsum.gif",fps=3)


# @userplot RangePlot
# @recipe function f(rp::RangePlot)
    
#     p,betaa,betab = RangePlot.args
#     p[:betaa] = betaa
#     p[:betab] = betab

#     dicts = dict_list(p)
#     for (i,d) in enumerate(dicts)
#         evolve!(d; trajectory=true)
#         global pl = boundary(d, overprint=true)
#     end
#     inset(dicts)    

# end


# dicts = dict_list(allparams)
# pl = Plots.plot()
# for (i, d) in enumerate(dicts)
#     # @argcheck d[:S] < logistic_stability_threshold(d[:r], d[:K], d[:μ], d[:σ])
#     # println("Running:\n$d")
#     evolve!(d; trajectory=true)
#     # println(d[:richness])
#     global pl = boundary(d, overprint=true)
# end

# inset(allparams)

# # plot!(xlims=(-2,0))
# display(pl)
# """