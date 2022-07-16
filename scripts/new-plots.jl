using DrWatson, Glob, Revise
@quickactivate
foreach(includet, glob("*.jl", srcdir()))

using Plots; gr(dpi = 500)
using LaTeXStrings

P(k, S, K, μ) = Dict{Symbol, Any}(
    :scaled => false,
    :threshold => 1,
    :b0 => 0.01,
    :z => .1,
    :S => S,
    :μ => μ,
    :σ => μ,
    :k => k,
    :λ => 0,
    :K => k == 1 ? K : Inf,
    :dist => "gamma",
    :N => 100,
    :seed => rand(UInt),
    :symm => false
);


# vary S
S_range = 5:5:100
z = .1
for μ ∈ (.01, .1, 1.)

    plot(
        xlabel = "Number of species " * L"S",
        ylabel = "Fraction of surviving species",
        title = "r = 1, z = $z, μ = σ = $μ (gamma)"
    )

    k = 1
    for K ∈ (.1, 100) 
        plot!(
            S_range,
            survivors.(P.(k, S_range, K, μ)) ./ S_range,
            label = "k = $k, K = $K",
            markers = :auto,
            color = k == 1 ? COLOR_LOG : COLOR_SUB,
        )
    end 

    k = .75
    plot!(
        S_range,
        survivors.(P.(k, S_range, K, μ)) ./ S_range,
        label = "k = $k",
        color = k == 1 ? COLOR_LOG : COLOR_SUB,
    )
    savefig(plotsdir("survivors_μ=σ=$μ.png"))
end

# vary μ = σ
S = 50
μ_range = .01:.05:1

plot(
    xlabel = "Interaction strength (μ = σ)",
    ylabel = "Fraction of surviving species",
    title = "r = 1, z = $z, S = 50"
)
k = 1
for K ∈ (.1, 100) 
    plot!(
        μ_range,
        survivors.(P.(k, S, K, μ_range)) ./ S,
        label = "k = $k, K = $K",
        markers = :auto,
        color = k == 1 ? COLOR_LOG : COLOR_SUB,
    )
end 

k = .75
plot!(
    μ_range,
    survivors.(P.(k, S, K, μ_range)) ./ S,
    label = "k = $k",
    color = k == 1 ? COLOR_LOG : COLOR_SUB,
)
savefig(plotsdir("survivors_S=$S.png"))


# vary μ = σ
S = 50
μ_range = .01:.1:1

plot(
    xlabel = "Interaction strength (μ = σ)",
    ylabel = "Total biomass",
    title = "r = 1, z = $z, S = 50"
)
k = 1
for K ∈ (.1, 100) 
    plot!(
        μ_range,
        biomass.(P.(k, S, K, μ_range)),
        label = "k = $k, K = $K",
        markers = :auto,
        color = k == 1 ? COLOR_LOG : COLOR_SUB,
    )
end 

k = .75
plot!(
    μ_range,
    biomass.(P.(k, S, K, μ_range)),
    label = "k = $k",
    yaxis = :log,
    color = k == 1 ? COLOR_LOG : COLOR_SUB,
)
savefig(plotsdir("biomass_S=$S.png"))



### trajectory

pp = P(1, 80, 100, 1.)
pp[:N] = 10
survivors(pp)

evolve!(pp; trajectory = true)
plot(pp[:trajectory], legend = false)

pp[:trajectory][end]
