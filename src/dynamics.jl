#= Dynamics and auxiliary functions are defined =#
# module Dynamics
using DifferentialEquations
using Random, Distributions
using LinearAlgebra
using Distributions
using ForwardDiff

# function production(x, p)
#     return (x > p[:b0] * p[:threshold] ? x^p[:k] : 0)
#     # if x > B0: x^k; else 0
# end
# function dproduction(x, p)
#     return (x > p[:b0] * p[:threshold] ? p[:k] * x^(p[:k] - 1) : 0)
# end

# function F!(f, x, p)
#     f .= p[:r] * p[:b0]^(1 - p[:k]) .* (production.(ppart.(x), Ref(p)) .- x .^ 2 ./ p[:K]) .- p[:z] * x .- x .* (p[:a] * x) .+ p[:λ]
# end      #      r*B0^(1-k)                        B^k if x>B0 else 0    - B^2 / K              - z * B       - x * A * B     + lambda (always 0?)
#          # logis: 1-k=0 --> r * 1             logis: k=1 --> B

# function betapdf(p,x::Real)
#     a = p[:betaa]
#     b = p[:betab]
#     beta_fn_value = Beta(a,b)
#     return ( x > p[:b0]*p[:threshold] ? x*pdf(beta_fn_value,x) : 0)
#     # return pdf(Beta(a,b),x)
# end

# function dbetapdfdx(p,x::Real)
#     return ForwardDiff.derivative(x->betapdf(p,x),x)
# end

# function dbetapdfdx(p,x::Real)
#     a = p[:betaa]
#     b = p[:betab]
#     dd = @. (1.0/beta(a,b)) * x^(a-2) * (1-x)^(b-2) * (a - x*(a+b-2) -1)
#     return ( x > p[:b0] * p[:threshold] ? dd : 0)
# end

# function Fbeta!(f,x,p)
#     f .= p[:r] * betapdf.(Ref(p),x) .- p[:z]*x .- x.*(p[:a]*x) .+ p[:λ]
# end

function vitdenom(p)
    return p[:b0] + p[:d]*(p[:z]/p[:r])^(1/(1-p[:k]))
end


function vit(p,x::Real)
    r = p[:r]
    z = p[:z]
    d = p[:d]
    k = p[:k]

    denom = vitdenom(p)
    return r * ((x+d)/denom)^(k-1) - z
end

function dvitdx(p,x::Real)
    return ForwardDiff.derivative(x->vit(p,x),x)
end

function FVit!(f,x,p)
    f.= x.*vit.(Ref(p),x) .- x.*(p[:a]*x) .+ p[:λ]
end

# function J!(j, x, p)
#     j = -x .* p[:a]
#     j[diagind(j)] .= p[:r] * p[:b0]^(1 - p[:k]) .* (dproduction.(x, Ref(p)) .- 2 * x ./ p[:K]) .- p[:z] .- p[:a] * x

#     # above_threshold = x .> p[:b0]
#     # j[diagind(j)[above_threshold]] .+= p[:k].*p[:r][above_threshold].*x[above_threshold].^(p[:k]-1)
# end

# function J(x, p)
#     j = -x .* p[:a]
#     j[diagind(j)] .= p[:r] * p[:b0]^(1 - p[:k]) .* (dproduction.(x, Ref(p)) .- 2 * x ./ p[:K]) .- p[:z] .- p[:a] * x
#     return j
# end

# function Jbeta(x, p)
#     j = -x .* p[:a]
#     j[diagind(j)] .= p[:r] * dbetapdfdx.(Ref(p),x)
#     return j
# end

function JVit(x,p)
    r = p[:r]
    z = p[:z]
    d = p[:d]
    k = p[:k]
        
    j = -x .* p[:a]
    
    denom = vitdenom(p)
    j[diagind(j)] .= -z .+ r.*((x.+d)/denom).^(k-1) + (k-1)*r.*x.*(((x.+d)/denom).^(k-2))/denom .- p[:a] * x
    return j

end


## solving

MAX_TIME = 1e3
MAX_ABUNDANCE = 1e3
TOL = 1e-3

converged(ϵ=TOL) = TerminateSteadyState(ϵ)
blowup(max_abundance=MAX_ABUNDANCE) = DiscreteCallback((u, t, integrator) -> maximum(u) > max_abundance, terminate!)

function evolve!(p; trajectory=false)

    if !haskey(p, :rng)
        p[:rng] = MersenneTwister(p[:seed])
    end
    if !haskey(p, :a)
        add_interactions!(p)
    end
    if !haskey(p, :r)
        add_growth_rates!(p)
    end
    if !haskey(p, :x0)
        add_initial_condition!(p)
    end

    pb = ODEProblem(
        ODEFunction(
            (f, x, p, t) -> FVit!(f, x, p); #in-place F faster
            # jac=(j, x, p, t) -> J!(j, x, p) #specify jacobian speeds things up
        ),
        p[:x0], #initial condition
        (0.0, MAX_TIME),
        p
    )

    sol = solve(pb,
        callback=CallbackSet(converged(), blowup()),
        save_on=trajectory #don't save whole trajectory, only endpoint
    )

    p[:equilibrium] = sol.retcode == SciMLBase.ReturnCode.Terminated ? sol.u[end] : NaN
    p[:converged] = (sol.retcode == SciMLBase.ReturnCode.Terminated && maximum(p[:equilibrium]) < MAX_ABUNDANCE)
    p[:richness] = sum(sol.u[end] .> p[:b0] * p[:threshold])
    p[:diversity] = p[:richness] == 0 ? 0 : Ω(sol.u[end] .* (sol.u[end] .> p[:b0] * p[:threshold]))

    if trajectory
        p[:trajectory] = sol
    end

end

function add_initial_condition!(p)
    p[:x0] = rand(p[:rng], Uniform(2, 10), p[:S])
end

function equilibria!(p)
    if !haskey(p, :rng)
        p[:rng] = MersenneTwister(p[:seed])
    end
    if !haskey(p, :a)
        add_interactions!(p)
    end
    if !haskey(p, :r)
        add_growth_rates!(p)
    end


    equilibria = Vector{Float64}[]
    sizehint!(equilibria, p[:N])

    Threads.@threads for _ in 1:p[:N]
        add_initial_condition!(p)
        pb = ODEProblem(
            ODEFunction(
                (f, x, p, t) -> Fbeta!(f, x, p); #in-place F faster
                # jac=(j, x, p, t) -> J!(j, x, p) #specify jacobian speeds things up
            ),
            p[:x0],
            (0.0, MAX_TIME),
            p
        )
        sol = solve(pb,
            callback=CallbackSet(TerminateSteadyState(1e-3), blowup()),
            save_on=false #don't save whole trajectory, only endpoint
        )
        push!(equilibria, sol.u[end])
    end
    p[:equilibria] = uniquetol(equilibria, atol=0.1)
    p[:num_equilibria] = length(p[:equilibria])
    p[:num_interior_equilibria] = sum(map(x -> all(x .> p[:b0] * p[:threshold]), p[:equilibria]))

    return p[:num_equilibria]
end



function add_interactions!(p)
    # add a random interaction matrix to p, the dict of parameters
    (m, s) = p[:scaled] ? (p[:μ] / p[:S], p[:σ] / sqrt(p[:S])) : (p[:μ], p[:σ])

    #choose the distribution
    if p[:dist] == "normal"
        dist = Normal(m, s)
    elseif p[:dist] == "uniform"
        dist = Uniform(max(0.0, m - s), min(2m, m + s))
    elseif p[:dist] == "gamma"
        dist = Gamma(m^2 / s^2, s^2 / m)
    end

    #fill the interaction matrix
    a = rand(p[:rng], dist, (p[:S], p[:S]))

    #implement connectance
    if haskey(p, :C) && p[:C] != 1
        c = rand(Binomial(1, p[:C]), (p[:S], p[:S]))
        println(c)
        a = a .* c
        println(a)
    end

    #self-interactions are implemented above
    a[diagind(a)] .= 0.0

    #implement eventual symmetry
    if p[:symm]
        for i in 1:p[:S], j in 1:p[:S]
            if i > j
                a[i, j] = a[j, i]
            end
        end
    end

    p[:a] = a
end

function add_growth_rates!(p)
    if haskey(p, :dist_r)
        p[:r] = rand(p[:rng], p[:dist_r], p[:S])
    else
        p[:r] = ones(p[:S])
    end
end

function stability!(p)
    # run N simulates and append results to p

    p[:rng] = MersenneTwister(p[:seed])
    stability = Vector{Bool}(undef, p[:N])
    diversity = Vector{Float64}(undef, p[:N])
    richness = Vector{Float64}(undef, p[:N])

    Threads.@threads for i in 1:p[:N]
        add_interactions!(p)
        add_growth_rates!(p)
        evolve!(p)
        stability[i] = p[:converged]
        diversity[i] = p[:diversity]
        richness[i] = p[:richness]
    end
    delete!(p, :converged)
    delete!(p, :rng)

    p[:richness] = mean(richness) / p[:S]
    p[:prob_stab] = mean(stability)


    p[:diversity] = mean(diversity)
    p[:diversity_se] = std(diversity) / sqrt(p[:N])
end


function diversity(p)
    # run N simulates and append results to p
    diversity = Vector{Float64}(undef, p[:N])
    p[:rng] = MersenneTwister()

    for i in 1:p[:N]
        add_interactions!(p)
        add_growth_rates!(p)
        add_initial_condition!(p)
        evolve!(p)
        diversity[i] = p[:diversity]
    end

    return mean(diversity)
end

function full_coexistence(p)
    # run N simulates and append results to p
    full_coexistence = Vector{Float64}(undef, p[:N])
    p[:rng] = MersenneTwister()

    for i in 1:p[:N]
        add_interactions!(p)
        add_growth_rates!(p)
        add_initial_condition!(p)
        evolve!(p)
        full_coexistence[i] = p[:converged] && p[:richness] == p[:S] ? 1 : 0
    end

    return mean(full_coexistence)
end

function ahmadian(p)
    # run N simulates and append results to p
    ahmadian = Vector{Bool}(undef, p[:N])
    p[:rng] = MersenneTwister()
    @unpack μ, σ, k = p
    for i in 1:p[:N]
        add_interactions!(p)
        add_growth_rates!(p)
        add_initial_condition!(p)
        evolve!(p)
        ahmadian[i] = sum(1 ./ (μ .- (1 - k) * ppart.(p[:equilibrium]) .^ (k - 2))) < 1 / σ^2
    end


    return mean(ahmadian)
end
# end