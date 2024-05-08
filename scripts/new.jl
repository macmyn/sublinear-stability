using OrdinaryDiffEq, Plots, LinearAlgebra, Random, Distributions, ForwardDiff, OMEinsum, DrWatson, ColorSchemes, Revise
include(srcdir("NonlinearStability.jl"))
gr()

tspan = (0.0, 1.0)
plots = plot(layout=(2,1))

a = 1
b = 2.0

all_params = Dict{Symbol,Any}(
:z => 5,
:r => 1,
:alpha => a,
:beta => b,
:N => [20,50,100],
:μ => 0.1,
:σ => 0.05,
:nothing => ""  # required so that types work correctly when calling dict_list
)               # (otherwise it tries to turn the array into a Real)

dicts = dict_list(all_params::Dict{Symbol,Any})
for (i, p) in Iterators.reverse(enumerate(dicts))
    
    # Set interation matrix (normal with zeros on diags)
    p[:A] = get_interaction_matrix(p)

    # Initial condition
    x0 = rand(rng, Uniform(1,5),p[:N])

    params = NamedTuple([pair for pair in p])  # ODEProblem only takes NamedTuple 🙄

    # Define problem and get solution
    prob = ODEProblem(new, x0, tspan, params,)
    sol = solve(prob, Tsit5())

    # Jacobian and eigen values
    final_state = sol[end]
    final_state_b = final_state.^(p[:beta]-1)  # N^{β-1} term (@ein doesn't like it in line below)
    @ein J[i,j] := p[:A][i,j]*final_state[i] * final_state_b[j]  # Build Jacobian from N_i, N_j
    J .*= -1 * p[:beta]  # Global factor (except on diags) of -β
    J[diagind(J)] = -sign(p[:alpha]) * p[:r] * p[:alpha] .* final_state.^p[:alpha]  # Diag terms
    eigvs = eigen(J).values
    
    ## Plot things ##
    label = "\$N = $(p[:N])\$"

    # Time series        
    plot!(sol[1:end],subplot=1,label=nothing,color=colors[i],alpha=0.5)
    plot!(sol[1],subplot=1,label=label,color=colors[i],alpha=0.5)  # Plot with label

    # Eigs plot
    scatter!(eigvs,subplot=2,label=label,color=colors[i])

end

plots
