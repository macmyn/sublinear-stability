using OrdinaryDiffEq, Plots, LinearAlgebra, Random, Distributions, ForwardDiff, OMEinsum, DrWatson, ColorSchemes, Revise, Infiltrator
include(srcdir("NonlinearStability.jl"))
gr()

tspan = (0.0, 1000.0)
plots = plot(layout=(2,1))

a = 2
b = 1.2

all_params = Dict{Symbol,Any}(
    :z => 1,
    # :r => 1,
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
        p[:r] = calculate_rfix(p)
        # p[:r] = 1
        println(p[:r])
        
        # Initial condition
        x0 = fill(0.1,p[:N])
        
        params = NamedTuple([pair for pair in p])  # ODEProblem only takes NamedTuple 🙄
        
        # Define problem and get solution
        prob = ODEProblem(general_interactions, x0, tspan, params,)
        sol = solve(prob, AutoTsit5(Rosenbrock23()))

        # Jacobian and values
        final_state = sol[end]
        final_state_b = final_state.^(p[:beta]-1)  # N^{β-1} term (@ein doesn't like it in line below)
        @ein J[i,j] := p[:A][i,j]*final_state[i] * final_state_b[j]  # Build Jacobian from N_i, N_j
        J .*= -1 * p[:beta]  # Global factor (except on diags) of -β
        J[diagind(J)] = -sign(p[:alpha]) * p[:r] * p[:alpha] .* final_state.^p[:alpha]  # Diag terms
        eigvs = eigen(J).values
        print(maximum(real(eigvs)))
        
        ## Plot things ##
        colors = palette(:tab10, length(all_params[:N]))
        label = "\$N = $(p[:N])\$"
        
        # Time series        
    # lol julia starts at 1 so this doesn't do anything...
    plot!(sol[2:end],subplot=1,label=nothing,color=colors[i],alpha=0.5)
    plot!(sol.t, sol[1,:], subplot=1,label=label,color=colors[i], alpha=0.5)
    # plot!(sol[1],subplot=1,label=label,color=colors[i],alpha=0.5)  # Plot with label
    plot!(yticks=[0,1,2,3,4,5],subplot=1)

    # Eigs plot
    scatter!(eigvs,subplot=2,label=label,color=colors[i])

end

plots
