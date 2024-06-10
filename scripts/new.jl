using OrdinaryDiffEq, Plots, LinearAlgebra, Random, Distributions, ForwardDiff, OMEinsum, DrWatson, ColorSchemes, Revise, Infiltrator, Debugger, Roots
include(srcdir("NonlinearStability.jl"))
pyplot()
gr()
using DelimitedFiles

# global AA = readdlm("a.txt")

MAXTIME = 100

plots = plot(layout=(2,1))

maxes = []

function main()

    all_params = Dict{Symbol,Any}(
        :alpha => -1,
        :beta => 1,
        :z => 1,
        :r => 1.0,
        :N => [20,50,100],
        :μ => 0.2,
        :σ => 0.02,
        :tspan => (0.0, 100.0),
        :init => "const",
        )

        dicts = dict_list(all_params::Dict{Symbol,Any})
        for (i, p) in Iterators.reverse(enumerate(dicts))
            
            # Set interation matrix (normal with zeros on diags)
            p[:A] = get_interaction_matrix(p)
            x0 = get_initial_condition(p)
            # p[:z] = calculate_rfix(p)

            params = NamedTuple([pair for pair in p])  # ODEProblem only takes NamedTuple 🙄

            # Define problem and get solution
            prob = ODEProblem(general_interactions, x0, p[:tspan], params,)
            sol = solve(prob, AutoTsit5(Rosenbrock23()))
            # sol = @timeout MAXTIME begin
            #     sol = solve(prob, AutoTsit5(Rosenbrock23()))
            # println(sol)
            # end NaN
            # sol = @timeout MAXTIME begin
            #     sol = solve(prob, Tsit5())
            # end NaN
            # Jacobian and values
            println("\n\n\nDONE HERE")
            eigvs = get_eigvs(sol, p)
            
            push!(maxes, maximum(real(eigvs)))
            
            ## Plot things ##
            colors = palette(:tab10, length(all_params[:N]))
            label = "\$N = $(p[:N])\$"

        plot_ts = 0:0.1:p[:tspan][2]

        # Time series        
        # lol julia starts at 1 so this doesn't do anything...
        plot!(sol[2:end](plot_ts),subplot=1,label=nothing,color=colors[i],alpha=0.5,markercolor =colors[i])
        plot!(plot_ts,sol(plot_ts)[1,:], subplot=1,label=label,color=colors[i], alpha=0.5)
        plot!(xlabel="Time", ylabel="Species abundance", subplot=1)

        # Eigs plot
        scatter!(eigvs,subplot=2,label=label,color=colors[i])

    end
    plot!(dpi=500)
    plots
    safesave(plotsdir(savename(all_params, "png")),plots)
end
# main()

# Debugger.@enter main()
main()