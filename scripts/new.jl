using OrdinaryDiffEq, Plots, LinearAlgebra, Random, Distributions, ForwardDiff, OMEinsum, DrWatson, ColorSchemes, Revise, Infiltrator, Debugger, Roots
include(srcdir("NonlinearStability.jl"))
gr()
using DelimitedFiles

# global AA = readdlm("a.txt")

MAXTIME = 100
macro timeout(seconds, expr, fail)
    quote
        tsk = @task $expr
        schedule(tsk)
        Timer($seconds) do timer
            istaskdone(tsk) || Base.throwto(tsk, InterruptException())
        end
        try
            fetch(tsk)
        catch _
            $fail
        end
    end
end

plots = plot(layout=(2,1))

maxes = []

function main()

    all_params = Dict{Symbol,Any}(
        :alpha => 1,
        :beta => 0.5,
        :z => 1,
        :r => 1.0,
        :N => [20],
        :μ => 0.1,
        :σ => 0.02,
        :tspan => (0.0, 10.0),
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
            # sol = solve(prob, AutoTsit5(Rosenbrock23()))
            # sol = @timeout MAXTIME begin
            #     sol = solve(prob, AutoTsit5(Rosenbrock23()))
            # println(sol)
            # end NaN
            sol = @timeout MAXTIME begin
                sol = solve(prob, Tsit5())
            end NaN
            println(sol)
            # println("solved. eiging...")
            # Jacobian and values
            println("\n\n\nDONE HERE")
            eigvs = get_eigvs(sol, p)
            # println(eigvs)
            
            push!(maxes, maximum(real(eigvs)))
            
            ## Plot things ##
            colors = palette(:tab10, length(all_params[:N]))
            label = "\$N = $(p[:N])\$"

        plot_ts = 0:0.1:10

            # Time series        
        # lol julia starts at 1 so this doesn't do anything...
        plot!(sol[2:end](plot_ts),subplot=1,label=nothing,color=colors[i],alpha=0.5)
        plot!(plot_ts,sol(plot_ts)[1,:], subplot=1,label=label,color=colors[i], alpha=0.5)
        # plot!(sol[1],subplot=1,label=label,color=colors[i],alpha=0.5)  # Plot with label
        # plot!(ylims=(0,0.3),subplot=1)

        # Eigs plot
        scatter!(eigvs,subplot=2,label=label,color=colors[i])

    end
    # println(maxes)
    # if maxes[end] < maxes[1]  # remember we're reversing
    #     println("Div does not ⟹ stab")
    # else
    #     println("Div ⟹ stab")
    # end

    plots
end
# main()

# Debugger.@enter main()
main()