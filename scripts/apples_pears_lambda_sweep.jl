using OrdinaryDiffEq, Plots, LinearAlgebra, Random, Distributions, ForwardDiff, OMEinsum, DrWatson, ColorSchemes, Revise, Infiltrator, Debugger, Roots, FileIO
include(srcdir("NonlinearStability.jl"))
pyplot()
gr()
using DelimitedFiles

# global AA = readdlm("a.txt")

MAXTIME = 100

plots = plot(layout=(2,1))

maxes = []
all_eigvs = []
function main()

    ####################################
    # dt(N) = N [sign(α)[r1 - r*N^α] + sign(β) * Σ[r2 - a * N^β]]
    ####################################


    all_params = Dict{Symbol,Any}(
        :type => "apples_pears",
        :correction => true,
        :alpha => -1,
        :beta => 1,
        :r => 1.0,
        :r1 => 1,
        :r2 => 0.1,
        :N => [10,20,50,200,500,1000,2000,4000],
        :μ => 0.2,
        :σ => 0.01,
        :tspan => (0.0, 100.0),
        :init => "const",
        )

        dicts = dict_list(all_params::Dict{Symbol,Any})
        for (i, p) in tqdm(Iterators.reverse(enumerate(dicts)))
        # for (i, p) in enumerate(dicts)
            
            # Set interation matrix (normal with zeros on diags)
            p[:A] = get_interaction_matrix(p)
            x0 = get_initial_condition(p)
            # p[:z] = calculate_rfix(p)

            params = NamedTuple([pair for pair in p])  # ODEProblem only takes NamedTuple 🙄

            # Define problem and get solution
            prob = ODEProblem(apple_pear_interactions, x0, p[:tspan], params)
            sol = solve(prob, AutoTsit5(Rosenbrock23()))

            # Jacobian and values
            # println("\n\n\nDONE HERE")
            global eigvs = get_eigvs_apples_pears(sol, p)
            eigvs = complex.(eigvs)
            push!(all_eigvs, eigvs)
            # eigvs = get_eigvs_PREVIOUS_WRONG(sol,p)
            
            push!(maxes, maximum(real(eigvs)))
            
            ## Plot things ##
            # colors = palette(:tab10, length(all_params[:N]))
            # label = "\$N = $(p[:N])\$"

            # plot_ts = 0:0.1:p[:tspan][2]

            # Time series        
            # lol julia starts at 1 so this doesn't do anything...
            # plot!(sol[2:end](plot_ts),subplot=1,label=nothing,color=colors[i],alpha=0.2,markercolor =colors[i])
            # plot!(plot_ts,sol(plot_ts)[1,:], subplot=1,label=label,color=colors[i], alpha=0.2)
            # plot!(xlabel="Time", ylabel="Species abundance", subplot=1)

            # Eigs plot
            # scatter!(eigvs,subplot=2,label=label,color=colors[i])
            
            # nstar = apples_pears_equilibrium(p, p[:correction])
            # hline!([nstar], subplot=1, color=colors[i])

            # # Estimated λ_m
            # global lambda_m = -1 * p[:r] * abs(p[:alpha]) * nstar^p[:alpha] + p[:μ]*abs(p[:beta])*nstar^p[:beta] + p[:σ]*nstar^p[:beta]*sqrt(p[:N])
            # vline!([lambda_m], subplot=2, color=colors[i])

            # # Outlier
            # global lambda_o = -1 * p[:r] * abs(p[:alpha]) *nstar^p[:alpha] - abs(p[:beta]) * p[:μ] * nstar^p[:beta] *(p[:N]-1)
            # vline!([lambda_o], subplot=2, color=colors[i], linestyle=:dash)

        end
    # plot!(xlim=(-6,-4),subplot=2)
    # plot!(xlim=(-0.8,0),subplot=2)
    safesave(datadir(savename(all_params,"jld2")),all_eigvs)
    # plot!(ylim=(-1,1),subplot=2)
    # plot!(dpi=500)
    # plots
    # safesave(plotsdir(savename(all_params, "png")),plots)
end
# main()

# Debugger.@enter main()
main()
