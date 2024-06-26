using OrdinaryDiffEq, PyPlot, LinearAlgebra, Random, Distributions, ForwardDiff, OMEinsum, DrWatson, ColorSchemes, Revise, Infiltrator, Debugger, Roots
include(srcdir("NonlinearStability.jl"))
# pyplot()
# gr()
using DelimitedFiles
# using StatsPlots
# pgfplots()
# global AA = readdlm("a.txt")

MAXTIME = 100
plt.rc("text", usetex=true)
rc("font", family="serif", weight="normal", size="18")
rc("axes", labelsize="10")
rc("xtick", labelsize="10")
rc("ytick", labelsize="10")

# rc("legend", size="10")
# rc("legendfontsize", size="10")
fig, axs = plt.subplots(2,1,figsize=(3,4), height_ratios=[4, 3])
# plots = plot(layout=(2,1),size=(400,290))
FS = 10
SAVE = false
# plot!(xtickfontsize=FS,ytickfontsize=FS,xguidefontsize=FS,yguidefontsize=FS,legendfontsize=FS,subplot=1)
# plot!(xtickfontsize=FS,ytickfontsize=FS,xguidefontsize=FS,yguidefontsize=FS,legendfontsize=FS,subplot=2)
maxes = []
cm = get_cmap(:tab10)
colorrange = (0:9)./10
function main()

    all_params = Dict{Symbol,Any}(
        :alpha => -1,
        :beta => 1.2,
        :z => 1,
        :r => 1.0,
        :N => [10,20,100],
        :μ => 0.2,
        :σ => 0.02,
        :tspan => (0.0, 150.0),
        :init => "const",
        :type => "sublinear"
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
            eigvs = get_eigvs_sublinear(sol, p)
            
            push!(maxes, maximum(real(eigvs)))
            
            ## Plot things ##
            # colors = palette(:tab10, length(all_params[:N]))
            
            
            # label = "\$N = $(p[:N])\$"

        plot_ts = 0:0.1:p[:tspan][2]

        # Time series        
        # lol julia starts at 1 so this doesn't do anything...
        # plot!(sol[2:end](plot_ts),subplot=1,label=nothing,color=colors[i],alpha=0.5,markercolor =colors[i])
        # plot!(plot_ts,sol(plot_ts)[1,:], subplot=1,label=label,color=colors[i], alpha=0.5)
        # plot!(xlabel="Time", ylabel="Abundances", subplot=1)
        axs[1].plot(plot_ts,sol(plot_ts)[1,:],alpha=0.5,color=cm(colorrange[i]),label=p[:N])
        for j in 2:p[:N]
        axs[1].plot(plot_ts,sol(plot_ts)[j,:],alpha=0.5,color=cm(colorrange[i]))
        end
        # end
        # axs[1].plot(sol[2:end](plot_ts),alpha=0.5)#,color=cm(colorrange[i]))
        # axs[1].plot(plot_ts,sol(plot_ts)[1,:],label=label,alpha=0.5)
        axs[1].set_xlabel(L"\textrm{Time}")
        axs[1].set_ylabel(L"\textrm{Abundances}")
        axs[2].set_xlabel(L"\Re(x)")
        axs[2].set_ylabel(L"\Im(x)")

        # Eigs plot
        # scatter!(eigvs,subplot=2,color=colors[i])
        for j in eachindex(eigvs)
        axs[2].scatter(real(eigvs[j]),imag(eigvs[j]),color=cm(colorrange[i]),marker=".",edgecolor="k",linewidth=0.1)
        end
        # Calculate N_* 
        nstar = vit_sublinear_equilibrium(p)
        axs[1].axhline([nstar],color=cm(colorrange[i]),alpha=0.5)

        pred_eig = vit_sublinear_eigs(p)
        axs[2].axvline([pred_eig],color=cm(colorrange[i]),alpha=0.5)
        println(pred_eig)

        
    end
    plt.tight_layout()
    leg = axs[1].legend(fontsize=10)
    plt.setp(leg.get_texts(), family="serif")
    # plot!(dpi=500)
    # plot!(xlim=(-0.5,0.1),subplot=2)
    if SAVE 
        savefig(plotsdir(savename(all_params,"pdf")),bbox_inches="tight")
    end
    # plots
end
# main()

# Debugger.@enter main()
main()