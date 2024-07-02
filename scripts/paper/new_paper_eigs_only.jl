using OrdinaryDiffEq, PyPlot, LinearAlgebra, Random, Distributions, ForwardDiff, OMEinsum, DrWatson, ColorSchemes, Revise, Infiltrator, Debugger, Roots, PyCall
include(srcdir("NonlinearStability.jl"))
# pyplot()
# gr()
using DelimitedFiles
# using StatsPlots
# pgfplots()
# global AA = readdlm("a.txt")

@pyimport matplotlib.animation as anim

MAXTIME = 100
plt.rc("text", usetex=true)
rc("font", family="serif", weight="normal", size="18")
rc("axes", labelsize="10")
rc("xtick", labelsize="10")
rc("ytick", labelsize="10")

# rc("legend", size="10")
# rc("legendfontsize", size="10")
fig, axs = plt.subplots(1,1,figsize=(4,1.5))
# plots = plot(layout=(2,1),size=(400,290))
FS = 10
SAVE = false
# plot!(xtickfontsize=FS,ytickfontsize=FS,xguidefontsize=FS,yguidefontsize=FS,legendfontsize=FS,subplot=1)
# plot!(xtickfontsize=FS,ytickfontsize=FS,xguidefontsize=FS,yguidefontsize=FS,legendfontsize=FS,subplot=2)
maxes = []
cm = get_cmap(:tab10)
colorrange = (0:9)./10
alphas = [-1,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,1,1.2,1.4,1.6,1.8,2.0]
# alphas = [-1,1.6,1.8,2.0]
function make_frame(k)
    if k >= length(alphas)
        nothing
    else
    all_params = Dict{Symbol,Any}(
        :alpha => alphas[k+1],
        :beta => 1,
        :z => 1,
        :r => 1,
        :N => [10,20,50],
        :μ => 0.2,
        :σ => 0.01,
        :tspan => (0.0, 150.0),
        :init => "const",
        :type => "sublinear_to_may"
        )

        axs.clear()

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
            eigvs = get_eigvs_sublinear(sol, p)[2:end]
            # println(eigvs)
            
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
        # axs[1].plot(plot_ts,sol(plot_ts)[1,:],alpha=0.5,color=cm(colorrange[i]),label=p[:N])
        # for j in 2:p[:N]
        # axs[1].plot(plot_ts,sol(plot_ts)[j,:],alpha=0.5,color=cm(colorrange[i]))
        # end
        # end
        # axs[1].plot(sol[2:end](plot_ts),alpha=0.5)#,color=cm(colorrange[i]))
        # axs[1].plot(plot_ts,sol(plot_ts)[1,:],label=label,alpha=0.5)
        # axs[1].set_xlabel(L"\textrm{Time}")
        # axs[1].set_ylabel(L"\textrm{Abundances}")
        axs.set_xlabel(L"\Re(x)")
        axs.set_ylabel(L"\Im(x)")

        # Eigs plot
        # scatter!(eigvs,subplot=2,color=colors[i])
        # for j in eachindex(eigvs)
        # axs.scatter(real(eigvs[j]),imag(eigvs[j]),color=cm(colorrange[i]),marker=".",edgecolor="k",linewidth=0.1,label=p[:N])
        # end
        axs.scatter(real.(eigvs),imag.(eigvs),color=cm(colorrange[i]),marker=".",edgecolor="k",linewidth=0.1,label="N="*string(p[:N]))
        # Calculate N_* 
        nstar = vit_sublinear_equilibrium(p)
        # axs[1].axhline([nstar],color=cm(colorrange[i]),alpha=0.5)

        pred_eig = vit_sublinear_eigs(p)
        # axs.axvline([pred_eig],color=cm(colorrange[i]),alpha=0.5)
        println(pred_eig)

    end
    end

    plt.tight_layout()
    # leg = axs[1].legend(fontsize=10)
    # plt.setp(leg.get_texts(), family="serif")
    # plot!(dpi=500)
    # plot!(xlim=(-0.5,0.1),subplot=2)
    plt.xlim(-4,0)
    plt.ylim(-0.01,0.01)
    plt.xlabel("Stability")
    # plt.legend(bbox_to_anchor=(1.4, 1.1,2,2),mode="expand",fontsize=15)
    plt.tight_layout()
    if k == 1
        plt.savefig("frame1.png",bbox_inches="tight",dpi=500)
    end
end
    # display(fig)
    # if SAVE 
    #     savefig(plotsdir("sublinear_to_may.pdf"),bbox_inches="tight")
    # end
    # axs.set_xscale("symlog")
    # plt.ylim(-0.005,0.005)
    # plots

# main()
function main()
    myanim = anim.FuncAnimation(fig, make_frame, frames=size(alphas,1)+5,interval=800)
    myanim.save("test.gif",dpi=500)
end

# Debugger.@enter main()
main()