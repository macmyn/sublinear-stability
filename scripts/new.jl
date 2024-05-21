using OrdinaryDiffEq, Plots, LinearAlgebra, Random, Distributions, ForwardDiff, OMEinsum, DrWatson, ColorSchemes, Revise, Infiltrator, Debugger
# include(srcdir("NonlinearStability.jl"))
gr()
using DelimitedFiles

global AA = readdlm("a.txt")
# AA[diagind(AA)] .= 0

plots = plot(layout=(2,1))

maxes = []

function main()

    all_params = Dict{Symbol,Any}(
        :alpha => 1,
        :beta => 1,
        :z => 1,
        # :r => 1.0,
        :N => [100],
        :μ => 0.0,
        :σ => 0.02,
        :tspan => (0.0, 10.0),
        :init => "const",
        )
        
        dicts = dict_list(all_params::Dict{Symbol,Any})
        for (i, p) in Iterators.reverse(enumerate(dicts))
            
            # Set interation matrix (normal with zeros on diags)
            # p[:A] = get_interaction_matrix(p)
            p[:A] = AA
            # x0 = get_initial_condition(p)
            x0 = fill(0.3,100)
            # println("found initial")
            # println("INITIAL: $x0")
            # println("typeof: $(typeof(x0))")

            p[:r] = calculate_rfix(p)
            
            # println("got r fix@")
            params = NamedTuple([pair for pair in p])  # ODEProblem only takes NamedTuple 🙄
            # @bp
            # Define problem and get solution
            prob = ODEProblem(general_interactions, x0, p[:tspan], params,)
            sol = solve(prob, AutoTsit5(Rosenbrock23()))
            println(sol[2])
            # println("solved. eiging...")
            # Jacobian and values
            # eigvs = get_eigvs(sol, p)
            # println(eigvs)
            
            # push!(maxes, maximum(real(eigvs)))
            
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
        # scatter!(eigvs,subplot=2,label=label,color=colors[i])

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





function get_interaction_matrix(p)
    d = Normal(p[:μ], p[:σ])
    A = rand(d, p[:N], p[:N])
    A[diagind(A)] .= 0
    return A
end

function calculate_rfix(p)
    r = sum(p[:A], dims=2)
    return r
end

function equilib_condition(N, p)
    return sign(p[:alpha])*p[:r]*N^p[:alpha] + (p[:N]-1)*p[:μ]*N^p[:beta] - sign(p[:alpha])*p[:z]
end

function solve_initial(p)
    fz = find_zero(f, 0.1)
    return fz    
end

function get_initial_condition(p)
    if p[:init] == "uniform"
        x0 = rand(rng, Uniform(1,5),p[:N])
    elseif p[:init] == "const"
        x0 = fill(0.2,p[:N])
    elseif p[:init] == "solve"
        root = solve_initial(p)
        x0 = fill(root, p[:N])
        # println("xo fill is: $x0")
    else 
        throw(ArgumentError("Incorrect input for p[:init]"))
    end
    return x0
end

function general_interactions(db, b, p, t)
    # @bp
    # b[b.<1e-3] .= 0
    # println(b)
    # @infiltrate

    # off_diags = p[:A] * (b.^p[:beta])
    # diags = sign(p[:alpha]) .* p[:r] .* (b.^p[:alpha])
    # z_term = sign(p[:alpha]) * p[:z]
    # bracket = z_term .- off_diags .- diags

    ## Ming ##
    # diags = b.^p[:alpha]
    
    # this apepars to work
    off_diags = (p[:A]-diagm(diag(p[:A]))) * (b.^p[:beta])
    diags = diag(p[:A]) .* (b.^p[:beta])
    bracket = p[:r] - (diags + off_diags)
    db .= b.*bracket
    # so does this
    # thing = p[:A] * (b.^p[:beta])
    # bracket = p[:r] - thing
    # db .= b .* bracket


    ## MING AGAIN ##
    # sum_off_diag = (AA - diagm(diagind(AA))) * (b.^p[:alpha])
    # sum_diag = diagm(diagind(AA)) * (b.^p[:beta])
    # sum_A = sum_off_diag + sum_diag
    # db = b.*(p[:r] - sum_A)
    # return db

    # return b .* bracket

    # bracket = sign(p[:alpha])*p[:z] .- sign(p[:alpha]) .* p[:r] .* (b.^p[:alpha]) - p[:A] * (b.^p[:beta])
                                                                                  # A[diagind] = 0 so mat mul is fine
    # db .= b .* bracket
end

function get_eigvs(sol, p)
    # @bp
    final_state = sol[end]
    final_state_b = final_state.^(p[:beta]-1)
    @ein J[i,j] := p[:A][i,j]*final_state[i] * final_state_b[j]  # Build Jacobian from final solution
    J .*= -1 * p[:beta]
    J[diagind(J)] = -sign(p[:alpha]) * p[:r] *p[:alpha] .* final_state.^p[:alpha]
    eigvs = eigen(J).values
    return eigvs
end


rng = MersenneTwister(42)

plot_font = "Computer Modern"
default(fontfamily=plot_font,
linewidth=2, framestyle=:box, label=nothing, grid=false)


# Debugger.@enter main()
main()