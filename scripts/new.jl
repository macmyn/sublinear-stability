using OrdinaryDiffEq, Plots, LinearAlgebra, Random, Distributions, ForwardDiff, OMEinsum, DrWatson, ColorSchemes, FileIO, JLD2
gr()
rng = MersenneTwister(42)

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false)

function new(db, b, p, t)
    bracket = sign(p[:alpha])*p[:z] .- sign(p[:alpha])*p[:r].*(b.^p[:alpha]) - p[:A] * (b.^p[:beta])
    db .= b .* bracket
end


tspan = (0.0, 1.0)
# plot()
plots = plot(layout=(2,1))

d_s = []


alphas = -2:0.1:2
betas = -2:0.1:2

for a in alphas, b in betas
# a = -1
# b = 1
    all_params = Dict{Symbol,Any}(
    :z => 5,
    :r => 1,
    :alpha => a,
    :beta => b,
    :N => [20,50,100],
    :μ => 0.1,
    :σ => 0.05,
    :nothing => ""
    )
    colors = palette(:tab10, length(all_params[:N]))

    dicts = dict_list(all_params::Dict{Symbol,Any})
    maximums = []
    for (i, p) in Iterators.reverse(enumerate(dicts))

        # println(i,p)
        d = Normal(p[:μ], p[:σ])
        A = rand(d, p[:N],p[:N])
        A[diagind(A)] .= 0
        p[:A] = A
        x0 = rand(rng, Uniform(1,5),p[:N])

        p = NamedTuple([pair for pair in p])
        prob = ODEProblem(new, x0, tspan, p,)

        try
        sol = solve(prob, Tsit5())

        label = "\$N = $(p[:N])\$"

        plot!(sol[1:end],subplot=1,label=nothing,color=colors[i],alpha=0.5)
        plot!(sol[1],subplot=1,label=label,color=colors[i],alpha=0.5)

        final_state = sol[end]
        final_state_b = final_state.^(p[:beta]-1)

        @ein J[i,j] := p[:A][i,j]*final_state[i] * final_state_b[j]
        J .*= -1 * p[:beta]
        J[diagind(J)] = -sign(p[:alpha]) * p[:r] *p[:alpha] .* final_state.^p[:alpha]

        eigvs = eigen(J).values
        push!(maximums, maximum(real(eigvs)))
        scatter!(eigvs,subplot=2,label=label,color=colors[i])

        catch
            println("CONFIGURATION: $a, $b failed")
            
            push!(maximums, NaN)
        end
        
    end

    # plots
    # # println(maximums)
    # if NaN in maximums
    #     div_stab = NaN
    # elseif maximums[end] < maximums[1]  # remember we're reversing
    #     div_stab = 'n'  # more species --> less stable
    # else
    #     div_stab = 'y'  # more species --> more stable
    # end

    push!(d_s, (all_params[:alpha], all_params[:beta], maximums))

end

FileIO.save("maximums_2.jld2", d_s)