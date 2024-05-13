using OrdinaryDiffEq, Plots, LinearAlgebra, Random, Distributions, ForwardDiff, OMEinsum, DrWatson, ColorSchemes, FileIO, JLD2
include(srcdir("NonlinearStability.jl"))
gr()

tspan = (0.0, 1.0)
plots = plot(layout=(2,1))

d_s = []

alphas = 0:0.1:1.4
betas = 1.0:0.1:1.5

all_params = Dict{Symbol,Any}(
    :z => 1,
    :r => 1,
    :N => [20,50,100],
    :μ => 0.1,
    :σ => 0.05,
    :nothing => ""
)

for a in alphas, b in betas

    dicts = dict_list(all_params::Dict{Symbol,Any})
    maximums = []
    for (i, p) in Iterators.reverse(enumerate(dicts))
        p[:alpha] = a
        p[:beta] = b

        # Set interaction matrix (normal with zeros)
        p[:A] = get_interaction_matrix(p)
        
        # Initial conditions
        x0 = rand(rng, Uniform(1,5),p[:N])

        p = NamedTuple([pair for pair in p])  # ODEProblem only takes NamedTuple 🙄
        prob = ODEProblem(general_interactions, x0, tspan, p,)

        try
        sol = solve(prob, AutoTsit5(Rosenbrock23()))

        # Jacobian and eigenvalues
        final_state = sol[end]
        final_state_b = final_state.^(p[:beta]-1)
        @ein J[i,j] := p[:A][i,j]*final_state[i] * final_state_b[j]  # Build Jacobian from final solution
        J .*= -1 * p[:beta]
        J[diagind(J)] = -sign(p[:alpha]) * p[:r] *p[:alpha] .* final_state.^p[:alpha]
        eigvs = eigen(J).values
        
        push!(maximums, maximum(real(eigvs)))
        println(maximum(real(eigvs)))

        catch e
            showerror(stdout, e)
            println("CONFIGURATION: $a, $b failed")
            
            push!(maximums, NaN)
        end
        
    end

    push!(d_s, (a, b, maximums))

end

DrWatson._wsave(s::String, v::Vector) = FileIO.save(s, "data", v)

name = savename(all_params, "jld2")
safesave(datadir(name), d_s)
println("Saved to $(datadir(name))")