using OrdinaryDiffEq, Plots, LinearAlgebra, Random, Distributions, ForwardDiff, OMEinsum, DrWatson, ColorSchemes, FileIO, JLD2, ProgressBars
include(srcdir("NonlinearStability.jl"))
gr()

plots = plot(layout=(2,1))

d_s = []

all_params = Dict{Symbol,Any}(
    :type => "apples_pears",
    :as => -3:0.1:3,
    :bs => 0:0.1:3,
    :z => 1.0,
    :r => 1.0,
    :r1 => 1.0,
    :r2 => 1.0,
    :N => [20,50,100],
    :μ => 0.2,
    :σ => 0.02,
    :tspan => (0.0, 100.0),
    :init => "const",  # "uniform"/"const"/"solve"  
)

for a in tqdm(all_params[:as]), b in all_params[:bs]
    dicts = dict_list(all_params::Dict{Symbol,Any})
    maximums = []
    # for (i, p) in tqdm(Iterators.reverse(enumerate(dicts)))
    for (i,p) in Iterators.reverse(enumerate(dicts))
        p[:alpha] = a
        p[:beta] = b

        p[:A] = get_interaction_matrix(p)
        x0 = get_initial_condition(p)
        # p[:z] = get_z(p)
        
        p = NamedTuple([pair for pair in p])  # ODEProblem only takes NamedTuple 🙄
        
        prob = ODEProblem(apple_pear_interactions, x0, p[:tspan], p)

        try
            # sol = @timeout MAXTIME begin
            #     solve(prob, Tsit5())
            # end NaN
        sol = solve(prob, AutoTsit5(Rosenbrock23()))
        # Jacobian and eigenvalues
        # eigvs = get_eigvs(sol,p)
        eigvs = get_eigvs_apples_pears(sol,p)
        
        push!(maximums, maximum(real(eigvs)))

        catch e
            showerror(stdout, e)
            println("CONFIGURATION: $a, $b failed")
            
            push!(maximums, NaN)
        end
        
    end

    push!(d_s, (a, b, maximums))

end

name = savename(all_params, "jld2")
safesave(datadir(name), d_s)
println("Saved to $(datadir(name))")