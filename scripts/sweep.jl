using OrdinaryDiffEq, Plots, LinearAlgebra, Random, Distributions, ForwardDiff, OMEinsum, DrWatson, ColorSchemes, FileIO, JLD2, ProgressBars, LaTeXStrings
include(srcdir("NonlinearStability.jl"))
gr()

plots = plot(layout=(2,1))
plot_font = "Computer Modern"
plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false)

d_s = []

all_params = Dict{Symbol,Any}(
    :as => -5:0.2:5,
    :bs => 0:0.2:5,
    :z => 1.0,
    :r => 1.0,
    :N => [10,20],
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
        p[:z] = get_z(p)
        
        p = NamedTuple([pair for pair in p])  # ODEProblem only takes NamedTuple 🙄
        
        prob = ODEProblem(general_interactions, x0, p[:tspan], p)

        try
            # sol = @timeout MAXTIME begin
            #     solve(prob, Tsit5())
            # end NaN
        sol = solve(prob, Tsit5())
        # Jacobian and eigenvalues
        eigvs = get_eigvs_sublinear(sol,p)
        
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