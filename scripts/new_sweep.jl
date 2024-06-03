using OrdinaryDiffEq, Plots, LinearAlgebra, Random, Distributions, ForwardDiff, OMEinsum, DrWatson, ColorSchemes, FileIO, JLD2, ProgressBars
include(srcdir("NonlinearStability.jl"))
gr()

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

d_s = []

all_params = Dict{Symbol,Any}(
<<<<<<< HEAD
    :as => -3:0.1:3,
    :bs => 0:0.1:3,
=======
    :as => -2.9:0.1:3,
    :bs => -2.9:0.1:3,
>>>>>>> copying_ming
    :z => 1.0,
    :r => 1.0,
    :N => [20,50,100],
    :μ => 0.2,
    :σ => 0.02,
    :tspan => (0.0, 100.0),
    :init => "const",  # "uniform"/"const"/"solve"  
)
global MAXTIME=30
<<<<<<< HEAD
for a in tqdm(all_params[:as]), b in all_params[:bs]
    dicts = dict_list(all_params::Dict{Symbol,Any})
    maximums = []
    # for (i, p) in tqdm(Iterators.reverse(enumerate(dicts)))
    for (i,p) in Iterators.reverse(enumerate(dicts))
=======
for a in all_params[:as], b in all_params[:bs]
    dicts = dict_list(all_params::Dict{Symbol,Any})
    maximums = []
    for (i, p) in tqdm(Iterators.reverse(enumerate(dicts)))
>>>>>>> copying_ming
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
        eigvs = get_eigvs(sol,p)
        
        push!(maximums, maximum(real(eigvs)))
        # println(maximum(real(eigvs)))

        catch e
            showerror(stdout, e)
            println("CONFIGURATION: $a, $b failed")
            
            push!(maximums, NaN)
        end
        
    end

    push!(d_s, (a, b, maximums))

end

DrWatson._wsave(s::String, v::Vector) = FileIO.save(s, "data", v)
DrWatson.default_allowed(::Dict) = (Real, String, Vector, Dict)
name = savename(all_params, "jld2")
safesave(datadir(name), d_s)
println("Saved to $(datadir(name))")