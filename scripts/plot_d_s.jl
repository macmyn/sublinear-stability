using Plots, CategoricalArrays, FileIO, LaTeXStrings

# maximums = FileIO.load("maximums_2.jld2")["d_s"]
file = datadir("nothing=_r=1_z=1_μ=0.1_σ=0.05.jld2")
maximums = FileIO.load(file)["data"]

# println(maximums)
global d = []
for m in maximums
    println(m)
    maxes = m[3]
    if any(isnan.(maxes))
        div_stab = missing
    elseif maxes[end] < maxes[1]  # remember we're reversing
        div_stab = "n"  # more species --> less stable
        println("less")
    else
        div_stab = "y"  # more species --> more stable
        println("more")

    end
    push!(d, div_stab)
end


coords = [(ds[1], ds[2]) for ds in maximums]
pos_inds = findall(x->x=="y",skipmissing(d))
neg_inds = findall(x->x=="n",skipmissing(d))
nan_inds = findall(ismissing,d)

pos_coords = coords[pos_inds]
neg_coords = coords[neg_inds]
nan_coords = coords[nan_inds]


scatter(xlabel=L"\alpha", ylabel=L"\beta")
scatter!(pos_coords, color="blue", label=L"\mathrm{div} = \mathrm{stab}")
scatter!(neg_coords, color="red", label=L"\mathrm{div} \neq \mathrm{stab}")
scatter!(nan_coords, color="white", label=L"\mathrm{error}")

# marker_map = Dict("y" => "red", "n" => "blue", NaN => "black")
# colours = [marker_map[ds] for ds in d]


