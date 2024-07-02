using Plots, Revise, FileIO, DrWatson, ColorSchemes, LaTeXStrings

# plot_font = "Computer Modern"
# default(fontfamily=plot_font,
# linewidth=2, framestyle=:box, label=nothing, grid=false)

fname = "N=[10, 20, 50, 200, 500, 1000, 2000, 4000]_alpha=-1_beta=1_correction=true_init=const_r=1.0_r1=1_r2=0.1_tspan=(0.0, 100.0)_type=apples_pears_μ=0.2_σ=0.01.jld2"
filename = datadir(fname)
params = parse_savename(fname)[2]

eigvs = FileIO.load(filename)["data"]

max_eigvs = maximum.(real.(eigvs))
max_eigvs = reverse(max_eigvs)

Ns = params["N"]
Ns = replace(Ns,"["=>"")
Ns = replace(Ns,"]"=>"")
Ns = [parse(Int,strip(n)) for n in split(Ns,",")]

# function get_ns_as_list(Ns::String)
#     Ns = replace(Ns, "[", "")
#     Ns = replace(Ns, "]", "")


p = plot(layout=grid(1,2,widths=(2/3,1/3)))
colors = palette(:thermal, length(Ns))

for (i, N) in enumerate(Iterators.reverse(Ns))  # calculated in reverse for historical plotting reasons
   scatter!(eigvs[i], label=N,color=colors[i], subplot=1) 
end

plot!(xlim=(-4,0),subplot=1)

scatter!(Ns,max_eigvs,subplot=2)
plot!(xlabel=L"N",ylabel=L"\lambda_\mathrm{max}")
plot!(dpi=500)

DrWatson._wsave(s::String, plot::Plots.Plot) = savefig(plot, s)
# savename = 
safesave(plotsdir(replace(fname,".jld2"=>".png")),p)