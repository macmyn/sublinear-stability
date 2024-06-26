using DifferentialEquations, Revise, DrWatson, OMEinsum, Random,PyPlot,PyCall, LaTeXStrings
# include(srcdir("NonlinearStability.jl"))

plt.rc("text", usetex=true)
# PyCall.PyDict(plt."rcParams")["text.latex.preamble"] = ["\\usepackage{amsmath}"]

## From pg 43 May Stability/Complexity

function lotka_volterra(y,p,t)
    dy1 = y[1]*(p[:a] - p[:α]*y[2])
    dy2 = y[2]*(-1*p[:b] + p[:β]*y[1])
    dy = [dy1, dy2]
    return dy
end

p = Dict(
    :a => 2,
    :b => 1,
    :α => 1,
    :β => 1,
)
TMAX = 30
u0 = [4,2]
p = NamedTuple([pair for pair in p])  # It doesn't like dicts for some reason
prob = ODEProblem(lotka_volterra, u0,(0,TMAX),p)
sol = solve(prob, Tsit5())

eq1 = p[:b]/p[:β]

plot_ts = 0:0.01:TMAX
sol_t = sol(plot_ts)
sol1 = sol_t[1,:]
sol2 = sol_t[2,:]

# pl = plot(layout=grid(1,2))
fig, axs = plt.subplots(1,2,figsize=(6,3))
axs[1].plot(plot_ts,sol1,label=L"S_1")
axs[1].plot(plot_ts,sol2,label=L"S_2")
axs[2].plot(sol1,sol2)
# display(pl)
xs1 = Array(0:0.25:5)
ys1 = Array(0:0.25:5)
meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))
xs, ys = meshgrid(0:0.5:6,0:0.5:6)

yvec = transpose(hcat(xs,ys))
# transpose!(yvec)

dy = [lotka_volterra(yvec[:,i], p, 0) for i in 1:size(yvec)[2]]
dy = hcat(dy...)
scale = 1
u = dy[1,:]./scale
v = dy[2,:]./scale
axs[2].quiver(xs, ys, u,v,angles="xy", scale_units="xy", scale=40)

axs[1].set_xlabel(L"t")
axs[1].set_ylabel(L"\textrm{Abundance}")
axs[2].set_xlabel(L"S_1")
axs[2].set_ylabel(L"S_2")

PyPlot.savefig(plotsdir("lv.png"),bbox_inches="tight",dpi=500)