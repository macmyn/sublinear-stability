using DifferentialEquations, Revise, DrWatson, OMEinsum, Random,PyPlot,PyCall, LaTeXStrings
# include(srcdir("NonlinearStability.jl"))

plt.rc("text", usetex=true)
# PyCall.PyDict(plt."rcParams")["text.latex.preamble"] = ["\\usepackage{amsmath}"]

## From pg 43 May Stability/Complexity

function lotka_volterra(y,p,t)
    dy1 = y[1]*p[:r1]*(1 - (y[1] + p[:a12]*y[2])/p[:K1])
    dy2 = y[2]*p[:r2]*(1 - (y[2] + p[:a21]*y[1])/p[:K2])
    
    dy = [dy1, dy2]
    return dy
end

p = Dict(
    :a12 => 0.6,
    :a21 => 0.5,
    :r1 => 1,
    :r2 => 1,
    :K1 => 1,
    :K2 => 1
)
TMAX = 30
u0 = [0.1,0.1]
p = NamedTuple([pair for pair in p])  # It doesn't like dicts for some reason
prob = ODEProblem(lotka_volterra, u0,(0,TMAX),p)
sol = solve(prob, Tsit5())

# eq1 = p[:b]/p[:β]

plot_ts = 0:0.01:TMAX
sol_t = sol(plot_ts)
sol1 = sol_t[1,:]
sol2 = sol_t[2,:]

# pl = plot(layout=grid(1,2))
fig, axs = plt.subplots(1,2,figsize=(6,3))
axs[1].plot(plot_ts,sol1,label=L"S_1")
axs[1].plot(plot_ts,sol2,label=L"S_2")
# axs[2].plot(sol1,sol2)
# display(pl)
meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))
OFFSET = 0.1
XMAX = 2.1
YMAX = 2.1
xs, ys = meshgrid(0+OFFSET:0.25:XMAX+OFFSET,0+OFFSET:0.25:YMAX+OFFSET)

yvec = transpose(hcat(xs,ys))
# transpose!(yvec)

dy = [lotka_volterra(yvec[:,i], p, 0) for i in 1:size(yvec)[2]]
dy = hcat(dy...)
scale = 10
norms = sqrt.(dy[1,:].^2 + dy[2,:].^2)
u = dy[1,:]./norms
v = dy[2,:]./norms
axs[2].quiver(xs, ys, u,v,angles="xy", scale_units="xy", scale=8)

y2(x) = (p[:K1] - x)/p[:a12]
y1(x) = p[:K2] - p[:a21]*x
axs[2].plot(xs,y1.(xs))
axs[2].plot(xs,y2.(xs))

axs[1].set_xlabel(L"t")
axs[1].set_ylabel(L"\textrm{Abundance}")
axs[2].set_xlabel(L"S_1")
axs[2].set_ylabel(L"S_2")
axs[2].set_xlim(0,XMAX)
axs[2].set_ylim(0,YMAX)
display(fig)
PyPlot.savefig(plotsdir("lv_competition.png"),bbox_inches="tight",dpi=500)