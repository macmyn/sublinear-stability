using DifferentialEquations, Revise, DrWatson, OMEinsum, Random, PyPlot, Distributions
## NOT used in paper
plt.rc("text", usetex=true)

function noisy_competition(y, p, t)
    k1_noise = p[:μ] + noise(p[:σ_squared])
    k2_noise = p[:μ] + noise(p[:σ_squared])
    dy1 = y[1] * (k1_noise - y[1] - p[:α]*y[2])
    dy2 = y[2] * (k2_noise - y[2] - p[:α]*y[1])
    dy = [dy1,dy2]
    return dy
end

function noise(σ_squared)
    dist = Normal(0,σ_squared)
    return rand(dist,1)[1]
end

p = Dict(
    # :k1 => 1,
    # :k2 => 1,
    :α => 0.8,
    :μ => 0.1,
    :σ_squared => 0.01,
)

TMAX = 5000
u0 = [1,1]
p = NamedTuple([pair for pair in p])
prob = ODEProblem(noisy_competition, u0, (0,TMAX), p)
sol = solve(prob, Tsit5())

# From Appendix B of May, Stability
nstar = p[:μ] / (1+p[:α])
mat = -nstar * [1 p[:α]; p[:α] 1]
eigs = eigen(mat).values
println(eigs)

plot_ts = 5:1:TMAX
sol_t = sol(plot_ts)
sol1 = sol_t[1,:]
sol2 = sol_t[2,:]

fig, axs = plt.subplots(1,2)
axs[1].plot(plot_ts,sol1)
axs[1].plot(plot_ts,sol2)
axs[2].scatter(sol1,sol2, marker=".")
axs[2].set_xlim(0,0.2)
axs[2].set_ylim(0,0.2)
display(fig)