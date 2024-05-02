using OrdinaryDiffEq, Plots
gr()

C_1 = 5.730

u_θ  = [2,4]

tspan = (0.0, 30.0)

function fn(du, u,p,t)
    a, b, c, d = p
    du[1] = a*u[1] - b * u[1]*u[2]
    du[2] = d*u[1]*u[2] - c*u[2]
    return du
end

radioactivadecay(u, p, t) = u*p

# Ks = [0.1,1,10]
plts = plot()
# for K in Ks
params = [1,1,1,1]
# println(K)
prob = ODEProblem(fn, u_θ, tspan, params)
sol = solve(prob, Tsit5())
plot!(sol, idxs=[1,2])
# end
