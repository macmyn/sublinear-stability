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

using DifferentialEquations
function f(u, p, t) 
    sleep(0.001)
    return 1.01 * u
end
u0 = 1 / 2
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)
@time sol1 = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
println("FIRST SOL\n\n\n")
println(sol1)
MaxTime = 100
sol2 = @timeout MaxTime begin
    sol2 = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
end NaN
println("SECOND SOL\n\n\n")
println(sol2)
println("Given 100 seconds, sol converges, so sol(0.1) = ", sol(0.1))

MaxTime = 0.01
sol3 = @timeout MaxTime begin
    sol3 = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
end NaN
println("THIRD SOL\n\n\n")
println(sol3)

println("\n\n\n")
println(sol2==sol1)

println("Given 0.01 seconds, sol doesn't converge, and sol = ", sol3)
