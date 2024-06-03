function get_interaction_matrix(p)
    d = Normal(p[:μ], p[:σ])
    A = rand(d, p[:N], p[:N])
    A[diagind(A)] .= 1
    return A
end

function calculate_rfix(p)
    r = sum(p[:A], dims=2)
    return r
end

function equilib_condition(N, p)
    return sign(p[:alpha])*p[:r]*N^p[:alpha] + (p[:N]-1)*p[:μ]*N^p[:beta] - sign(p[:alpha])*p[:z]
end

function solve_initial(p)
    fz = find_zero(f, 0.1)
    return fz    
end

function get_initial_condition(p)
    if p[:init] == "uniform"
        x0 = rand(rng, Uniform(1,5),p[:N])
    elseif p[:init] == "const"
        x0 = fill(0.1,p[:N])
    elseif p[:init] == "solve"
        root = solve_initial(p)
        x0 = fill(root, p[:N])
        # println("xo fill is: $x0")
    else 
        throw(ArgumentError("Incorrect input for p[:init]"))
    end
    return x0
end

function get_z(p)
    if p[:z] == "r-fix"
        return calculate_rfix(p)
    elseif typeof(p[:z]) == Float64
        return p[:z]
    else
        throw(ArgumentError("Incorrect input for p[:z]"))
    end
end

function general_interactions(db, b, p, t)
    # @bp
    # bracket = sign(p[:alpha])*p[:z] .- sign(p[:alpha]) .* p[:r] .* (b.^p[:alpha]) - p[:A] * (b.^p[:beta])
    #                                                                               # A[diagind] = 0 so mat mul is fine
    # db .= b .* bracket

    z_term = sign(p[:alpha]) .* p[:z]  # z is the 'r-fix' term if applicable
    diag_term = sign(p[:alpha]) .* p[:r] .* (b.^p[:alpha])  # we're using r as A_ii
    off_diag_mat = p[:A] - diagm(diag(p[:A]))  # set diag(A) to 0
    off_diag_term = off_diag_mat * (b.^p[:beta])

    bracket = (z_term .- diag_term .- off_diag_term)
    # println(bracket)
    # println(b)
    db .= b .* bracket
end

function get_eigvs(sol, p)
    final_state = sol[end]
    final_state_b = final_state.^(p[:beta]-1)
    # println("got here")
    @ein J[i,j] := p[:A][i,j]*final_state[i] * final_state_b[j]  # Build Jacobian from final solution
    J .*= -1 * p[:beta]
    J[diagind(J)] = -sign(p[:alpha]) * p[:r] *p[:alpha] .* final_state.^p[:alpha]
    eigvs = eigen(J).values
    return eigvs
end


rng = MersenneTwister(42)

plot_font = "Computer Modern"
default(fontfamily=plot_font,
linewidth=2, framestyle=:box, label=nothing, grid=false)


