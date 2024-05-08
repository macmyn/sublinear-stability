function get_interaction_matrix(p)
    d = Normal(p[:μ], p[:σ])
    A = rand(d, p[:N], p[:N])
    A[diagind(A)] .= 0
    return A
end

function general_interactions(db, b, p, t)
    bracket = sign(p[:alpha])*p[:z] .- sign(p[:alpha]) * p[:r] .* (b.^p[:alpha]) - p[:A] * (b.^p[:beta])
                                                                                  # A[diagind] = 0 so mat mul is fine
    db .= b .* bracket
end

rng = MersenneTwister(42)

plot_font = "Computer Modern"
default(fontfamily=plot_font,
linewidth=2, framestyle=:box, label=nothing, grid=false)
colors = palette(:tab10, length(all_params[:N]))

