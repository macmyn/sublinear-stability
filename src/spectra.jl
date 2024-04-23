#= Plots of the spectrum of the community matrix and countour lines are defined. =#

using Revise, LaTeXStrings

function spectrum(p)
    @assert p[:converged]
    return eigen(Jbeta(p[:equilibrium], p)).values
end

T(x, y, n, p) = (q = (p[:μ] - 1 / p[:K]^(2 - p[:k]));
sum(
    n .^ 2 ./ ((x .- q * n .- (p[:k] - 1) .* n .^ (p[:k] - 1)) .^ 2 .+ y^2)
)
)

function boundary(p; overprint=false)

    if !haskey(p, :a)
        random_interactions!(p)
    end
    if !haskey(p, :equilibrium)
        evolve!(p)
    end
    if p[:scaled]
        p[:σ] = p[:σ] / sqrt(p[:S])
    end
    if p[:scaled]
        p[:μ] = p[:μ] / p[:S]
    end

    diversity = p[:S]

    s = spectrum(p)
    println(s)

    n = p[:equilibrium]
    if p[:converged]
        s = spectrum(p)
        println("SPECTRUM:\n$s")
        println(minimum(real.(s)))

        s = s[real.(s).>1*minimum(real.(s))]
        println("made it here")
        println("s is now: $s")
        X = range(1.2 * minimum(real.(s)), 0.8 * maximum(real.(s)); length=100)
        println("and here")
        Y = range(1.2 * minimum(imag.(s)), 1.2 * maximum(imag.(s)); length=100)
        if overprint
            pl=scatter!(s, alpha=0.3, legend=false,
                aspect_ratio=1, grid=false, label="1")#, color=COLOR_SUB35)
        else
            scatter(s, alpha=0.3, legend=false,
                aspect_ratio=1, grid=false, label="2")#, color=COLOR_SUB35)
        end

        # if p[:k] == 1.0  # logisic growth
        #     annotate!(minimum(minimum(real.(s))), maximum(maximum(imag.(s)))+0.05, "S=$diversity")
        #     annotate!(-0.65,-0.17, L"\longrightarrow")
        #     annotate!(-0.65,-0.19, "Increasing diversity")
        #     title!("Logistic growth")

        # elseif p[:k] == 0.75  # sublinear growth
        #     annotate!(minimum(minimum(real.(s))), maximum(maximum(imag.(s)))+0.02, "S=$diversity")
        #     annotate!(-0.15,-0.08, L"\longleftarrow")
        #     annotate!(-0.15,-0.09, "Increasing diversity")
        #     title!("Sublinear growth")
        
        if haskey(p,:betaa)
            annotate!(minimum(minimum(real.(s))), maximum(maximum(imag.(s)))+0.02, "S=$diversity")
        
        end
        
        pl= Plots.contour!(X, Y, (x, y) -> T(x, y, n, p), levels=[1 / (p[:σ])^2],
            linewidth=2, color=:auto, colorbar=false)
                    
        else
        print("The system is not feasible")
    end
    return pl
end

function inset(p)
    if haskey(p,:betaa)
        xs = range(0,1,100)
        ys = betapdf.(p[:betaa],p[:betab],xs)
        plot!(xs,ys,inset=(1,bbox(0.65,0.08,0.25,0.25)),subplot=2,xticks=0:1:1)
    end

end
