using DrWatson, Glob
@quickactivate
foreach(include, glob("*.jl", srcdir()))

using DataFrames, StatsPlots, DataFramesMeta, LaTeXStrings


df = collect_results(datadir("fig2"))
X = @subset(df, :σ .< :μ)

@df X scatter(
    :μ,
    :σ,
    marker_z = :richness
)