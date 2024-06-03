using Test
using DelimitedFiles
using OMEinsum
using Random
using Plots
using LinearAlgebra
include("..//src/NonlinearStability.jl")

N_SPECIES = 100
N_SPECIES = 100
EXP_INTRA = 1
EXP_INTER = 1
K = 1
a = readdlm("test/a.txt")
N_INITIAL = fill(0.1,N_SPECIES)
R_CONST = fill(1,N_SPECIES)

p = Dict(
    :z => R_CONST,
    :r => 1,
    :A => a,
    :alpha => EXP_INTRA,
    :beta => EXP_INTER
)


general_interactions(1,N_INITIAL,p,1)
# println(d)