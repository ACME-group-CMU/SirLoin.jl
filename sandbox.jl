# based off of https://github.com/rosswhitfield/ase/blob/master/ase/dft/bee.py

using Pkg; Pkg.activate("./")
using Libxc
using LinearAlgebra
using Random, Distributions
using AtomIO
using DFTK
using DFTK: LibxcFunctional, parse_system
using Unitful, UnitfulAtomic

include("params_mbeef.jl")
using .mBEEFParams 

function get_mbeef_coefs(omega, ens_size=2000; seed=123)
    # for the 64 x 64 (biggest) case, this line takes ~0.6s on my desktop...we can decide if it makes more sense to just store and precompute these SVD's
    _, S, V = svd(omega) 

    Random.seed!(seed)
    draw = rand(Normal(), (length(S), ens_size))
    return sqrt(2) .* (V * sqrt.(diagm(S))) * draw
end
