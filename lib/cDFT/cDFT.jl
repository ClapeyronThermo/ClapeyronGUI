module cDFT

using LinearAlgebra
using ForwardDiff, NLSolvers
using Clapeyron
using Clapeyron: d, N_A, k_B, R̄
using Clapeyron: @comps
import Clapeyron: Solvers, log, sqrt, log1p, PackedVofV, sparse, SparseMatrixCSC
using Clapeyron: assoc_similar, assoc_matrix_solve, assoc_options, assoc_pair_length
using Clapeyron: issite, compute_index, complement_index, indices
using Clapeyron: SingleComp
using GCIdentifier, ChemicalIdentifiers
using StaticArrays
#if !isdefined(Clapeyron,Symbol("@sum"))
    include("utils/sum.jl")
#else
    #using Clapeyron: @sum
#end

include("base/base.jl")
include("fields/fields.jl")
include("models/models.jl")
include("structure/structure.jl")

include("utils/base.jl")
include("utils/connectivity.jl")
include("utils/boundary_conditions.jl")
include("utils/profiles.jl")
include("utils/integrals.jl")
include("utils/anderson.jl") #TODO, port this to Clapeyron

include("propagator/propagator.jl")

include("models/BasicIdeal.jl")

include("models/FMT.jl")
include("models/PCSAFT.jl")
include("models/PPCSAFT.jl")
include("models/gcPPCSAFT.jl")
include("models/hetero_gcPPCSAFT.jl")
include("models/QPPCSAFT.jl")
include("models/SAFTVRMie.jl")

include("methods/converge.jl")
include("methods/surface_tension.jl")
include("methods/interfacial_tension.jl")

println("All imported")


end