module MarginalTreatmentEffects

    using MarginalTreatmentEffectsWithMultipleInstruments
    using DataFrames
    using LinearAlgebra
    using Suppressor
    using Mustache
    using Printf

    include("weights.jl")
    include("utils.jl")
    include("plotting.jl")
    include("mst2018econometrica.jl")
    include("toplevel.jl")

end
