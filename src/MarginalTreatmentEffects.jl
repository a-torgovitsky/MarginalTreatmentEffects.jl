module MarginalTreatmentEffects

    using MarginalTreatmentEffectsWithMultipleInstruments
    using DataFrames
    using LinearAlgebra
    using Suppressor
    using Mustache
    using Printf

    include("utils.jl")
    include("mst2018econometrica.jl")
    include("toplevel.jl")

end
