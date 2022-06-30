module MarginalTreatmentEffects

    using MarginalTreatmentEffectsWithMultipleInstruments
    using DataFrames
    using LinearAlgebra
    using Suppressor

    include("utils.jl")
    include("toplevel.jl")

end
