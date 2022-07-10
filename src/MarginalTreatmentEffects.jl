module MarginalTreatmentEffects

    using MarginalTreatmentEffectsWithMultipleInstruments
    using DataFrames
    using LinearAlgebra
    using Suppressor
    using Mustache
    using Printf

    include("weights.jl")
    include("utils.jl")
    include("plotting-aux.jl")
    include("plotting-main.jl")
    include("mst2018econometrica.jl")
    include("mt2018review.jl")
    include("toplevel.jl")

end
