module MarginalTreatmentEffects

    using JuMP
    using Clp
    using DataFrames
    using Parameters
    using LinearAlgebra
    using CSV
    using Suppressor
    using NamedArrays
    using DataStructures
    using Mustache
    using Printf

    include("basis.jl")
    include("mtr.jl")
    include("dgp.jl")
    include("target_parameters.jl")
    include("ivlike.jl")
    include("mutual_consistency.jl")
    include("bounds.jl")
    include("weights.jl")
    include("utils.jl")
    include("plotting-aux.jl")
    include("plotting-main.jl")

    include("mst2018econometrica.jl")
    include("mt2018review.jl")
    include("mtw2021econometrics.jl")
    include("toplevel.jl")

end
