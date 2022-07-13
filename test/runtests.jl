using MarginalTreatmentEffects
using Test
using CSV
using DataFrames

@testset "Numerical illustrations in the paper" begin
    @testset "Mutual consistency (Figure 2)" begin
        fn = joinpath(@__DIR__, "illustrate-mc.csv")
        results_save = CSV.read(fn, DataFrame)
        results_recreate = illustrate_mc()
        @test isapprox(results_save, results_recreate)
    end

    @testset "ATT (Figure 4)" begin
        fn = joinpath(@__DIR__, "simulation-att.csv")
        results_save = CSV.read(fn, DataFrame)
        results_recreate = simulation_att()
        @test isapprox(results_save, results_recreate)
    end

    @testset "LATE(+20) (Figure 5)" begin
        fn = joinpath(@__DIR__, "simulation-prte.csv")
        results_save = CSV.read(fn, DataFrame)
        results_recreate = simulation_prte()
        @test isapprox(results_save, results_recreate)
    end

    @testset "PRTE misspecification (Figure 6)" begin
        fn = joinpath(@__DIR__, "prte-misspecification.csv")
        results_save = CSV.read(fn, DataFrame)
        results_recreate = prte_misspecification()
        @test isapprox(results_save[:, 2:end], results_recreate[:, 2:end])
    end
end
