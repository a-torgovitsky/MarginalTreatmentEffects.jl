import MarginalTreatmentEffects as MTE
using Test
using CSV
using DataFrames

@testset "Numerical illustrations in MST (2018, Econometrica)" begin
    dgp = MTE.dgp_econometrica()
    tp = MTE.late(dgp, 0.35, 0.9)
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    np_bases = [(MTE.constantspline_basis(knots),
                 MTE.constantspline_basis(knots))]
    polybases = degree -> [(MTE.bernstein_basis(degree),
                            MTE.bernstein_basis(degree))]
    baseline = Dict{Symbol, Any}(:lb => 0, :ub => 1, :saturated => false)
    @testset "LATE Bounds w/ IV Slope (Figure 2)" begin
        assumptions = copy(baseline)
        assumptions[:ivslope] = true
        result = MTE.compute_bounds(tp, np_bases, assumptions, dgp)
        @test isapprox(round(result[:lb], digits = 3), -0.421)
        @test isapprox(round(result[:ub], digits = 3),  0.500)
    end
    @testset "LATE Bounds w/ IV & OLS Slopes (Figure 3)" begin
        assumptions = copy(baseline)
        assumptions[:ivslope] = true
        assumptions[:olsslope] = true
        result = MTE.compute_bounds(tp, np_bases, assumptions, dgp)
        @test isapprox(round(result[:lb], digits = 3), -0.411)
        @test isapprox(round(result[:ub], digits = 3),  0.500)
    end
    @testset "LATE Bounds w/ Nonparametric IV Slopes (Figure 4)" begin
        assumptions = copy(baseline)
        assumptions[:ivslopeind] = [1, 2]
        result = MTE.compute_bounds(tp, np_bases, assumptions, dgp)
        @test isapprox(round(result[:lb], digits = 3), -0.320)
        @test isapprox(round(result[:ub], digits = 3),  0.407)
    end
    @testset "Sharp LATE Bounds (Figure 5)" begin
        assumptions = copy(baseline)
        assumptions[:saturated] = true
        result = MTE.compute_bounds(tp, np_bases, assumptions, dgp)
        @test isapprox(round(result[:lb], digits = 3), -0.138)
        @test isapprox(round(result[:ub], digits = 3),  0.407)
    end
    @testset "Sharp LATE Bounds w/ Decr. MTRs (Figure 6)" begin
        assumptions = copy(baseline)
        assumptions[:saturated] = true
        assumptions[:decreasing_level] = [(1, 0), (1, 1)]
        result = MTE.compute_bounds(tp, np_bases, assumptions, dgp)
        @test isapprox(round(result[:lb], digits = 3), -0.095)
        @test isapprox(round(result[:ub], digits = 3),  0.077)
    end
    @testset "Sharp LATE Bounds w/ Decr., 9th degree MTRs (Figure 7)" begin
        assumptions = copy(baseline)
        assumptions[:saturated] = true
        assumptions[:decreasing_level] = [(1, 0), (1, 1)]
        result = MTE.compute_bounds(tp, polybases(9), assumptions, dgp)
        @test isapprox(round(result[:lb], digits = 3), -0.000)
        @test isapprox(round(result[:ub], digits = 3),  0.067)
    end
end

@testset "Numerical illustrations in MT (2018, Ann. Review of Economics)" begin
    dgp = MTE.dgp_review()
    tp = MTE.att(dgp)
    knots = vcat(0, 1, dgp.pscore)
    np_bases = [(MTE.constantspline_basis(knots),
                 MTE.constantspline_basis(knots))]
    polybases = degree -> [(MTE.bernstein_basis(degree),
                            MTE.bernstein_basis(degree))]
    baseline = Dict{Symbol, Any}(:lb => 0, :ub => 1, :saturated => false)
    @testset "ATT Bounds w/ 4th degree MTRs (Figure 4)" begin
        assumptions = copy(baseline)
        assumptions[:ivslope] = true
        assumptions[:tslsslopeind] = true
        result = MTE.compute_bounds(tp, polybases(4), assumptions, dgp)
        @test isapprox(round(result[:lb], digits = 3), -0.494)
        @test isapprox(round(result[:ub], digits = 3), -0.073)
    end
    @testset "ATT Bounds w/ 9th degree MTRs (Figure 5)" begin
        assumptions = copy(baseline)
        assumptions[:ivslope] = true
        assumptions[:tslsslopeind] = true
        result = MTE.compute_bounds(tp, polybases(9), assumptions, dgp)
        @test isapprox(round(result[:lb], digits = 3), -0.537)
        @test isapprox(round(result[:ub], digits = 3),  0.049)
    end
    @testset "ATT Bounds w/ Nonparametric MTRs (Figure 7)" begin
        assumptions = copy(baseline)
        assumptions[:ivslope] = true
        assumptions[:tslsslopeind] = true
        result = MTE.compute_bounds(tp, np_bases, assumptions, dgp)
        @test isapprox(round(result[:lb], digits = 3), -0.587)
        @test isapprox(round(result[:ub], digits = 3),  0.154)
    end
    @testset "ATT Bounds w/ Decr., 9th degree MTRs (Figure 8)" begin
        assumptions = copy(baseline)
        assumptions[:ivslope] = true
        assumptions[:tslsslopeind] = true
        assumptions[:decreasing_level] = [(1, 0), (1, 1)]
        result = MTE.compute_bounds(tp, polybases(9), assumptions, dgp)
        @test isapprox(round(result[:lb], digits = 3), -0.467)
        @test isapprox(round(result[:ub], digits = 3), -0.159)
    end
    @testset "ATT Bounds w/ more IV-like estimands (Figure 9)" begin
        assumptions = copy(baseline)
        assumptions[:ivslope] = true
        assumptions[:tslsslopeind] = true
        assumptions[:wald] = [(2, 4)]
        assumptions[:olsslope] = true
        assumptions[:decreasing_level] = [(1, 0), (1, 1)]
        result = MTE.compute_bounds(tp, polybases(9), assumptions, dgp)
        @test isapprox(round(result[:lb], digits = 3), -0.414)
        @test isapprox(round(result[:ub], digits = 3), -0.275)
    end
end

@testset "Numerical illustrations in MTW (Forthcoming, JOE)" begin
    @testset "Mutual consistency (Figure 2)" begin
        fn = joinpath(@__DIR__, "illustrate-mc.csv")
        results_save = CSV.read(fn, DataFrame)
        results_recreate = MTE.illustrate_mc()
        @test isapprox(results_save, results_recreate)
    end

    @testset "ATT (Figure 4)" begin
        fn = joinpath(@__DIR__, "simulation-att.csv")
        results_save = CSV.read(fn, DataFrame)
        results_recreate = MTE.simulation_att()
        @test isapprox(results_save, results_recreate)
    end

    @testset "LATE(+20) (Figure 5)" begin
        fn = joinpath(@__DIR__, "simulation-prte.csv")
        results_save = CSV.read(fn, DataFrame)
        results_recreate = MTE.simulation_prte()
        @test isapprox(results_save, results_recreate)
    end

    @testset "PRTE misspecification (Figure 6)" begin
        fn = joinpath(@__DIR__, "prte-misspecification.csv")
        results_save = CSV.read(fn, DataFrame)
        results_recreate = MTE.prte_misspecification()
        @test isapprox(results_save[:, 2:end], results_recreate[:, 2:end])
    end
end
