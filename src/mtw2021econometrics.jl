function default_mtrs()
    # Note: Coefficients here are in the Bernstein basis
    # In the paper I wrote it using the standard basis
    mtrbasis = bernstein_basis(1)
    mtr0 = MTR(mtrbasis, hcat(.5, .4))
    mtr1 = MTR(mtrbasis, hcat(.8, .4))
    return (mtr0, mtr1)
end
export default_mtrs

function default_suppZ()
    return [0 0 ;
            0 1 ;
            1 0 ;
            1 1  ]
end

function default_pscore()
    return [.3, .5, .6, .7]
end

function illustration_dgp()
    DGP(
        suppZ = default_suppZ(),
        densZ = fill(.25, 4),
        pscore = default_pscore(),
        mtrs = default_mtrs()
    )
end
export illustration_dgp

function simulation_dgp()
    DGP(
        suppZ = default_suppZ(),
        densZ = [.4, .3, .1, .2],
        pscore = default_pscore(),
        mtrs = default_mtrs()
    )
end
export simulation_dgp

function prte_dgp(; ℓ_gen = 1)
                    # E[Y | G = g]    Pr[G = g]
    n = NamedArray([    .6                .1 ;   # at
                        .6                .1 ;   # ec
                        .1                .4 ;   # 1c
                        .6                .2 ;   # 2c
                        .2                .1 ;   # rc
                        .1                .1  ], # nt
                   (OrderedDict("at" => 1, "ec" => 2, "1c" => 3,
                                "2c" => 4, "rc" => 5, "nt" => 6),
                    OrderedDict("ExpY" => 1, "Pr" => 2)),
                   ("Rows", "Cols"))
    @assert sum(n[:, "Pr"]) == 1
    # default, but writing it out here so construction below is clear
    # suppZ = [0 0 ;
    #          0 1 ;
    #          1 0 ;
    #          1 1  ]
    pscore = [n["at", "Pr"],
              n["at", "Pr"] + n["ec", "Pr"] + n["2c", "Pr"],
              n["at", "Pr"] + n["ec", "Pr"] + n["1c", "Pr"],
              1 - n["nt", "Pr"]]
    knots = vcat(0, pscore, 1)
    @assert length(knots) == 6 # otherwise construction below will fail
    @assert pscore == sort(pscore) # also needed for construction below

    θ = [fill(NaN, (2, length(knots) - 1)) for ℓ in 1:2]
    ############################################################################
    # Mapping from choice groups and their average effects to MTEs averaged
    # over one of the five knot regions.
    #
    # e.g. ℓ = 1, given Z2 = 0
    # if you are in (0, p(0,0)]
    # then you must be an always taker
    # if you are in (p(0,0), p(1,0)]
    # then you could be an eager complier or Z1-complier
    # if you are in (p(1,0), p(1,1)]
    # you could be a reluctant complier, Z2-complier, or never-taker
    #
    # e.g. ℓ = 1, given Z2 = 1
    # if you are in (0, p(0,1)]
    # then you could be an always-taker, eager complier, or Z2-complier
    # if you are in (p(0,1), p(1,1)]
    # then you could be a reluctant complier or a Z1-complier
    # if you are in (p(1,1), 1]
    # then you must be a never-taker
    ############################################################################
    # ℓ = 1
    #
    # Z2 = 0
    θ[1][1, 1]    = n["at", "ExpY"]
    θ[1][1, 2:3] .= group_average(n, ["ec", "1c"])
    θ[1][1, 4:5] .= group_average(n, ["rc", "2c", "nt"])
    # Z2 = 1
    θ[1][2, 1:2] .= group_average(n, ["at", "2c", "ec"])
    θ[1][2, 3:4] .= group_average(n, ["1c", "rc"])
    θ[1][2, 5]    = n["nt", "ExpY"]
    # ℓ = 2
    #
    # Z1 = 0
    θ[2][1, 1]    = n["at", "ExpY"]
    θ[2][1, 2]    = group_average(n, ["ec", "2c"])
    θ[2][1, 3:5] .= group_average(n, ["rc", "1c", "nt"])
    # Z1 = 1
    θ[2][2, 1:3] .= group_average(n, ["at", "1c", "ec"])
    θ[2][2, 4]    = group_average(n, ["2c", "rc"])
    θ[2][2, 5]    = n["nt", "ExpY"]

    # Assume Y(0) = 0 to make life easy
    notℓ = 2 - (ℓ_gen - 1) # = 2 if ℓ_gen is 1, and 1 if ℓ_gen is 2
    mtr0 = MTR(switching_constantspline_basis(knots, notℓ = notℓ),
               zeros(size(θ[1])))
    mtr1 = MTR(switching_constantspline_basis(knots, notℓ = notℓ), θ[ℓ_gen])

    dgp = DGP(suppZ = default_suppZ(),
              pscore = pscore,
              densZ = [.25, 0, .5, .25], # note misnomer in "support" -- oh well
              mtrs = (mtr0, mtr1))
    return dgp
end
export prte_dgp

function group_average(n, glist)
    dot(n[glist, "ExpY"], n[glist, "Pr"])/sum(n[glist, "Pr"])
end

# Illustration of mutual consistency (Figure 2)
function run_illustrate_mc(savedir::String)
    results = illustrate_mc()
    savedir = joinpath(savedir, project)
    fn_results = joinpath(savedir, "illustrate-mc.csv")
    CSV.write(fn_results, results)
    texfn1 = joinpath(savedir, "FigureMTRConsistent.tex")
    texfn2 = joinpath(savedir, "FigureMTRInConsistent.tex")
    return results, texfn1, texfn2
end
export run_illustrate_mc

# Numerical illustration for ATT (Figure 4)
function run_simulation_att(savedir::String)
    results = simulation_att()
    savedir = joinpath(savedir, project)
    fn_results = joinpath(savedir, "simulation-att.csv")
    CSV.write(fn_results, results)
    texfn = joinpath(savedir, "FigureATT.tex")
    return results, texfn
end
export run_simulation_att

# Numerical illustration for LATE(+20) (Figure 5)
function run_simulation_prte(savedir::String)
    results = simulation_prte()
    savedir = joinpath(savedir, project)
    fn_results = joinpath(savedir, "simulation-prte.csv")
    CSV.write(fn_results, results)
    texfn = joinpath(savedir, "FigurePRTE.tex")
    return results, texfn
end
export run_simulation_prte

# PRTE misspecification example (Figure 6)
function run_prte_misspecification(savedir::String)
    results = prte_misspecification()
    savedir = joinpath(savedir, project)
    fn_results = joinpath(savedir, "prte-misspecification.csv")
    CSV.write(fn_results, results)
    texfn = joinpath(savedir, "FigurePRTEMisspecification.tex")
    return results, texfn
end
export run_prte_misspecification
