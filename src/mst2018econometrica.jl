# Produce DGP for MST (2018).
function dgp_econometrica()
    return DGP(suppZ = reshape([0, 1, 2], :, 1),
               densZ = [0.5, 0.4, 0.1],
               pscore = [0.35, 0.6, 0.7],
               mtrs = (MTR(bernstein_basis(2), hcat(0.6, 0.4, 0.3)),
                       MTR(bernstein_basis(2), hcat(0.75, 0.5, 0.25))))
end

# Produce default aesthetics for MST (2018).
function defaults_econometrica()
    settings = Dict(
        :axisheight => "2in",
        :axiswidth => "3in",
        :fontsize => "\\large",
        :legendcols => 4,
        :linewidth => "1.5pt",
        :linewidthmtr => "2.1pt",
        :markrepeat => "25",
        :markphase => "12",
        :title => nothing,
        :titlesuffix => "",
        :titlevspace => "5pt",
        :xmin => "0",
        :xmax => "1",
        :weightymin => "-8",
        :weightymax => "8",
        :mtrlegendtext => nothing,
        :mtrymin => "-.01",
        :mtrymax => "1.01",
        :mtrylabel => "MTR",
        :xlabel => "\$u\$",
        :ylabelbuffer => "1.20",
        :ylabeltextwidth => "1in",
        :ylabelweights =>  "Weights (where \$\\neq 0\$)"
    )

    colors = ["gray", "darkgray", "lightgray", "darkgray", "lightgray",
              "darkgray", "lightgray"]
    marks = ["*", "x", "diamond*", "square*", "|", "triangle*", "pentagon*"]
    marksize = ["1pt", "2.5pt", "1.5pt", "1.25pt", "3pt", "1.25pt", "1.25pt"]
    linetype = repeat(["solid"], 7)

    return settings, colors, marks, marksize, linetype
end

# DGP MTRs and Weights for LATE and IV Slope (Figure 1)
function run_np_ivs_notitle(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    bases = [(constantspline_basis(knots), constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(:lb => 0,
                                    :ub => 1,
                                    :saturated => false,
                                    :ivslope => true)
    opts = defaults_econometrica()
    opts[1][:title] = "~" # no title
    texfn = mtrs_and_weights(savedir,
                             "np-ivs-no-title";
                              dgp = dgp,
                              tp = late(dgp, 0.35, 0.9),
                              bases = bases,
                              assumptions = assumptions,
                              mtroption = "truth",
                              opts = opts)
    if compile
        compile_latex(texfn)
    end
end

# LATE Bounds w/ IV Slope (Figure 2)
function run_np_ivs(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    bases = [(constantspline_basis(knots), constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(:lb => 0,
                                    :ub => 1,
                                    :saturated => false,
                                    :ivslope => true)
    opts = defaults_econometrica()
    opts[1][:title] = "Nonparametric bounds"
    texfn = mtrs_and_weights(savedir,
                             "np-ivs";
                             dgp = dgp,
                             tp = late(dgp, 0.35, 0.9),
                             bases = bases,
                             assumptions = assumptions,
                             mtroption = "max",
                             opts = opts,
                             attributes = Dict("LogLevel" => 0,
                                               "SolveType" => 1))
    if compile
        compile_latex(texfn)
    end
end

# LATE Bounds w/ IV & OLS Slopes (Figure 3)
function run_np_ivs_olss(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    bases = [(constantspline_basis(knots), constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(:lb => 0,
                                    :ub => 1,
                                    :saturated => false,
                                    :ivslope => true,
                                    :olsslope => true,)
    opts = defaults_econometrica()
    opts[1][:title] = "Nonparametric bounds"
    texfn = mtrs_and_weights(savedir,
                             "np-ivs-olss";
                             dgp = dgp,
                             tp = late(dgp, 0.35, 0.9),
                             bases = bases,
                             assumptions = assumptions,
                             mtroption = "max",
                             opts = opts)
    if compile
        compile_latex(texfn)
    end
end

# LATE Bounds w/ Nonparametric IV Slopes (Figure 4)
function run_np_ivnps(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    bases = [(constantspline_basis(knots), constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(:lb => 0,
                                    :ub => 1,
                                    :saturated => false,
                                    :ivslopeind => [1, 2])

    # The MTRs that achieve the upper and lower bounds on the target parameter
    # need not be unique. As a result, the plot of the maximizing MTRs may be
    # different when running the code with different solvers or on different
    # machines.
    # In this example, Figure 4 in MST (2018) shows a maximizing MTR for d = 1
    # that takes on the value of 0 on the interval [0, 0.35).
    # However, this is not the only MTR that achieves the upper bound of
    # LATE(0.35, 0.90).
    # In particular, there are no nontrivial weights on this interval, meaning
    # the maximizing MTR can take on any value on this interval without
    # affecting the bounds.
    # Running the same code with a different solver or different machine may
    # result in a different maximizing MTR.
    # Therefore, I produce two versions of the plot.
    # The first version solves the programs without any intervention.
    # The second version adds a constraint to the programs to ensure the
    # maximizing MTR for d = 1 on [0, 0.35) takes on the value of 1.
    # This ensures that the second version matches the original Figure 4 on the
    # the first partition [0, 0.35).
    # More constraints can be added to ensure the maximizing MTR is the same as
    # the one in the original Figure 4.

    # version 1: no additional constraints
    opts = defaults_econometrica()
    opts[1][:title] = "Nonparametric bounds"
    texfn, v1lb, v1ub = mtrs_and_weights(savedir,
                                         "np-ivnps-v1";
                                          dgp = dgp,
                                          tp = late(dgp, 0.35, 0.9),
                                          bases = bases,
                                          assumptions = assumptions,
                                          mtroption = "max",
                                          opts = opts,
                                          attributes = Dict("LogLevel" => 0,
                                                            "SolveType" => 1),
                                          return_bounds = true)
    if compile
        compile_latex(texfn)
    end

    # version 2: add constraint to match original Figure 4
    opts = defaults_econometrica()
    opts[1][:title] = "Nonparametric bounds"
    fixdf = DataFrame(ℓ = 1, d = 1, j = 1, k = 1, fix = 1.0)
    texfn, v2lb, v2ub = mtrs_and_weights(savedir,
                                         "np-ivnps-v2";
                                         dgp = dgp,
                                         tp = late(dgp, 0.35, 0.9),
                                         bases = bases,
                                         assumptions = assumptions,
                                         mtroption = "max",
                                         opts = opts,
                                         attributes = Dict("LogLevel" => 0,
                                                           "SolveType" => 1),
                                         fixdf = fixdf,
                                         return_bounds = true)
    if compile
        compile_latex(texfn)
    end

    # The two versions produce the same upper and lower bounds.
    @assert v1lb == v2lb
    @assert v1ub == v2ub
end

# Sharp LATE Bounds (Figure 5)
function run_np_sharp(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    bases = [(constantspline_basis(knots), constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(:lb => 0,
                                    :ub => 1,
                                    :saturated => true)
    opts = defaults_econometrica()
    opts[1][:title] = "Nonparametric bounds"
    texfn = mtrs_and_weights(savedir,
                             "np-sharp";
                             dgp = dgp,
                             tp = late(dgp, 0.35, 0.9),
                             bases = bases,
                             assumptions = assumptions,
                             mtroption = "max",
                             opts = opts)
    if compile
        compile_latex(texfn)
    end
end

# Sharp LATE Bounds w/ Decr. MTRs (Figure 6)
function run_np_sharp_decr(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    bases = [(constantspline_basis(knots), constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(:lb => 0,
                                    :ub => 1,
                                    :saturated => true,
                                    :decreasing_level => [(1, 0), (1, 1)])
    opts = defaults_econometrica()
    opts[1][:title] = "Nonparametric bounds, MTRs decreasing"
    texfn = mtrs_and_weights(savedir,
                             "np-sharp-decr";
                             dgp = dgp,
                             tp = late(dgp, 0.35, 0.9),
                             bases = bases,
                             assumptions = assumptions,
                             mtroption = "max",
                             opts = opts)
    if compile
        compile_latex(texfn)
    end
end

# Sharp LATE Bounds w/ Decr., 9th degree MTRs (Figure 7)
function run_np_sharp_decr_k9(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    bases = [(bernstein_basis(9), bernstein_basis(9))]
    assumptions = Dict{Symbol, Any}(:lb => 0,
                                    :ub => 1,
                                    :saturated => true,
                                    :decreasing_level => [(1, 0), (1, 1)])
    opts = defaults_econometrica()
    opts[1][:title] = "9th degree polynomial bounds, MTRs decreasing"
    texfn = mtrs_and_weights(savedir,
                             "np-sharp-decr-k9";
                             dgp = dgp,
                             tp = late(dgp, 0.35, 0.9),
                             bases = bases,
                             assumptions = assumptions,
                             mtroption = "max",
                             opts = opts)
    if compile
        compile_latex(texfn)
    end
end

# Bounds on Family of PRTEs (Figure 8)
function run_tikz_extrapolate(savedir::String, compile::Bool = false)
    texfn = tikz_extrapolate(savedir, "tikz-extrapolate")
    if compile
        compile_latex(texfn)
    end
end
