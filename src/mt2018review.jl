# Produce DGP for MT (2018).
function dgp_review()
    return DGP(suppZ = reshape([1, 2, 3, 4], :, 1),
               densZ = [0.25, 0.25, 0.25, 0.25],
               pscore = [0.12, 0.29, 0.48, 0.78],
               mtrs = (MTR(bernstein_basis(2), hcat(0.9, 0.35, 0.1)),
                       MTR(bernstein_basis(2), hcat(0.35, 0.2, 0))))
end

# Produce default aesthetics for MT (2018).
function defaults_review()
    settings = Dict(:axisheight => "2in",
                    :axiswidth => "3in",
                    :fontsize => "\\large",
                    :legendcols => 4,
                    :linewidth => "1.7pt",
                    :linewidthmtr => "2.1pt",
                    # :markrepeat => "25",
                    # :markphase => "12",
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
                    :mteylabel => "MTE",
                    :xlabel => "\$u\$",
                    :ylabelbuffer => "1.20",
                    :ylabeltextwidth => "1in",
                    :ylabelweights =>  "Weights (where \$\\neq 0\$)")

    colors = ["gray", "teal", "orange", "cyan!60!black", "red!70!white",
              "lime!80!black", "red", "yellow!60!black", "pink",
              "magenta", "teal"]
    marks = repeat([""], 7)
    marksize = repeat([""], 7)
    linetype = vcat("dotted", repeat(["solid"], 6))

    return settings, colors, marks, marksize, linetype
end

# DGP MTRs and MTE (Figure 1)
function run_tikz_mtr(savedir::String)
    return mtr_mte(savedir, "tikz-mtr")
end

# Weights for Conventional Target Parameters (Figure 2)
function run_tikz_weights(savedir::String)
    return conventional_weights(savedir, "figure2")
end

# LATE Extrapolation (Figure 3)
function run_tikz_late_extrap(savedir::String)
    return late_extrap(savedir, "figure3")
end

# ATT Bounds w/ 4th degree MTRs (Figure 4)
function run_k4(savedir::String)
    dgp = dgp_review()
    bases = [(bernstein_basis(4), bernstein_basis(4))]
    assumptions = Dict{Symbol, Any}(:lb => 0,
                                    :ub => 1,
                                    :saturated => false,
                                    :ivslope => true,
                                    :tslsslopeind => true)
    opts = defaults_review()
    opts[1][:title] = "Bounds"
    opts[1][:titlesuffix] = " -- shown at upper bound"
    return mtrs_and_weights(savedir,
                            "figure4";
                            dgp = dgp,
                            tp = att(dgp),
                            bases = bases,
                            assumptions = assumptions,
                            mtroption = "max",
                            opts = opts)
end

# ATT Bounds w/ 9th degree MTRs (Figure 5)
function run_k9(savedir::String)
    dgp = dgp_review()
    bases = [(bernstein_basis(9), bernstein_basis(9))]
    assumptions = Dict{Symbol, Any}(:lb => 0,
                                    :ub => 1,
                                    :saturated => false,
                                    :ivslope => true,
                                    :tslsslopeind => true)
    opts = defaults_review()
    opts[1][:title] = "Bounds"
    opts[1][:titlesuffix] = " -- shown at upper bound"
    return mtrs_and_weights(savedir,
                            "figure5";
                            dgp = dgp,
                            tp = att(dgp),
                            bases = bases,
                            assumptions = assumptions,
                            mtroption = "max",
                            opts = opts)
end

# ATT Bounds w/ Different MTR Assumptions (Figure 6)
function run_kbounds(savedir::String)
    dgp = dgp_review()
    return kbounds(savedir, "figure6"; dgp = dgp)
end

# ATT Bounds w/ Nonparametric MTRs (Figure 7)
function run_np(savedir::String)
    dgp = dgp_review()
    knots = vcat(0, 1, dgp.pscore)
    bases = [(constantspline_basis(knots), constantspline_basis(knots))];
    assumptions = Dict{Symbol, Any}(:lb => 0,
                                    :ub => 1,
                                    :saturated => false,
                                    :ivslope => true,
                                    :tslsslopeind => true)
    opts = defaults_review()
    opts[1][:title] = "Bounds"
    opts[1][:titlesuffix] = " -- shown at upper bound"
    return mtrs_and_weights(savedir,
                            "figure7";
                            dgp = dgp,
                            tp = att(dgp),
                            bases = bases,
                            assumptions = assumptions,
                            mtroption = "max",
                            opts = opts,
                            attributes = Dict("LogLevel" => 0,
                                              "SolveType" => 1))
end

# ATT Bounds w/ Decr., 9th degree MTRs (Figure 8)
function run_k9_decr(savedir::String)
    dgp = dgp_review()
    bases = [(bernstein_basis(9), bernstein_basis(9))]
    assumptions = Dict{Symbol, Any}(:lb => 0,
                                    :ub => 1,
                                    :saturated => false,
                                    :ivslope => true,
                                    :tslsslopeind => true,
                                    :decreasing_level => [(1, 0), (1, 1)])
    opts = defaults_review()
    opts[1][:title] = "Bounds"
    opts[1][:titlesuffix] = " -- shown at upper bound"
    return mtrs_and_weights(savedir,
                            "figure8";
                            dgp = dgp,
                            tp = att(dgp),
                            bases = bases,
                            assumptions = assumptions,
                            mtroption = "max",
                            opts = opts)
end

# ATT Bounds w/ more IV-like estimands (Figure 9)
function run_k9_decr_add_more(savedir::String)
    dgp = dgp_review()
    bases = [(bernstein_basis(9), bernstein_basis(9))]
    assumptions = Dict{Symbol, Any}(:lb => 0,
                                    :ub => 1,
                                    :saturated => false,
                                    :ivslope => true,
                                    :tslsslopeind => true,
                                    :wald => [(2, 4)],
                                    :olsslope => true,
                                    :decreasing_level => [(1, 0), (1, 1)])
    opts = defaults_review()
    opts[1][:title] = "Bounds"
    opts[1][:titlesuffix] = " -- shown at upper bound"
    return mtrs_and_weights(savedir,
                            "figure9";
                            dgp = dgp,
                            tp = att(dgp),
                            bases = bases,
                            assumptions = assumptions,
                            mtroption = "max",
                            opts = opts)
end

# LATE Bounds w/ Different Information Sets (Figure 10)
function run_late_bounds_information(savedir::String)
    return late_information(savedir, "figure10")
end

# LATE Bounds w/ Different MTR Assumptions (Figure 11)
function run_late_bounds_assumptions(savedir::String)
    return late_assumptions(savedir, "figure11")
end
