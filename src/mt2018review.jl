function dgp_review()
    DGP(
      suppZ = reshape([1, 2, 3, 4], :, 1), # reshape returns a 2-dim object
      densZ = [0.25, 0.25, 0.25, 0.25],
      pscore = [0.12, 0.29, 0.48, 0.78],
      mtrs = (
        MTR(bernstein_basis(2), hcat(0.9, 0.35, 0.1)),
        MTR(bernstein_basis(2), hcat(0.35, 0.2, 0))
      )
    )
end

function defaults_review()
    settings = Dict(
        :axisheight => "2in",
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
        :ylabelweights =>  "Weights (where \$\\neq 0\$)"
    )

    colors = ["gray", "teal", "orange", "cyan!60!black", "red!70!white",
              "lime!80!black", "red", "yellow!60!black", "pink",
              "magenta", "teal"]
    marks = repeat([""], 7)
    marksize = repeat([""], 7)
    linetype = vcat("dotted", repeat(["solid"], 6))

    return settings, colors, marks, marksize, linetype
end

# Figure 1: MTRs and MTE
function run_tikz_mtr(savedir::String, compile::Bool = false)
    dgp = dgp_review()
    texfn = mtr_mte(savedir, "tikz-mtr"; dgp = dgp)
    if compile
        compile_latex(texfn)
    end
end

# Figure 2: weights for conventional target parameters
function run_tikz_weights(savedir::String, compile::Bool = false)
    dgp = dgp_review()
    texfn = conventional_weights(savedir, "tikz-weights"; dgp = dgp)
    if compile
        compile_latex(texfn)
    end
end

# Figure 3: LATE extrapolation
function run_tikz_late_extrap(savedir::String, compile::Bool = false)
    dgp = dgp_review()
    texfn = late_extrap(savedir, "tikz-late-extrap"; dgp = dgp)
    if compile
        compile_latex(texfn)
    end
end

# Figure 4: MTRs and Weights for ATT with IV Slope and TSLS Slope, 4th degree
function run_k4(savedir::String, compile::Bool = false)
    dgp = dgp_review()
    basis = [(bernstein_basis(4), bernstein_basis(4))]
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslope => true,
        :tslsslopeind => true
    );
    opts = defaults_review()
    opts[1][:title] = "Bounds"
    opts[1][:titlesuffix] = " -- shown at upper bound"
    texfn = mtrs_and_weights(savedir, "k4";
        dgp = dgp,
        tp = att(dgp),
        basis = basis,
        assumptions = assumptions,
        mtroption = "max",
        opts = opts
    )
    if compile
        compile_latex(texfn)
    end
end

# Figure 5: MTRs and Weights for ATT with IV Slope and TSLS Slope, 9th degree
function run_k9(savedir::String, compile::Bool = false)
    dgp = dgp_review()
    basis = [(bernstein_basis(9), bernstein_basis(9))]
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslope => true,
        :tslsslopeind => true
    );
    opts = defaults_review()
    opts[1][:title] = "Bounds"
    opts[1][:titlesuffix] = " -- shown at upper bound"
    texfn = mtrs_and_weights(savedir, "k9";
        dgp = dgp,
        tp = att(dgp),
        basis = basis,
        assumptions = assumptions,
        mtroption = "max",
        opts = opts
    )
    if compile
        compile_latex(texfn)
    end
end

# Figure 6: Bounds on ATT for different polynomial degrees
function run_kbounds(savedir::String, compile::Bool = false)
    dgp = dgp_review()
    texfn = kbounds(savedir, "kbounds"; dgp = dgp)
    if compile
        compile_latex(texfn)
    end
end

# Figure 7: Nonparametric Bounds on ATT
function run_np(savedir::String, compile::Bool = false)
    dgp = dgp_review()
    knots = vcat(0, 1, dgp.pscore)
    basis = [(constantspline_basis(knots), constantspline_basis(knots))];
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslope => true,
        :tslsslopeind => true
    )
    opts = defaults_review()
    opts[1][:title] = "Bounds"
    opts[1][:titlesuffix] = " -- shown at upper bound"
    texfn = mtrs_and_weights(savedir, "np";
        dgp = dgp,
        tp = att(dgp),
        basis = basis,
        assumptions = assumptions,
        mtroption = "max",
        opts = opts,
        attributes = Dict(
            "LogLevel" => 0,
            "SolveType" => 1
        )
    )
    if compile
        compile_latex(texfn)
    end
end

# Figure 8: 9th-degree, decreasing polynomial bounds on ATT
function run_k9_decr(savedir::String, compile::Bool = false)
    dgp = dgp_review()
    basis = [(bernstein_basis(9), bernstein_basis(9))]
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslope => true,
        :tslsslopeind => true,
        :decreasing_level => [(1, 0), (1, 1)]
    );
    opts = defaults_review()
    opts[1][:title] = "Bounds"
    opts[1][:titlesuffix] = " -- shown at upper bound"
    texfn = mtrs_and_weights(savedir, "k9-decr";
        dgp = dgp,
        tp = att(dgp),
        basis = basis,
        assumptions = assumptions,
        mtroption = "max",
        opts = opts
    )
    if compile
        compile_latex(texfn)
    end
end

# Figure 9: 9th-degree, decreasing polynomial bounds on ATT
function run_k9_decr_add_more(savedir::String, compile::Bool = false)
    dgp = dgp_review()
    basis = [(bernstein_basis(9), bernstein_basis(9))]
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslope => true,
        :tslsslopeind => true,
        :wald => [(2, 4)],
        :olsslope => true,
        :decreasing_level => [(1, 0), (1, 1)]
    );
    opts = defaults_review()
    opts[1][:title] = "Bounds"
    opts[1][:titlesuffix] = " -- shown at upper bound"
    texfn = mtrs_and_weights(savedir, "k9-decr-add-more";
        dgp = dgp,
        tp = att(dgp),
        basis = basis,
        assumptions = assumptions,
        mtroption = "max",
        opts = opts
    )
    if compile
        compile_latex(texfn)
    end
end

# Figure 10: Bounds on LATE⁺₂₃(α) under different information sets
function run_late_bounds_information(savedir::String, compile::Bool = false)
    dgp = dgp_review()
    texfn = late_information(savedir, "late-bounds-information"; dgp = dgp)
    if compile
        compile_latex(texfn)
    end
end

# Figure 11: Bounds on LATE⁺₂₃(α) under different assumptions
function run_late_bounds_assumptions(savedir::String, compile::Bool = false)
    dgp = dgp_review()
    texfn = late_assumptions(savedir, "late-bounds-assumptions"; dgp = dgp)
    if compile
        compile_latex(texfn)
    end
end
