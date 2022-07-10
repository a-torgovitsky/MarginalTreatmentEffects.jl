function dgp_econometrica()
    DGP(
        suppZ = reshape([0, 1, 2], :, 1), # reshape returns a 2-dim object,
        densZ = [0.5, 0.4, 0.1],
        pscore = [0.35, 0.6, 0.7],
        mtrs = (
            MTR(bernstein_basis(2), hcat(0.6, 0.4, 0.3)),
            MTR(bernstein_basis(2), hcat(0.75, 0.5, 0.25))
        )
    )
end

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

# Figure 1: DGP MTRs and weights for LATE(0.35, 0.90) and IV Slope
function run_np_ivs_notitle(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    basis = [(constantspline_basis(knots),
              constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslope => true
    )
    opts = defaults_econometrica()
    opts[1][:title] = "~"
    texfn = mtrs_and_weights(savedir, "np-ivs-no-title";
        dgp = dgp,
        tp = late(dgp, 0.35, 0.9),
        basis = basis,
        assumptions = assumptions,
        mtroption = "truth",
        opts = opts
    )
    if compile
        compile_latex(texfn)
    end
end

# Figure 2: maximizing MTRs for LATE(0.35, 0.90) with IV Slope Estimand
function run_np_ivs(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    basis = [(constantspline_basis(knots),
              constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslope => true
    )
    opts = defaults_econometrica()
    opts[1][:title] = "Nonparametric bounds"
    texfn = mtrs_and_weights(savedir, "np-ivs";
        dgp = dgp,
        tp = late(dgp, 0.35, 0.9),
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

# Figure 3: maximizing MTRs for LATE(0.35, 0.90) with IV and OLS Slope Estimand
function run_np_ivs_olss(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    basis = [(constantspline_basis(knots),
              constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslope => true,
        :olsslope => true,
    )
    opts = defaults_econometrica()
    opts[1][:title] = "Nonparametric bounds"
    texfn = mtrs_and_weights(savedir, "np-ivs-olss";
        dgp = dgp,
        tp = late(dgp, 0.35, 0.9),
        basis = basis,
        assumptions = assumptions,
        mtroption = "max",
        opts = opts
    )
    if compile
        compile_latex(texfn)
    end
end

# Figure 4: maximizing MTRs for LATE(0.35, 0.90) with Nonparametric IV Slope
function run_np_ivnps(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    basis = [(constantspline_basis(knots),
              constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslopeind => [2, 3],
    )
    opts = defaults_econometrica()
    opts[1][:title] = "Nonparametric bounds"
    texfn = mtrs_and_weights(savedir, "np-ivnps";
        dgp = dgp,
        tp = late(dgp, 0.35, 0.9),
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

# Figure 5: maximizing MTRs for sharp LATE(0.35, 0.90) bounds
function run_np_sharp(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    basis = [(constantspline_basis(knots),
              constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => true
    )
    opts = defaults_econometrica()
    opts[1][:title] = "Nonparametric bounds"
    texfn = mtrs_and_weights(savedir, "np-sharp";
        dgp = dgp,
        tp = late(dgp, 0.35, 0.9),
        basis = basis,
        assumptions = assumptions,
        mtroption = "max",
        opts = opts
    )
    if compile
        compile_latex(texfn)
    end
end

# Figure 6: maximizing, decr. MTRs for sharp LATE(0.35, 0.90) bounds
function run_np_sharp_decr(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    basis = [(constantspline_basis(knots),
              constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => true,
        :decreasing_level => [(1, 0), (1, 1)]
    )
    opts = defaults_econometrica()
    opts[1][:title] = "Nonparametric bounds, MTRs decreasing"
    texfn = mtrs_and_weights(savedir, "np-sharp-decr";
        dgp = dgp,
        tp = late(dgp, 0.35, 0.9),
        basis = basis,
        assumptions = assumptions,
        mtroption = "max",
        opts = opts
    )
    if compile
        compile_latex(texfn)
    end
end

# Figure 7: maximizing, decr. 9th order MTRs for sharp LATE(0.35, 0.90) bounds
function run_np_sharp_decr_k9(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    basis = [(bernstein_basis(9),
              bernstein_basis(9))];
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => true,
        :decreasing_level => [(1, 0), (1, 1)]
    )
    opts = defaults_econometrica()
    # BUG: Order is 10, but degree is 9!
    # NOTE: in appendix, K refers to degree, not order
    opts[1][:title] = "Order 9 polynomial bounds, MTRs decreasing"
    texfn = mtrs_and_weights(savedir, "np-sharp-decr-k9";
        dgp = dgp,
        tp = late(dgp, 0.35, 0.9),
        basis = basis,
        assumptions = assumptions,
        mtroption = "max",
        opts = opts
    )
    if compile
        compile_latex(texfn)
    end
end

# Figure 8: Bounds on LATE(0.35, ̄u)
function run_tikz_extrapolate(savedir::String, compile::Bool = false)
    texfn = tikz_extrapolate(savedir, "tikz-extrapolate")
    if compile
        compile_latex(texfn)
    end
end

"""
    tikz_extrapolate(savedir::String, filename::String)

This function produces the tex file used to create Figure 8 in MST (2018).
It can't be used to create other/similar figures.
The file path of this tex file will be returned.

# Parameters

- savedir: directory where the tex file will be stored
- filename: name of the tex file
"""
function tikz_extrapolate(savedir::String, filename::String)
    # setup
    dgp = dgp_econometrica()
    basis = [(bernstein_basis(9),
              bernstein_basis(9))];
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => true,
        :decreasing_level => [(1, 0), (1, 1)]
    )

    # initialize
    lbsegments = Vector{Dict}() # keeps track of segments for lower bound
    ubsegments = Vector{Dict}() # keeps track of segments for upper bound
    truesegments = Vector{Dict}() # keeps track of segments for actual value

    α = 0.005
    results = DataFrame(u = (0.35 + α):α:1)
    results[:, "LB"] .= NaN
    results[:, "UB"] .= NaN
    results[:, "Truth"] .= NaN

    # compute results
    for u_idx in 1:nrow(results)
        tp = late(dgp, 0.35, results[u_idx, :u])
        r = compute_bounds(tp, basis, assumptions, dgp);
        results[u_idx, 2:3] = round.([r[:lb], r[:ub]], digits = 4)
        results[u_idx, 4] = round.(eval_tp(tp, [dgp.mtrs], dgp), digits = 4)
    end

    # convert to coordinates
    lbcoordinates = df_to_coordinates(results, :u, "LB"; tol = Inf)
    for coordinate_idx in 1:length(lbcoordinates)
        segment = Dict(
            "coordinates" => lbcoordinates[coordinate_idx]
        )
        push!(lbsegments, segment)
    end
    ubcoordinates = df_to_coordinates(results, :u, "UB"; tol = Inf)
    for coordinate_idx in 1:length(ubcoordinates)
        segment = Dict(
            "coordinates" => ubcoordinates[coordinate_idx]
        )
        push!(ubsegments, segment)
    end
    truecoordinates = df_to_coordinates(results, :u, "Truth"; tol = Inf)
    for coordinate_idx in 1:length(truecoordinates)
        segment = Dict(
            "coordinates" => truecoordinates[coordinate_idx]
        )
        push!(truesegments, segment)
    end

    # create tex file
    templatefn = joinpath(savedir, "mst2018econometrica",
                          "tikz-template-extrapolate.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(
        template;
        lbsegments = lbsegments,
        ubsegments = ubsegments,
        truesegments = truesegments
    )
    texfn = joinpath(dirname(templatefn), filename * ".tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    return texfn
end
