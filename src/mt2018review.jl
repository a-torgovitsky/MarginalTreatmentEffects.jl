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
        opts = opts
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

# Plot MTRs and MTE
function mtr_mte(
    savedir::String,
    filename::String;
    dgp::DGP
)
    # initialize
    m0segments = Vector{Dict}() # keep track of MTR segments for d = 0
    m1segments = Vector{Dict}() # keep track of MTR segments for d = 1
    mtesegments = Vector{Dict}() # keep track of MTE segments

    # get data
    step = 500
    ugrid = (1/step):(1/step):(1 - 1/step)
    ev = DataFrame(z = 1, u = ugrid)
    ev[:, "m0"] = evaluate_mtr(dgp.mtrs[1], ev)
    ev[:, "m1"] = evaluate_mtr(dgp.mtrs[2], ev)
    ev[:, "mte"] = ev[:, "m1"] - ev[:, "m0"]

    m0coordinates = df_to_coordinates(ev, :u, "m0")
    for segment_idx in 1:length(m0coordinates)
        segment = Dict(
            "coordinates" => m0coordinates[segment_idx]
        )
        push!(m0segments, segment)
    end
    m1coordinates = df_to_coordinates(ev, :u, "m1")
    for segment_idx in 1:length(m1coordinates)
        segment = Dict(
            "coordinates" => m1coordinates[segment_idx]
        )
        push!(m1segments, segment)
    end
    mtecoordinates = df_to_coordinates(ev, :u, "mte")
    for segment_idx in 1:length(mtecoordinates)
        segment = Dict(
            "coordinates" => mtecoordinates[segment_idx]
        )
        push!(mtesegments, segment)
    end

    # create tex file
    templatefn = joinpath(savedir, "mt2018review", "tikz-template-mtr.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(
        template;
        m0segments = m0segments,
        m1segments = m1segments,
        mtesegments = mtesegments
    )
    texfn = joinpath(dirname(templatefn), filename * ".tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    return texfn
end

# Plot weights for conventional target parameters
function conventional_weights(
    savedir::String,
    filename::String;
    dgp::DGP
)
    # initialize
    weights = Vector{Dict}()

    # setup
    u₁ = dgp.pscore[findall(dgp.suppZ .== 1)][1]
    u₂ = dgp.pscore[findall(dgp.suppZ .== 2)][1]
    u₃ = dgp.pscore[findall(dgp.suppZ .== 3)][1]
    u₄ = dgp.pscore[findall(dgp.suppZ .== 4)][1]

    # p-score labels on top of plot
    pscore = replace(string(dgp.pscore), "[" => "", "]" => "")
    pscorelabel = join("\$p(" .* string.(collect(1:length(dgp.pscore))) .* ")\$", ",")

    # weights for ATE
    segment_info = df_to_coordinates(
        compute_average_weights(ate(dgp)),
        :u,
        2,
        steps = 1/500
    )
    segments = Vector{Dict}()
    for segment_idx in 1:length(segment_info)
        push!(segments, Dict(
            "coordinates" => segment_info[segment_idx],
            "opts" => ""
        ))
    end
    tp_ate = Dict(
        "legendtitle" => "\$\\text{ATE}\$",
        "segments" => segments
    )
    push!(weights, tp_ate)

    # weights for ATT
    segment_info = df_to_coordinates(
        compute_average_weights(att(dgp), dgp),
        :u,
        2,
        steps = 1/500
    )
    segments = Vector{Dict}()
    for segment_idx in 1:length(segment_info)
        push!(segments, Dict(
            "coordinates" => segment_info[segment_idx],
            "opts" => ifelse(
                segment_idx == length(segment_info),
                "",
                "forget plot"
            )
        ))
    end
    tp_att = Dict(
        "legendtitle" => "\$\\text{ATT}\$",
        "segments" => segments
    )
    push!(weights, tp_att)

    # weights for ATU
    segment_info = df_to_coordinates(
        compute_average_weights(atu(dgp), dgp),
        :u,
        2,
        steps = 1/500
    )
    segments = Vector{Dict}()
    for segment_idx in 1:length(segment_info)
        push!(segments, Dict(
            "coordinates" => segment_info[segment_idx],
            "opts" => ifelse(
                segment_idx == length(segment_info),
                "",
                "forget plot"
            )
        ))
    end
    tp_atu = Dict(
        "legendtitle" => "\$\\text{ATU}\$",
        "segments" => segments
    )
    push!(weights, tp_atu)

    # weights for LATE₁₂
    segment_info = df_to_coordinates(
        compute_average_weights(late(dgp, u₁, u₂)),
        :u,
        2,
        steps = 1/500
    )
    segments = Vector{Dict}()
    for segment_idx in 1:length(segment_info)
        push!(segments, Dict(
            "coordinates" => segment_info[segment_idx],
            "opts" => ifelse(segment_idx == length(segment_info), "", "forget plot")
        ))
    end
    tp_late₁₂ = Dict(
        "legendtitle" => "\$\\text{LATE}_{1 \\rightarrow 2}\$",
        "segments" => segments
    )
    push!(weights, tp_late₁₂)

    # weights for LATE₂₃
    segment_info = df_to_coordinates(
        compute_average_weights(late(dgp, u₂, u₃)),
        :u,
        2,
        steps = 1/500
    )
    segments = Vector{Dict}()
    for segment_idx in 1:length(segment_info)
        push!(segments, Dict(
            "coordinates" => segment_info[segment_idx],
            "opts" => ifelse(segment_idx == length(segment_info), "", "forget plot")
        ))
    end
    tp_late₂₃ = Dict(
        "legendtitle" => "\$\\text{LATE}_{2 \\rightarrow 3}\$",
        "segments" => segments
    )
    push!(weights, tp_late₂₃)

    # weights for LATE₃₄
    segment_info = df_to_coordinates(
        compute_average_weights(late(dgp, u₃, u₄)),
        :u,
        2,
        steps = 1/500
    )
    segments = Vector{Dict}()
    for segment_idx in 1:length(segment_info)
        push!(segments, Dict(
            "coordinates" => segment_info[segment_idx],
            "opts" => ifelse(segment_idx == length(segment_info), "", "forget plot")
        ))
    end
    tp_late₃₄ = Dict(
        "legendtitle" => "\$\\text{LATE}_{3 \\rightarrow 4}\$",
        "segments" => segments
    )
    push!(weights, tp_late₃₄)

    # create tex file
    templatefn = joinpath(savedir, "mt2018review", "tikz-template-weights.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(
        template;
        pscore = pscore,
        pscorelabel = pscorelabel,
        weights = weights
    )
    texfn = joinpath(dirname(templatefn), filename * ".tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    return texfn
end

# Plot extrapolation curves
function late_extrap(
    savedir::String,
    filename::String;
    dgp::DGP
)
    # initialize
    topticks = Vector{Dict}()
    curves = Vector{Dict}()

    # setup
    u₁ = dgp.pscore[findall(dgp.suppZ .== 1)][1]
    u₂ = dgp.pscore[findall(dgp.suppZ .== 2)][1]
    u₃ = dgp.pscore[findall(dgp.suppZ .== 3)][1]
    u₄ = dgp.pscore[findall(dgp.suppZ .== 4)][1]
    
    assumptions = Dict{Symbol, Any}(:lb => 0, :ub => 1, :saturated => true)
    basis = [(bernstein_basis(2), bernstein_basis(2))]

    # get data
    α = collect(0:0.01:0.58)
    push!(α, u₂ - u₁)
    push!(α, u₄ - u₃)
    sort!(unique!(α))

    # Since LATE's are point-identified, we can use the basis from the DGP and
    # saturated IV-like estimands to obtain the LATE's we want.
    results = DataFrame(α = α)
    results[:, "UB: LATE⁻₂₃(α)"] .= NaN
    results[:, "LB: LATE⁻₂₃(α)"] .= NaN
    results[:, "UB: LATE⁺₂₃(α)"] .= NaN
    results[:, "LB: LATE⁺₂₃(α)"] .= NaN
    results[:, "UB: LATEᵖᵐ₂₃(α)"] .= NaN
    results[:, "LB: LATEᵖᵐ₂₃(α)"] .= NaN
    for row in 1:nrow(results)
        lb = u₂ - results[row, :α]
        ub = u₃
        if (0 < lb < 1) & (0 < ub < 1)
            tp = late(dgp, lb, ub)
            r = compute_bounds(tp, basis, assumptions, dgp)
            results[row, 2] = r[:ub]
            results[row, 3] = r[:lb]
        end
        lb = u₂
        ub = u₃ + results[row, :α]
        if (0 < lb < 1) & (0 < ub < 1)
            tp = late(dgp, lb, ub)
            r = compute_bounds(tp, basis, assumptions, dgp)
            results[row, 4] = r[:ub]
            results[row, 5] = r[:lb]
        end
        lb = u₂ - results[row, :α] / 2
        ub = u₃ + results[row, :α] / 2
        # TODO: do I need the conditional if I ensure α is in the right range?
        if (0 < lb < 1) & (0 < ub < 1)
            tp = late(dgp, lb, ub)
            r = compute_bounds(tp, basis, assumptions, dgp)
            results[row, 6] = r[:ub]
            results[row, 7] = r[:lb]
        end
    end

    # late⁻₂₃(α)
    non_nan = findall(.!isnan.(results[:, 2]))
    coordinates = df_to_coordinates(results[non_nan, :], :α, 2; tol = Inf)
    segments = Vector{Dict}()
    for coordinate_idx in 1:length(coordinates)
        segment = Dict(
            "opts" => ifelse(coordinate_idx == length(coordinates), "", "forget plot"),
            "coordinates" => coordinates[coordinate_idx]
        )
        push!(segments, segment)
    end
    push!(curves, Dict(
        "segments" => segments,
        "legendtitle" => "\$\\text{LATE}^{-}_{2 \\rightarrow 3}(\\alpha)\$"
    ))

    # late⁺₂₃(α)
    non_nan = findall(.!isnan.(results[:, 4]))
    coordinates = df_to_coordinates(results[non_nan, :], :α, 4; tol = Inf)
    segments = Vector{Dict}()
    for coordinate_idx in 1:length(coordinates)
        segment = Dict(
            "opts" => ifelse(coordinate_idx == length(coordinates), "", "forget plot"),
            "coordinates" => coordinates[coordinate_idx]
        )
        push!(segments, segment)
    end
    push!(curves, Dict(
        "segments" => segments,
        "legendtitle" => "\$\\text{LATE}^{+}_{2 \\rightarrow 3}(\\alpha)\$"
    ))

    # lateᵖᵐ₂₃(α)
    non_nan = findall(.!isnan.(results[:, 6]))
    coordinates = df_to_coordinates(results[non_nan, :], :α, 6; tol = Inf)
    segments = Vector{Dict}()
    for coordinate_idx in 1:length(coordinates)
        segment = Dict(
            "opts" => ifelse(coordinate_idx == length(coordinates), "", "forget plot"),
            "coordinates" => coordinates[coordinate_idx]
        )
        push!(segments, segment)
    end
    push!(curves, Dict(
        "segments" => segments,
        "legendtitle" => "\$\\text{LATE}^{\\pm}_{2 \\rightarrow 3}(\\alpha)\$"
    ))

    # top ticks
    # BUG: in the original figure, the xlabel is p(1) - p(2), but it should be
    # p(2) - p(1). For now, I am making the mistake in the original figure.
    late₁₃ = Dict(
        "xpos" => round(u₂ - u₁, digits = 2),
        "ypos" => round(results[findall(results[:,:α] .== u₂ - u₁), 2][1],
                        digits = 4),
        "xlabel" => "\$p(1) - p(2)\$",
        "nodelabel" => "late13",
        "xlabelpos" => .09,
        "ylabelpos" => -.37,
        "label" => "\$\\text{LATE}_{1 \\rightarrow 3}\$"
    )
    push!(topticks, late₁₃)
    late₂₄ = Dict(
        "xpos" => round(u₄ - u₃, digits = 2),
        "ypos" => round(results[findall(results[:,:α] .== u₄ - u₃), 4][1],
                        digits = 4),
        "xlabel" => "\$p(4) - p(3)\$",
        "nodelabel" => "late24",
        "xlabelpos" => .20,
        "ylabelpos" => -.20,
        "label" => "\$\\text{LATE}_{2 \\rightarrow 4}\$"
    )
    push!(topticks, late₂₄)

    # xpos and xlabel
    xpos = Vector()
    xlabel = Vector()
    for tick in topticks
        push!(xpos, tick["xpos"])
        push!(xlabel, tick["xlabel"])
    end
    xpos = join(xpos, ",")
    xlabel = join(xlabel, ",")

    # create tex file
    templatefn = joinpath(savedir, "mt2018review",
                          "tikz-template-late-bounds.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(
        template;
        topticks = topticks,
        curves = curves,
        xpos = xpos,
        xlabel = xlabel
    )
    texfn = joinpath(dirname(templatefn), filename * ".tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    return texfn
end

# Plot LATE bounds under different information
function late_information(
    savedir::String,
    filename::String;
    dgp::DGP
)
    # initialize
    topticks = Vector{Dict}()
    curves = Vector{Dict}()

    # setup
    u₁ = dgp.pscore[findall(dgp.suppZ .== 1)][1]
    u₂ = dgp.pscore[findall(dgp.suppZ .== 2)][1]
    u₃ = dgp.pscore[findall(dgp.suppZ .== 3)][1]
    u₄ = dgp.pscore[findall(dgp.suppZ .== 4)][1]
    basis = [(bernstein_basis(9), bernstein_basis(9))]
    first = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslope => true,
        :tslsslopeind => true,
        :decreasing_level => [(1, 0), (1, 1)]
    )
    second = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslope => true,
        :tslsslopeind => true,
        :wald => [(2, 4)],
        :olsslope => true,
        :decreasing_level => [(1, 0), (1, 1)]
    )
    sharp = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => true,
        :decreasing_level => [(1, 0), (1, 1)]
    )

    # get data
    α = collect(0:0.01:0.52)
    push!(α, u₄ - u₃)
    sort!(unique!(α))
    results = DataFrame(α = α)
    results[:, "UB: LATE⁺₂₃(α) First 𝒮"] .= NaN
    results[:, "LB: LATE⁺₂₃(α) First 𝒮"] .= NaN
    results[:, "UB: LATE⁺₂₃(α) Second 𝒮"] .= NaN
    results[:, "LB: LATE⁺₂₃(α) Second 𝒮"] .= NaN
    results[:, "UB: LATE⁺₂₃(α) Sharp 𝒮"] .= NaN
    results[:, "LB: LATE⁺₂₃(α) Sharp 𝒮"] .= NaN

    for row in 1:nrow(results)
        lb = u₂
        ub = u₃ + results[row, :α]
        tp = late(dgp, lb, ub)
        if (0 < lb < 1) & (0 < ub < 1)
            r = compute_bounds(tp, basis, first, dgp)
            results[row, 2] = r[:ub]
            results[row, 3] = r[:lb]
            r = compute_bounds(tp, basis, second, dgp)
            results[row, 4] = r[:ub]
            results[row, 5] = r[:lb]
            r = compute_bounds(tp, basis, sharp, dgp)
            results[row, 6] = r[:ub]
            results[row, 7] = r[:lb]
        end
    end

    # sharp information set
    segments = Vector{Dict}()
    non_nan = findall(.!isnan.(results[:, 6]))
    coordinates = df_to_coordinates(results[non_nan, :], :α, 6; tol = Inf)
    for coordinate_idx in 1:length(coordinates)
        segment = Dict(
            "opts" => "forget plot",
            "coordinates" => coordinates[coordinate_idx]
        )
        push!(segments, segment)
    end
    non_nan = findall(.!isnan.(results[:, 7]))
    coordinates = df_to_coordinates(results[non_nan, :], :α, 7; tol = Inf)
    for coordinate_idx in 1:length(coordinates)
        segment = Dict(
            "opts" => ifelse(coordinate_idx == length(coordinates), "", "forget plot"),
            "coordinates" => coordinates[coordinate_idx]
        )
        push!(segments, segment)
    end
    push!(curves, Dict(
        "segments" => segments,
        "legendtitle" => "Sharp information set"
    ))

    # second information set
    segments = Vector{Dict}()
    non_nan = findall(.!isnan.(results[:, 4]))
    coordinates = df_to_coordinates(results[non_nan, :], :α, 4; tol = Inf)
    for coordinate_idx in 1:length(coordinates)
        segment = Dict(
            "opts" => "forget plot",
            "coordinates" => coordinates[coordinate_idx]
        )
        push!(segments, segment)
    end
    non_nan = findall(.!isnan.(results[:, 5]))
    coordinates = df_to_coordinates(results[non_nan, :], :α, 5; tol = Inf)
    for coordinate_idx in 1:length(coordinates)
        segment = Dict(
            "opts" => ifelse(coordinate_idx == length(coordinates), "", "forget plot"),
            "coordinates" => coordinates[coordinate_idx]
        )
        push!(segments, segment)
    end
    push!(curves, Dict(
        "segments" => segments,
        "legendtitle" => "Second information set"
    ))

    # first information set
    segments = Vector{Dict}()
    non_nan = findall(.!isnan.(results[:, 2]))
    coordinates = df_to_coordinates(results[non_nan, :], :α, 2; tol = Inf)
    for coordinate_idx in 1:length(coordinates)
        segment = Dict(
            "opts" => "forget plot",
            "coordinates" => coordinates[coordinate_idx]
        )
        push!(segments, segment)
    end
    non_nan = findall(.!isnan.(results[:, 3]))
    coordinates = df_to_coordinates(results[non_nan, :], :α, 3; tol = Inf)
    for coordinate_idx in 1:length(coordinates)
        segment = Dict(
            "opts" => ifelse(coordinate_idx == length(coordinates), "", "forget plot"),
            "coordinates" => coordinates[coordinate_idx]
        )
        push!(segments, segment)
    end
    push!(curves, Dict(
        "segments" => segments,
        "legendtitle" => "First information set"
    ))

    # top ticks
    push!(topticks, Dict(
        "xpos" => round(u₄ - u₃, digits = 2),
        "ypos" => round(results[findall(results[:,:α] .== u₄ - u₃), 5][1], digits = 5),
        "xlabel" => "\$p(4) - p(3)\$",
        "nodelabel" => "late24",
        "xlabelpos" => 0.15,
        "ylabelpos" => -0.1,
        "label" => "\$\\text{LATE}_{2 \\rightarrow 4}\$"
    ))

    # xpos and xlabel
    xpos = Vector()
    xlabel = Vector()
    for tick in topticks
        push!(xpos, tick["xpos"])
        push!(xlabel, tick["xlabel"])
    end
    xpos = join(xpos, ",")
    xlabel = join(xlabel, ",")

    # create tex file
    templatefn = joinpath(savedir, "mt2018review",
                          "tikz-template-late-bounds.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(
        template;
        topticks = topticks,
        curves = curves,
        xpos = xpos,
        xlabel = xlabel
    )
    texfn = joinpath(dirname(templatefn), filename * ".tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    return texfn
end

# Plot LATE bounds under different assumptions
function late_assumptions(
    savedir::String,
    filename::String;
    dgp::DGP
)
    # initialize
    topticks = Vector{Dict}()
    curves = Vector{Dict}()

    # setup
    u₁ = dgp.pscore[findall(dgp.suppZ .== 1)][1]
    u₂ = dgp.pscore[findall(dgp.suppZ .== 2)][1]
    u₃ = dgp.pscore[findall(dgp.suppZ .== 3)][1]
    u₄ = dgp.pscore[findall(dgp.suppZ .== 4)][1]

    nondecr = Dict{Symbol, Any}(:lb => 0, :ub => 1, :saturated => true)
    decr = copy(nondecr)
    decr[:decreasing_level] = [(1, 0), (1, 1)]

    k9 = [(bernstein_basis(9), bernstein_basis(9))]

    # get data
    α = collect(0:0.01:0.52)
    push!(α, u₄ - u₃)
    sort!(unique!(α))
    results = DataFrame(α = α)
    results[:, "UB: NP"] .= NaN
    results[:, "LB: NP"] .= NaN
    results[:, "UB: NP, decr"] .= NaN
    results[:, "LB: NP, decr"] .= NaN
    results[:, "UB: K = 9, decr"] .= NaN
    results[:, "LB: K = 9, decr"] .= NaN

    for row in 1:nrow(results)
        lower = u₂
        upper = u₃ + results[row, :α]
        tp = late(dgp, lower, upper)
        if (0 < lower < 1) & (0 < upper < 1)
            knots = vcat(0, 1, dgp.pscore, lower, upper)
            np = [(constantspline_basis(knots), constantspline_basis(knots))]
            r = compute_bounds(tp, np, nondecr, dgp)
            results[row, 2] = r[:ub]
            results[row, 3] = r[:lb]
            r = compute_bounds(tp, np, decr, dgp)
            results[row, 4] = r[:ub]
            results[row, 5] = r[:lb]
            r = compute_bounds(tp, k9, decr, dgp)
            results[row, 6] = r[:ub]
            results[row, 7] = r[:lb]
        end
    end

    # nonparametric
    segments = Vector{Dict}()
    non_nan = findall(.!isnan.(results[:, 2]))
    coordinates = df_to_coordinates(results[non_nan, :], :α, 2; tol = Inf)
    for coordinate_idx in 1:length(coordinates)
        segment = Dict(
            "opts" => "forget plot",
            "coordinates" => coordinates[coordinate_idx]
        )
        push!(segments, segment)
    end
    non_nan = findall(.!isnan.(results[:, 3]))
    coordinates = df_to_coordinates(results[non_nan, :], :α, 3; tol = Inf)
    for coordinate_idx in 1:length(coordinates)
        segment = Dict(
            "opts" => ifelse(coordinate_idx == length(coordinates), "", "forget plot"),
            "coordinates" => coordinates[coordinate_idx]
        )
        push!(segments, segment)
    end
    push!(curves, Dict(
        "segments" => segments,
        "legendtitle" => "Nonparametric"
    ))

    # nonparametric, decreasing
    segments = Vector{Dict}()
    non_nan = findall(.!isnan.(results[:, 4]))
    coordinates = df_to_coordinates(results[non_nan, :], :α, 4; tol = Inf)
    for coordinate_idx in 1:length(coordinates)
        segment = Dict(
            "opts" => "forget plot",
            "coordinates" => coordinates[coordinate_idx]
        )
        push!(segments, segment)
    end
    non_nan = findall(.!isnan.(results[:, 5]))
    coordinates = df_to_coordinates(results[non_nan, :], :α, 5; tol = Inf)
    for coordinate_idx in 1:length(coordinates)
        segment = Dict(
            "opts" => ifelse(coordinate_idx == length(coordinates), "", "forget plot"),
            "coordinates" => coordinates[coordinate_idx]
        )
        push!(segments, segment)
    end
    push!(curves, Dict(
        "segments" => segments,
        "legendtitle" => "Nonparametric, decreasing"
    ))

    # 9th degree, decreasing
    segments = Vector{Dict}()
    non_nan = findall(.!isnan.(results[:, 6]))
    coordinates = df_to_coordinates(results[non_nan, :], :α, 6; tol = Inf)
    for coordinate_idx in 1:length(coordinates)
        segment = Dict(
            "opts" => "forget plot",
            "coordinates" => coordinates[coordinate_idx]
        )
        push!(segments, segment)
    end
    non_nan = findall(.!isnan.(results[:, 7]))
    coordinates = df_to_coordinates(results[non_nan, :], :α, 7; tol = Inf)
    for coordinate_idx in 1:length(coordinates)
        segment = Dict(
            "opts" => ifelse(coordinate_idx == length(coordinates), "", "forget plot"),
            "coordinates" => coordinates[coordinate_idx]
        )
        push!(segments, segment)
    end
    push!(curves, Dict(
        "segments" => segments,
        "legendtitle" => "Ninth-degree polynomial, decreasing"
    ))

    # top ticks
    push!(topticks, Dict(
        "xpos" => round(u₄ - u₃, digits = 2),
        "ypos" => round(results[findall(results[:,:α] .== u₄ - u₃), 5][1], digits = 5),
        "xlabel" => "\$p(4) - p(3)\$",
        "nodelabel" => "late24",
        "xlabelpos" => 0.15,
        "ylabelpos" => 0.05,
        "label" => "\$\\text{LATE}_{2 \\rightarrow 4}\$",
        "bendopts" => "[bend left=10]"
    ))

    # xpos and xlabel
    xpos = Vector()
    xlabel = Vector()
    for tick in topticks
        push!(xpos, tick["xpos"])
        push!(xlabel, tick["xlabel"])
    end
    xpos = join(xpos, ",")
    xlabel = join(xlabel, ",")

    # create tex file
    templatefn = joinpath(savedir, "mt2018review",
                          "tikz-template-late-bounds.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(
        template;
        topticks = topticks,
        curves = curves,
        xpos = xpos,
        xlabel = xlabel
    )
    texfn = joinpath(dirname(templatefn), filename * ".tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    return texfn
end

# Plot bounds for different polynomial degrees
function kbounds(
    savedir::String,
    filename::String;
    dgp::DGP
)
    # initialize
    curves = Vector{Dict}()

    # setup
    nondecr = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslope => true,
        :tslsslopeind => true
    )
    decr = copy(nondecr)
    decr[:decreasing_level] = [(1, 0), (1, 1)]

    # ypos
    tp = att(dgp)
    basis = [(bernstein_basis(2), bernstein_basis(2))]
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => true,
    )
    r = compute_bounds(tp, basis, assumptions, dgp)
    ypos = round(r[:lb], digits = 4)

    # get data
    results = DataFrame(degree = 1:1:19)
    results[:, "LB Polynomial"] .= NaN
    results[:, "UB Polynomial"] .= NaN
    results[:, "LB Nonparametric"] .= NaN
    results[:, "UB Nonparametric"] .= NaN
    results[:, "LB Polynomial, Decreasing"] .= NaN
    results[:, "UB Polynomial, Decreasing"] .= NaN
    results[:, "LB Nonparametric, Decreasing"] .= NaN
    results[:, "UB Nonparametric, Decreasing"] .= NaN
    for row in 1:nrow(results)
        basis = [(bernstein_basis(results[row, :degree]),
                  bernstein_basis(results[row, :degree]))]
        r = compute_bounds(tp, basis, nondecr, dgp)
        results[row, 2:3] = [r[:lb], r[:ub]]
        r = compute_bounds(tp, basis, decr, dgp)
        results[row, 6:7] = [r[:lb], r[:ub]]

        knots = vcat(0, 1, dgp.pscore)
        basis = [(constantspline_basis(knots),
                  constantspline_basis(knots))]
        r = compute_bounds(tp, basis, nondecr, dgp)
        results[row, 4:5] = [r[:lb], r[:ub]]
        r = compute_bounds(tp, basis, decr, dgp)
        results[row, 8:9] = [r[:lb], r[:ub]]
    end

    # Polynomial
    poly_lb = df_to_coordinates(results, :degree, 2; tol = Inf)
    for coordinate_idx in 1:length(poly_lb)
        push!(curves, Dict(
            "opts" => "+[forget plot]",
            "coordinates" => poly_lb[coordinate_idx]
        ))
    end
    poly_ub = df_to_coordinates(results, :degree, 3; tol = Inf)
    for coordinate_idx in 1:length(poly_ub)
        push!(curves, Dict(
            "opts" => "",
            "coordinates" => poly_ub[coordinate_idx],
            "shift" => -1
        ))
    end

    # Nonparametric
    np_lb = df_to_coordinates(results, :degree, 4; tol = Inf)
    for coordinate_idx in 1:length(np_lb)
        push!(curves, Dict(
            "opts" => "+[dashed, forget plot]",
            "coordinates" => np_lb[coordinate_idx]
        ))
    end
    np_ub = df_to_coordinates(results, :degree, 5; tol = Inf)
    for coordinate_idx in 1:length(np_ub)
        push!(curves, Dict(
            "opts" => "+[dashed]",
            "coordinates" => np_ub[coordinate_idx]
        ))
    end

    # Polynomial, Decreasing
    poly_lb_decr = df_to_coordinates(results, :degree, 6; tol = Inf)
    for coordinate_idx in 1:length(poly_lb_decr)
        push!(curves, Dict(
            "opts" => "+[forget plot]",
            "coordinates" => poly_lb_decr[coordinate_idx]
        ))
    end
    poly_ub_decr = df_to_coordinates(results, :degree, 7; tol = Inf)
    for coordinate_idx in 1:length(poly_ub_decr)
        push!(curves, Dict(
            "opts" => "",
            "coordinates" => poly_ub_decr[coordinate_idx],
            "shift" => -2
        ))
    end

    # Nonparametric, Decreasing
    np_lb_decr = df_to_coordinates(results, :degree, 8; tol = Inf)
    for coordinate_idx in 1:length(np_lb_decr)
        push!(curves, Dict(
            "opts" => "+[dashed, forget plot]",
            "coordinates" => np_lb_decr[coordinate_idx]
        ))
    end
    np_ub_decr = df_to_coordinates(results, :degree, 9; tol = Inf)
    for coordinate_idx in 1:length(np_ub_decr)
        push!(curves, Dict(
            "opts" => "+[dashed]",
            "coordinates" => np_ub_decr[coordinate_idx]
        ))
    end

    # Truth
    push!(curves, Dict(
        "opts" => "[dotted, mark=star, forget plot]",
        "coordinates" => join("(" .* string.(collect(1:1:19)) .* ", $ypos)")
    ))

    # create tex file
    templatefn = joinpath(savedir, "mt2018review", "tikz-template-kbounds.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(
        template;
        curves = curves,
        ypos = ypos,
        ylabel = "ATT"
    )
    texfn = joinpath(dirname(templatefn), filename * ".tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    return texfn
end
