"""
    mtrs_and_weights(savedir::String,
                     filename::String;
                     dgp::DGP,
                     tp::TargetParameter,
                     bases::Vector{Tuple{MTRBasis, MTRBasis}},
                     assumptions::Dict,
                     mtroption::String,
                     opts::Tuple{Dict,
                                 Vector{String},
                                 Vector{String},
                                 Vector{String},
                                 Vector{String}},
                     attributes::Dict = Dict("LogLevel" => 0),
                     startdf::DataFrame= DataFrame(ℓ = Int64[],
                                                   d = Int64[],
                                                   j = Int64[],
                                                   k = Int64[],
                                                   start = Float64[]),
                     fixdf::DataFrame= DataFrame(ℓ = Int64[],
                                                 d = Int64[],
                                                 j = Int64[],
                                                 k = Int64[],
                                                 fix = Float64[]),
                     return_bounds::Bool = false)

This function produces the tex file used to reproduce MST (2018) and MT (2018).
The file path of this tex file will be returned if `return_bounds` is false.
Otherwise, return a tuple with the file path of the tex file, the lower bound,
and upper bound.

# Parameters

- savedir: directory where the tex file will be stored
- filename: name of the tex file
- dgp: data generating process
- tp: target parameter
- bases: vector of tuples, where each tuple is a pair of MTRBasis structs
- assumptions: dictionary of assumptions, including IV-like estimand
- mtroption: either "truth", "max", or "min" to plot true DGP MTRs, maximizing
    MTRs, or minimizing MTRs, resp.
- opts: aesthetic information organized in
    1. settings::Dict
    2. colors::Vector{String}
    3. marks::Vector{String}
    4. marksize::Vector{String}
    5. linetype::Vector{String}
- attributes: set attributes associated with Clp.jl
- startdf: data frame to warm-start decision variables
- fixdf: data frame to fix decision variables to a specified value
- return_bounds: if true, also return lower and upper bounds; default is false
"""
function mtrs_and_weights(savedir::String,
                          filename::String;
                          dgp::DGP,
                          tp::TargetParameter,
                          bases::Vector{Tuple{MTRBasis, MTRBasis}},
                          assumptions::Dict,
                          mtroption::String,
                          opts::Tuple{Dict,
                                      Vector{String},
                                      Vector{String},
                                      Vector{String},
                                      Vector{String}},
                          attributes::Dict = Dict("LogLevel" => 0),
                          startdf::DataFrame= DataFrame(ℓ = Int64[],
                                                        d = Int64[],
                                                        j = Int64[],
                                                        k = Int64[],
                                                        start = Float64[]),
                          fixdf::DataFrame= DataFrame(ℓ = Int64[],
                                                      d = Int64[],
                                                      j = Int64[],
                                                      k = Int64[],
                                                      fix = Float64[]),
                          return_bounds::Bool = false)
    # Initialize Mustache.jl variables.
    m0segments = Vector{Dict}() # keeps track of MTR segments for d = 0
    m1segments = Vector{Dict}() # keeps track of MTR segments for d = 1
    d0weights = Vector{Dict}()  # keeps track of weight segments for d = 0
    d1weights = Vector{Dict}()  # keeps track of weight segments for d = 1
    legend = Vector{Dict}()     # keeps track of legend entries
    ivlike_weights = Vector()   # used to compute y-axis limits

    # Setup
    settings, colors, marks, marksize, linetype = opts # aesthetic information
    aesthetic_counter = 1       # keeps track of colors, marks, etc.

    # Compute bounds
    result = compute_bounds(tp,
                            bases,
                            assumptions,
                            dgp,
                            attributes,
                            startdf,
                            fixdf)

    # Collect data for MTRs
    if mtroption == "truth"
        mtr0 = dgp.mtrs[1]
        mtr1 = dgp.mtrs[2]
        settings[:mtrlegendtext] = "DGP MTRs"
        settings[:title] = settings[:title] * settings[:titlesuffix]
    elseif mtroption == "max"
        mtr0 = result[:mtr_ub][1][1]
        mtr1 = result[:mtr_ub][1][2]
        settings[:title] = settings[:title] * parse_bounds(result) *
            settings[:titlesuffix]
        settings[:mtrlegendtext] = "Maximizing MTRs"
    elseif mtroption == "min"
        mtr0 = result[:mtr_lb][1][1]
        mtr1 = result[:mtr_lb][1][2]
        settings[:title] = settings[:title] * parse_bounds(result) *
            settings[:titlesuffix]
        settings[:mtrlegendtext] = "Minimizing MTRs"
    else
        @error "Unsupported value for mtroption." mtroption
    end
    step = 500
    ugrid = (1/step):(1/step):(1 - 1/step)
    ev = DataFrame(z = 1, u = ugrid)
    mtr_results = DataFrame(u = (1/step):(1/step):(1 - 1/step))
    mtr_results[:, "mtr0"] = evaluate_mtr(mtr0, ev)
    mtr_results[:, "mtr1"] = evaluate_mtr(mtr1, ev)

    # Store MTR data in dictionary for Mustache.jl
    m0coordinates = df_to_coordinates(mtr_results, :u, :mtr0)
    for coordinate_idx in 1:length(m0coordinates)
        segment = Dict("pathname" => "mtr0" * string(coordinate_idx),
                       "coordinates" => m0coordinates[coordinate_idx])
        push!(m0segments, segment)
    end
    m1coordinates = df_to_coordinates(mtr_results, :u, :mtr1)
    for coordinate_idx in 1:length(m1coordinates)
        segment = Dict("pathname" => "mtr1" * string(coordinate_idx),
                       "coordinates" => m1coordinates[coordinate_idx])
        push!(m1segments, segment)
    end

    # Collect data for target parameter
    push!(legend, Dict("color" => colors[aesthetic_counter],
                       "mark" => marks[aesthetic_counter],
                       "marksize" => marksize[aesthetic_counter],
                       "legendtitle" => legendtitle(tp),
                       "linetype" => linetype[aesthetic_counter]))
    tp_weights = compute_average_weights(tp, dgp)
    tp_d0_coord = df_to_coordinates(tp_weights, :u, 3, steps = 1/step)
    for coordinate_idx in 1:length(tp_d0_coord)
        segment = Dict("pathname" => "d0" * pathtitle(tp) *
                           string(coordinate_idx),
                       "color" => colors[aesthetic_counter],
                       "mark" => marks[aesthetic_counter],
                       "marksize" => marksize[aesthetic_counter],
                       "coordinates" => tp_d0_coord[coordinate_idx],
                       "linetype" => linetype[aesthetic_counter])
        push!(d0weights, segment)
    end
    tp_d1_coord = df_to_coordinates(tp_weights, :u, 2, steps = 1/step)
    for coordinate_idx in 1:length(tp_d1_coord)
        segment = Dict("pathname" => "d1" * pathtitle(tp) *
                           string(coordinate_idx),
                       "color" => colors[aesthetic_counter],
                       "mark" => marks[aesthetic_counter],
                       "marksize" => marksize[aesthetic_counter],
                       "coordinates" => tp_d1_coord[coordinate_idx],
                       "linetype" => linetype[aesthetic_counter])
        push!(d1weights, segment)
    end
    aesthetic_counter += 1

    # Local function to update `legend`, `d1weights`, and `d0weights`
    function update_ivlike_weights(ivlike::IVLike,
                                   legendtitle::String = legendtitle(ivlike),
                                   pathtitle::String = pathtitle(ivlike))
        # Set up legend.
        push!(legend, Dict("color" => colors[aesthetic_counter],
                           "mark" => marks[aesthetic_counter],
                           "marksize" => marksize[aesthetic_counter],
                           "legendtitle" => legendtitle,
                           "linetype" => linetype[aesthetic_counter]))

        # Update `ivlike_weights` to find max & min weights for y-axis limits.
        s_weights = compute_average_weights(ivlike, dgp)
        push!(ivlike_weights, s_weights[:, 2]...)
        push!(ivlike_weights, s_weights[:, 3]...)

        # Update coordinates to plot curves.
        ivlike_d0_coord = Vector()
        d0_coordinates = df_to_coordinates(s_weights, :u, 3, steps = 1/step)
        push!(ivlike_d0_coord, d0_coordinates...)
        ivlike_d1_coord = Vector()
        d1_coordinates = df_to_coordinates(s_weights, :u, 2, steps = 1/step)
        push!(ivlike_d1_coord, d1_coordinates...)

        # Store data about d = 0 weights in dictionary for Mustache.jl.
        for coordinate_idx in 1:length(ivlike_d0_coord)
            segment = Dict(
                "pathname" => "d0" * pathtitle * string(coordinate_idx),
                "color" => colors[aesthetic_counter],
                "mark" => marks[aesthetic_counter],
                "marksize" => marksize[aesthetic_counter],
                "coordinates" => ivlike_d0_coord[coordinate_idx],
                "linetype" => linetype[aesthetic_counter]
            )
            push!(d0weights, segment)
        end

        # Store data about d = 1 weights in dictionary for Mustache.jl.
        for coordinate_idx in 1:length(ivlike_d1_coord)
            segment = Dict(
                "pathname" => "d1" * pathtitle * string(coordinate_idx),
                "color" => colors[aesthetic_counter],
                "mark" => marks[aesthetic_counter],
                "marksize" => marksize[aesthetic_counter],
                "coordinates" => ivlike_d1_coord[coordinate_idx],
                "linetype" => linetype[aesthetic_counter]
            )
            push!(d1weights, segment)
        end
        aesthetic_counter += 1 # update aesthetic counter
    end

    # Collect data for IV-like estimands.
    if haskey(assumptions, :ivslope) && assumptions[:ivslope]
        update_ivlike_weights(ivslope(dgp))
    end
    if haskey(assumptions, :tslsslopeind) && assumptions[:tslsslopeind]
        update_ivlike_weights(tslsslope_indicator(dgp))
    end
    if haskey(assumptions, :wald)
        for p in assumptions[:wald]
            update_ivlike_weights(wald(dgp; z₀ = p[1], z₁ = p[2]))
        end
    end
    if haskey(assumptions, :olsslope) && assumptions[:olsslope]
        update_ivlike_weights(olsslope(dgp))
    end
    if haskey(assumptions, :ivslopeind)
        s_ivlike = ivslope_indicator(dgp; support = assumptions[:ivslopeind])
        for i in 1:length(s_ivlike.params[:support])
            zind = s_ivlike.params[:support][i]
            s = IVLike(
                "IV Slope for 𝟙(Z = " * string(dgp.suppZ[zind]) * ")",
                [s_ivlike.s[i]],
                Dict(:support => [dgp.suppZ[zind]])
            )
            update_ivlike_weights(s,
                                  legendtitle(s_ivlike)[i],
                                  pathtitle(s_ivlike)[i])
        end
    end
    if haskey(assumptions, :saturated) && assumptions[:saturated]
        s_ivlike = make_slist(dgp.suppZ)
        # TODO(omkarakatta): is this generalizable to other DGPs?
        legend_order = [1, 3, 5, 2, 4, 6] # ensures correct legend aesthetics
        for saturated_idx in 1:length(s_ivlike.s)
            s = IVLike(
                "Saturated; index $saturated_idx",
                [s_ivlike.s[legend_order[saturated_idx]]],
                nothing
            )
            update_ivlike_weights(s,
                                  legendtitle(s_ivlike, dgp)[saturated_idx],
                                  pathtitle(s_ivlike)[saturated_idx])
        end
    end

    # Find y-axis limits using weights.
    settings[:weightymax] = ceil(max(tp_weights[:, 2]...,
                                     tp_weights[:, 3]...,
                                     ivlike_weights...)) + 1
    settings[:weightymin] = -1 * settings[:weightymax]

    # Create tex file.
    # NOTE: `project` is defined in global scope
    templatefn = joinpath(savedir, project, "tikz-template.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(template;
                 settings...,
                 m0segments = m0segments,
                 m1segments = m1segments,
                 d0weights = d0weights,
                 d1weights = d1weights,
                 legend = legend)
    texfn = joinpath(dirname(templatefn), filename * ".tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    if return_bounds
        return texfn, result[:lb], result[:ub]
    else
        return texfn
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
    # Initialize Mustache.jl variables.
    lbsegments = Vector{Dict}()   # keeps track of segments for lower bound
    ubsegments = Vector{Dict}()   # keeps track of segments for upper bound
    truesegments = Vector{Dict}() # keeps track of segments for actual value

    # Setup
    dgp = dgp_econometrica()
    bases = [(bernstein_basis(9), bernstein_basis(9))];
    assumptions = Dict{Symbol, Any}(:lb => 0,
                                    :ub => 1,
                                    :saturated => true,
                                    :decreasing_level => [(1, 0), (1, 1)])

    # Get data to produce plots.
    α = 0.005
    results = DataFrame(u = (0.35 + α):α:1)
    results[:, "LB"] .= NaN
    results[:, "UB"] .= NaN
    results[:, "Truth"] .= NaN
    for u_idx in 1:nrow(results)
        tp = late(dgp, 0.35, results[u_idx, :u])
        r = compute_bounds(tp, bases, assumptions, dgp);
        results[u_idx, 2:3] = round.([r[:lb], r[:ub]], digits = 4)
        results[u_idx, 4] = round.(eval_tp(tp, [dgp.mtrs], dgp), digits = 4)
    end

    # Convert to coordinates.
    lbcoordinates = df_to_coordinates(results, :u, "LB"; disctol = Inf)
    for coordinate_idx in 1:length(lbcoordinates)
        segment = Dict("coordinates" => lbcoordinates[coordinate_idx])
        push!(lbsegments, segment)
    end
    ubcoordinates = df_to_coordinates(results, :u, "UB"; disctol = Inf)
    for coordinate_idx in 1:length(ubcoordinates)
        segment = Dict("coordinates" => ubcoordinates[coordinate_idx])
        push!(ubsegments, segment)
    end
    truecoordinates = df_to_coordinates(results, :u, "Truth"; disctol = Inf)
    for coordinate_idx in 1:length(truecoordinates)
        segment = Dict("coordinates" => truecoordinates[coordinate_idx])
        push!(truesegments, segment)
    end

    # Create tex file.
    templatefn = joinpath(savedir, project, "tikz-template-extrapolate.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(template;
                 lbsegments = lbsegments,
                 ubsegments = ubsegments,
                 truesegments = truesegments)
    texfn = joinpath(dirname(templatefn), filename * ".tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    return texfn
end

"""
    mtr_and_mte(savedir::String, filename::String)

This function plots the DGP MTRs and the corresponding MTE.
The file path of this tex file will be returned.

# Parameters

- savedir: directory where the tex file will be stored
- filename: name of the tex file
"""
function mtr_mte(savedir::String, filename::String)
    # Initialize Mustache.jl variables.
    m0segments = Vector{Dict}()  # keep track of MTR segments for d = 0
    m1segments = Vector{Dict}()  # keep track of MTR segments for d = 1
    mtesegments = Vector{Dict}() # keep track of MTE segments

    # Setup
    dgp = dgp_review()

    # Get data to produce plots.
    step = 500
    ugrid = (1/step):(1/step):(1 - 1/step)
    ev = DataFrame(z = 1, u = ugrid)
    ev[:, "m0"] = evaluate_mtr(dgp.mtrs[1], ev)
    ev[:, "m1"] = evaluate_mtr(dgp.mtrs[2], ev)
    ev[:, "mte"] = ev[:, "m1"] - ev[:, "m0"]

    m0coordinates = df_to_coordinates(ev, :u, "m0")
    for segment_idx in 1:length(m0coordinates)
        segment = Dict("coordinates" => m0coordinates[segment_idx])
        push!(m0segments, segment)
    end
    m1coordinates = df_to_coordinates(ev, :u, "m1")
    for segment_idx in 1:length(m1coordinates)
        segment = Dict("coordinates" => m1coordinates[segment_idx])
        push!(m1segments, segment)
    end
    mtecoordinates = df_to_coordinates(ev, :u, "mte")
    for segment_idx in 1:length(mtecoordinates)
        segment = Dict("coordinates" => mtecoordinates[segment_idx])
        push!(mtesegments, segment)
    end

    # Create tex file.
    templatefn = joinpath(savedir, project, "tikz-template-mtr.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(template;
                 m0segments = m0segments,
                 m1segments = m1segments,
                 mtesegments = mtesegments)
    texfn = joinpath(dirname(templatefn), filename * ".tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    return texfn
end

"""
    conventional_weights(savedir::String, filename::String)

This function produces the tex file used to create Figure 2 in MT (2018).
It can't be used to create other/similar figures.
The file path of this tex file will be returned.

This function plots the weights of conventional target parameters:
    1. ATE
    2. ATT
    3. ATU
    4. LATE₁₂
    5. LATE₂₃
    6. LATE₃₄

# Parameters

- savedir: directory where the tex file will be stored
- filename: name of the tex file
"""
function conventional_weights(savedir::String, filename::String)
    # Initialize Mustache.jl variables.
    weights = Vector{Dict}() # keep track of weights for target parameters

    # Setup
    dgp = dgp_review()
    u₁ = dgp.pscore[findall(dgp.suppZ .== 1)][1]
    u₂ = dgp.pscore[findall(dgp.suppZ .== 2)][1]
    u₃ = dgp.pscore[findall(dgp.suppZ .== 3)][1]
    u₄ = dgp.pscore[findall(dgp.suppZ .== 4)][1]
    steps = 500

    # p-score labels on top of plot
    pscore = replace(string(dgp.pscore), "[" => "", "]" => "")
    pscorelabel = join("\$p(" .* string.(collect(1:length(dgp.pscore))) .*
                       ")\$", ",")

    # Update `weights` with coordinate information and legend title.
    # To ensure all segments associated with a target parameter is drawn with
    # the same aesthetics, all but the last segment must be drawn with the
    # "forget plot" option.
    function update_tp_weights(tp::TargetParameter,
                               legendtitle::String = legendtitle(tp))
        segments = Vector{Dict}()
        weight_info = compute_average_weights(tp, dgp)
        segment_info = df_to_coordinates(weight_info, :u, 2, steps = 1/steps)
        for i in 1:length(segment_info)
            opts = ifelse(i == length(segment_info), "", "forget plot")
            segment = Dict("coordinates" => segment_info[i], "opts" => opts)
            push!(segments, segment)
        end
        tp_info = Dict("legendtitle" => legendtitle, "segments" => segments)
        push!(weights, tp_info)
    end

    update_tp_weights(ate(dgp))
    update_tp_weights(att(dgp))
    update_tp_weights(atu(dgp))
    update_tp_weights(late(dgp, u₁, u₂),
                      "\$\\text{LATE}_{1 \\rightarrow 2}\$")
    update_tp_weights(late(dgp, u₂, u₃),
                      "\$\\text{LATE}_{2 \\rightarrow 3}\$")
    update_tp_weights(late(dgp, u₃, u₄),
                      "\$\\text{LATE}_{3 \\rightarrow 4}\$")

    # Create tex file.
    templatefn = joinpath(savedir, project, "tikz-template-weights.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(template;
                 pscore = pscore,
                 pscorelabel = pscorelabel,
                 weights = weights)
    texfn = joinpath(dirname(templatefn), filename * ".tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    return texfn
end

"""
    late_extrap(savedir::String, filename::String)

This function produces the tex file used to create Figure 3 in MT (2018).
It can't be used to create other/similar figures.
The file path of this tex file will be returned.

This function plots the following target parameters:
    1. LATE⁻₂₃(α)
    2. LATE⁺₂₃(α)
    3. LATEᵖᵐ₂₃(α)
where α is the degree of extrapolation.

# Parameters

- savedir: directory where the tex file will be stored
- filename: name of the tex file
"""
function late_extrap(savedir::String, filename::String)
    # Initialize Mustache.jl variables.
    topticks = Vector{Dict}() # keep track of information in top x-axis
    curves = Vector{Dict}()   # keep track of different LATE curves

    # setup
    dgp = dgp_review()
    u₁ = dgp.pscore[findall(dgp.suppZ .== 1)][1]
    u₂ = dgp.pscore[findall(dgp.suppZ .== 2)][1]
    u₃ = dgp.pscore[findall(dgp.suppZ .== 3)][1]
    u₄ = dgp.pscore[findall(dgp.suppZ .== 4)][1]

    # Get data to produce plots.
    α = collect(0:0.01:0.58)
    push!(α, u₂ - u₁)
    push!(α, u₄ - u₃)
    sort!(unique!(α))

    # Since LATE's are point-identified, we can use the bases from the DGP and
    # saturated IV-like estimands to obtain the LATE's we want.
    assumptions = Dict{Symbol, Any}(:lb => 0, :ub => 1, :saturated => true)
    bases = [(bernstein_basis(2), bernstein_basis(2))]
    results = DataFrame(α = α)
    results[:, "UB: LATE⁻₂₃(α)"] .= NaN
    results[:, "LB: LATE⁻₂₃(α)"] .= NaN
    results[:, "UB: LATE⁺₂₃(α)"] .= NaN
    results[:, "LB: LATE⁺₂₃(α)"] .= NaN
    results[:, "UB: LATEᵖᵐ₂₃(α)"] .= NaN
    results[:, "LB: LATEᵖᵐ₂₃(α)"] .= NaN
    for row in 1:nrow(results)
        # Update `results` with upper and lower bounds.
        function update_results(row, lbindex, ubindex)
            if (0 < lower < 1) & (0 < upper < 1)
                tp = late(dgp, lower, upper)
                r = compute_bounds(tp, bases, assumptions, dgp)
                results[row, ubindex] = r[:ub]
                results[row, lbindex] = r[:lb]
            end
        end
        lower = u₂ - results[row, :α]
        upper = u₃
        update_results(row, 3, 2)
        lower = u₂
        upper = u₃ + results[row, :α]
        update_results(row, 5, 4)
        lower = u₂ - results[row, :α] / 2
        upper = u₃ + results[row, :α] / 2
        update_results(row, 7, 6)
    end

    # Update `curves` with legend title and segment-specific information.
    # `disctol = Inf` ⟹ all segments will be drawn on one curve.
    function update_curves(column, legendtitle)
        segments = Vector{Dict}()
        non_nan = findall(.!isnan.(results[:, column]))
        coordinates = df_to_coordinates(results[non_nan, :],
                                        :α,
                                        column;
                                        disctol = Inf)
        for i in 1:length(coordinates)
            opts = ifelse(i == length(coordinates), "", "forget plot")
            segment = Dict("opts" => opts,
                           "coordinates" => coordinates[i])
            push!(segments, segment)
        end
        push!(curves, Dict(
            "segments" => segments,
            "legendtitle" => legendtitle
        ))
    end
    update_curves(2, "\$\\text{LATE}^{-}_{2 \\rightarrow 3}(\\alpha)\$")
    update_curves(4, "\$\\text{LATE}^{+}_{2 \\rightarrow 3}(\\alpha)\$")
    update_curves(6, "\$\\text{LATE}^{\\pm}_{2 \\rightarrow 3}(\\alpha)\$")

    # Set up top ticks.
    xpos = round(u₂ - u₁, digits = 2)
    ypos = round(results[findall(results[:,:α] .== u₂ - u₁), 2][1], digits = 4)
    late₁₃ = Dict("xpos" => xpos,
                  "ypos" => ypos,
                  "xlabel" => "\$p(2) - p(1)\$",
                  "nodelabel" => "late13",
                  "xlabelpos" => .09,
                  "ylabelpos" => -.37,
                  "label" => "\$\\text{LATE}_{1 \\rightarrow 3}\$")
    push!(topticks, late₁₃)
    xpos = round(u₄ - u₃, digits = 2)
    ypos = round(results[findall(results[:,:α] .== u₄ - u₃), 4][1], digits = 4)
    late₂₄ = Dict("xpos" => xpos,
                  "ypos" => ypos,
                  "xlabel" => "\$p(4) - p(3)\$",
                  "nodelabel" => "late24",
                  "xlabelpos" => .20,
                  "ylabelpos" => -.20,
                  "label" => "\$\\text{LATE}_{2 \\rightarrow 4}\$")
    push!(topticks, late₂₄)

    # Set up `xpos` and `xlabel`.
    xpos = Vector()
    xlabel = Vector()
    for tick in topticks
        push!(xpos, tick["xpos"])
        push!(xlabel, tick["xlabel"])
    end
    xpos = join(xpos, ",")
    xlabel = join(xlabel, ",")

    # Create tex file.
    templatefn = joinpath(savedir, project, "tikz-template-late-bounds.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(template;
                 topticks = topticks,
                 curves = curves,
                 xpos = xpos,
                 xlabel = xlabel,
                 cycleshift = 0)
    texfn = joinpath(dirname(templatefn), filename * ".tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    return texfn
end

# Plot bounds for different polynomial degrees
"""
    kbounds(savedir::String, filename::String; dgp::DGP)

This function produces the tex file used to plot the bounds on the ATT under
different assumptions on the MTR functions:
    1. Polynomial basis of degree 1, 2, 3, …, 19
    2. Nonparametric basis
    3. Decreasing, Polynomial basis of degree 1, 2, 3, …, 19
    4. Decreasing, Nonparametric basis

# Parameters

- savedir: directory where the tex file will be stored
- filename: name of the tex file
- dgp: data generating process
"""
function kbounds(savedir::String, filename::String; dgp::DGP)
    # Initialize Mustache.jl variables.
    curves = Vector{Dict}() # Keep track of curves

    # Setup
    nondecr = Dict{Symbol, Any}(:lb => 0,
                                :ub => 1,
                                :saturated => false,
                                :ivslope => true,
                                :tslsslopeind => true)
    decr = copy(nondecr)
    decr[:decreasing_level] = [(1, 0), (1, 1)]

    # Solve program to get `ypos`.
    tp = att(dgp)
    bases = [(bernstein_basis(2), bernstein_basis(2))]
    assumptions = Dict{Symbol, Any}(:lb => 0, :ub => 1, :saturated => true)
    r = compute_bounds(tp, bases, assumptions, dgp)
    ypos = round(r[:lb], digits = 4)

    # Get data to produce plots.
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
        # bounds for polynomial basis
        bases = [(bernstein_basis(results[row, :degree]),
                  bernstein_basis(results[row, :degree]))]
        r = compute_bounds(tp, bases, nondecr, dgp)
        results[row, 2:3] = [r[:lb], r[:ub]]
        r = compute_bounds(tp, bases, decr, dgp)
        results[row, 6:7] = [r[:lb], r[:ub]]
        # bounds for nonparametric basis
        knots = vcat(0, 1, dgp.pscore)
        bases = [(constantspline_basis(knots),
                  constantspline_basis(knots))]
        r = compute_bounds(tp, bases, nondecr, dgp)
        results[row, 4:5] = [r[:lb], r[:ub]]
        r = compute_bounds(tp, bases, decr, dgp)
        results[row, 8:9] = [r[:lb], r[:ub]]
    end

    # Polynomial
    poly_lb = df_to_coordinates(results, :degree, 2; disctol = Inf)
    for coordinate_idx in 1:length(poly_lb)
        push!(curves, Dict("opts" => "+[forget plot]",
                           "coordinates" => poly_lb[coordinate_idx]))
    end
    poly_ub = df_to_coordinates(results, :degree, 3; disctol = Inf)
    for coordinate_idx in 1:length(poly_ub)
        push!(curves, Dict("opts" => "",
                           "coordinates" => poly_ub[coordinate_idx],
                           "shift" => -1))
    end

    # Nonparametric
    np_lb = df_to_coordinates(results, :degree, 4; disctol = Inf)
    for coordinate_idx in 1:length(np_lb)
        push!(curves, Dict("opts" => "+[dashed, forget plot]",
                           "coordinates" => np_lb[coordinate_idx]))
    end
    np_ub = df_to_coordinates(results, :degree, 5; disctol = Inf)
    for coordinate_idx in 1:length(np_ub)
        push!(curves, Dict("opts" => "+[dashed]",
                           "coordinates" => np_ub[coordinate_idx]))
    end

    # Polynomial, Decreasing
    poly_lb_decr = df_to_coordinates(results, :degree, 6; disctol = Inf)
    for coordinate_idx in 1:length(poly_lb_decr)
        push!(curves, Dict("opts" => "+[forget plot]",
                           "coordinates" => poly_lb_decr[coordinate_idx]))
    end
    poly_ub_decr = df_to_coordinates(results, :degree, 7; disctol = Inf)
    for coordinate_idx in 1:length(poly_ub_decr)
        push!(curves, Dict("opts" => "",
                           "coordinates" => poly_ub_decr[coordinate_idx],
                           "shift" => -2))
    end

    # Nonparametric, Decreasing
    np_lb_decr = df_to_coordinates(results, :degree, 8; disctol = Inf)
    for coordinate_idx in 1:length(np_lb_decr)
        push!(curves, Dict("opts" => "+[dashed, forget plot]",
                           "coordinates" => np_lb_decr[coordinate_idx]))
    end
    np_ub_decr = df_to_coordinates(results, :degree, 9; disctol = Inf)
    for coordinate_idx in 1:length(np_ub_decr)
        push!(curves, Dict("opts" => "+[dashed]",
                           "coordinates" => np_ub_decr[coordinate_idx]))
    end

    # Truth
    opts = "[dotted, mark=star, forget plot]"
    coordinates = join("(" .* string.(collect(1:1:19)) .* ", $ypos)")
    push!(curves, Dict("opts" => opts, "coordinates" => coordinates))

    # Create tex file.
    templatefn = joinpath(savedir, project, "tikz-template-kbounds.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(template;
                 curves = curves,
                 ypos = ypos,
                 ylabel = "ATT")
    texfn = joinpath(dirname(templatefn), filename * ".tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    return texfn
end

"""
    late_information(savedir::String, filename::String; dgp::DGP)

This function produces the tex file used to create Figure 10 in MT (2018).
It can't be used to create other/similar figures.
The file path of this tex file will be returned.

This function extrapolates and plots the upper and lower bounds of LATE⁺₂₃(α)
under the assumption of decreasing, 9th degree MTRs using three information
sets:
    1. 1st Information Set: IV Slope, Nonparametric TSLS Slope
    2. 2nd Information Set: IV Slope, Nonparametric TSLS Slope, OLS Slope, Wald
    3. Sharp Information Set: see MST (2018)

# Parameters

- savedir: directory where the tex file will be stored
- filename: name of the tex file
"""
function late_information(savedir::String, filename::String)
    # Initialize Mustache.jl variables.
    topticks = Vector{Dict}() # keep track of information in top x-axis
    curves = Vector{Dict}()   # keep track of different LATE curves

    # Setup
    dgp = dgp_review()
    u₁ = dgp.pscore[findall(dgp.suppZ .== 1)][1]
    u₂ = dgp.pscore[findall(dgp.suppZ .== 2)][1]
    u₃ = dgp.pscore[findall(dgp.suppZ .== 3)][1]
    u₄ = dgp.pscore[findall(dgp.suppZ .== 4)][1]
    bases = [(bernstein_basis(9), bernstein_basis(9))]
    first = Dict{Symbol, Any}(:lb => 0,
                              :ub => 1,
                              :saturated => false,
                              :ivslope => true,
                              :tslsslopeind => true,
                              :decreasing_level => [(1, 0), (1, 1)])
    second = Dict{Symbol, Any}(:lb => 0,
                               :ub => 1,
                               :saturated => false,
                               :ivslope => true,
                               :tslsslopeind => true,
                               :wald => [(2, 4)],
                               :olsslope => true,
                               :decreasing_level => [(1, 0), (1, 1)])
    sharp = Dict{Symbol, Any}(:lb => 0,
                              :ub => 1,
                              :saturated => true,
                              :decreasing_level => [(1, 0), (1, 1)])

    # Get data to produce plots.
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
    lower = u₂
    for row in 1:nrow(results)
        upper = u₃ + results[row, :α]
        tp = late(dgp, lower, upper)
        if (0 < lower < 1) & (0 < upper < 1)
            r = compute_bounds(tp, bases, first, dgp)
            results[row, 2] = r[:ub]
            results[row, 3] = r[:lb]
            r = compute_bounds(tp, bases, second, dgp)
            results[row, 4] = r[:ub]
            results[row, 5] = r[:lb]
            r = compute_bounds(tp, bases, sharp, dgp)
            results[row, 6] = r[:ub]
            results[row, 7] = r[:lb]
        end
    end

    # Update `curves` with legend title and segment-specific information.
    function update_curves(ubindex, lbindex, legendtitle)
        segments = Vector{Dict}()
        non_nan = findall(.!isnan.(results[:, ubindex]))
        coordinates = df_to_coordinates(results[non_nan, :],
                                        :α,
                                        ubindex;
                                        disctol = Inf)
        for i in 1:length(coordinates)
            segment = Dict("opts" => "forget plot",
                           "coordinates" => coordinates[i])
            push!(segments, segment)
        end
        non_nan = findall(.!isnan.(results[:, lbindex]))
        coordinates = df_to_coordinates(results[non_nan, :],
                                        :α,
                                        lbindex;
                                        disctol = Inf)
        for i in 1:length(coordinates)
            opts = ifelse(i == length(coordinates), "", "forget plot")
            segment = Dict("opts" => opts,
                           "coordinates" => coordinates[i])
            push!(segments, segment)
        end
        push!(curves, Dict("segments" => segments,
                           "legendtitle" => legendtitle))
    end
    update_curves(6, 7, "Sharp information set")
    update_curves(4, 5, "Second information set")
    update_curves(2, 3, "First information set")

    # Set up top ticks.
    push!(topticks, Dict(
        "xpos" => round(u₄ - u₃, digits = 2),
        "ypos" => round(results[findall(results[:,:α] .== u₄ - u₃), 5][1],
                        digits = 5),
        "xlabel" => "\$p(4) - p(3)\$",
        "nodelabel" => "late24",
        "xlabelpos" => 0.15,
        "ylabelpos" => -0.1,
        "label" => "\$\\text{LATE}_{2 \\rightarrow 4}\$"
    ))

    # Set up `xpos` and `xlabel`.
    xpos = Vector()
    xlabel = Vector()
    for tick in topticks
        push!(xpos, tick["xpos"])
        push!(xlabel, tick["xlabel"])
    end
    xpos = join(xpos, ",")
    xlabel = join(xlabel, ",")

    # create tex file
    templatefn = joinpath(savedir, project, "tikz-template-late-bounds.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(template;
                 topticks = topticks,
                 curves = curves,
                 xpos = xpos,
                 xlabel = xlabel,
                 cycleshift = 2)
    texfn = joinpath(dirname(templatefn), filename * ".tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    return texfn
end

"""
    late_assumptions(savedir::String, filename::String)

This function produces the tex file used to create Figure 11 in MT (2018).
It can't be used to create other/similar figures.
The file path of this tex file will be returned.

This function extrapolates and plots the sharp bounds for LATE⁺₂₃(α)
under increasingly restrictive assumptions on the MTRs:
    1. Nonparametric
    2. Nonparametric, decreasing
    3. 9th degree polynomial, decreasing

# Parameters

- savedir: directory where the tex file will be stored
- filename: name of the tex file
"""
function late_assumptions(savedir::String, filename::String)
    # Initialize Mustache.jl variables.
    topticks = Vector{Dict}() # keep track of information in top x-axis
    curves = Vector{Dict}()   # keep track of different LATE curves

    # Setup
    dgp = dgp_review()
    u₁ = dgp.pscore[findall(dgp.suppZ .== 1)][1]
    u₂ = dgp.pscore[findall(dgp.suppZ .== 2)][1]
    u₃ = dgp.pscore[findall(dgp.suppZ .== 3)][1]
    u₄ = dgp.pscore[findall(dgp.suppZ .== 4)][1]
    nondecr = Dict{Symbol, Any}(:lb => 0, :ub => 1, :saturated => true)
    decr = copy(nondecr)
    decr[:decreasing_level] = [(1, 0), (1, 1)]
    k9 = [(bernstein_basis(9), bernstein_basis(9))]

    # Get data to produce plots.
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
    lower = u₂
    for row in 1:nrow(results)
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

    # Update `curves` with legend title and segment-specific information.
    function update_curves(ubindex, lbindex, legendtitle)
        segments = Vector{Dict}()
        non_nan = findall(.!isnan.(results[:, ubindex]))
        coordinates = df_to_coordinates(results[non_nan, :],
                                        :α,
                                        ubindex;
                                        disctol = Inf)
        for i in 1:length(coordinates)
            segment = Dict("opts" => "forget plot",
                           "coordinates" => coordinates[i])
            push!(segments, segment)
        end
        non_nan = findall(.!isnan.(results[:, lbindex]))
        coordinates = df_to_coordinates(results[non_nan, :],
                                        :α,
                                        lbindex;
                                        disctol = Inf)
        for i in 1:length(coordinates)
            opts = ifelse(i == length(coordinates), "", "forget plot")
            segment = Dict("opts" => opts,
                           "coordinates" => coordinates[i])
            push!(segments, segment)
        end
        push!(curves,
              Dict("segments" => segments, "legendtitle" => legendtitle))
    end
    update_curves(2, 3, "Nonparametric")
    update_curves(4, 5, "Nonparametric, decreasing")
    update_curves(6, 7, "Ninth-degree polynomial, decreasing")

    # Set up top ticks.
    xpos = round(u₄ - u₃, digits = 2)
    ypos = round(results[findall(results[:,:α] .== u₄ - u₃), 5][1], digits = 5)
    push!(topticks, Dict(
        "xpos" => xpos,
        "ypos" => ypos,
        "xlabel" => "\$p(4) - p(3)\$",
        "nodelabel" => "late24",
        "xlabelpos" => 0.15,
        "ylabelpos" => 0.05,
        "label" => "\$\\text{LATE}_{2 \\rightarrow 4}\$",
        "bendopts" => "[bend left=10]"
    ))

    # Set up `xpos` and `xlabel`.
    xpos = Vector()
    xlabel = Vector()
    for tick in topticks
        push!(xpos, tick["xpos"])
        push!(xlabel, tick["xlabel"])
    end
    xpos = join(xpos, ",")
    xlabel = join(xlabel, ",")

    # Create tex file.
    templatefn = joinpath(savedir, project, "tikz-template-late-bounds.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(template;
                 topticks = topticks,
                 curves = curves,
                 xpos = xpos,
                 xlabel = xlabel,
                 cycleshift = 0)
    texfn = joinpath(dirname(templatefn), filename * ".tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    return texfn
end

"""
    illustrate_mc()

This function produces the data to plot Figure 2 of MTW (2021).
It illustrates the concept of mutual consistency.
It can't be used to create other/similar figures.
"""
function illustrate_mc()
    dgp = illustration_dgp()
    assumptions = Dict{Symbol, Any}(:lb => 0, :ub => 1)
    basis1 = (interacted_bernstein_basis(1, notℓ = 2),
              interacted_bernstein_basis(1, notℓ = 2))
    basis2 = (interacted_bernstein_basis(3, notℓ = 1),
              interacted_bernstein_basis(3, notℓ = 1))
    bases = [basis1, basis2]
    tp = ey1(dgp, ℓ = 2)

    # for plotting
    ev = DataFrame(u = 0:.01:1, d = 1)
    results = DataFrame(u = ev.u)

    # Impose mutual consistency
    assumptions[:mutually_consistent] = true
    r = compute_bounds(tp, bases, assumptions, dgp)
    ev.z = fill([0,0], nrow(ev))
    results[:, "MTR (l=1 d=1 z2=0)"] = evaluate_mtr(r[:mtr_ub][1][2], ev)
    ev.z = fill([0,1], nrow(ev))
    results[:, "MTR (l=1 d=1 z2=1)"] = evaluate_mtr(r[:mtr_ub][1][2], ev)
    ev.z = fill([0,0], nrow(ev))
    results[:, "MTR (l=2 d=1 z1=0) con"] = evaluate_mtr(r[:mtr_ub][2][2], ev)
    ev.z = fill([1,0], nrow(ev))
    results[:, "MTR (l=2 d=1 z1=1) con"] = evaluate_mtr(r[:mtr_ub][2][2], ev)

    # Do not impose mutual consistency. Keep target parameter defined in terms
    # of the second model. Then the max MTR will be a second model MTR that is
    # inconsistent with the first model MTR.
    assumptions[:mutually_consistent] = false
    r = compute_bounds(tp, bases, assumptions, dgp)
    ev.z = fill([0,0], nrow(ev))
    results[:, "MTR (l=2 d=1 z1=0) incon"] = evaluate_mtr(r[:mtr_ub][2][2], ev)
    ev.z = fill([1,0], nrow(ev))
    results[:, "MTR (l=2 d=1 z1=1) incon"] = evaluate_mtr(r[:mtr_ub][2][2], ev)
    return results
end
export illustrate_mc

"""
    simulation_att()

This function produces the data to plot Figure 4 of MTW (2021).
It shows the impact of mutual consistency for the ATT, which is an
instrument-invariant parameter.
It can't be used to create other/similar figures.
"""
function simulation_att()
    dgp = simulation_dgp()

    results = DataFrame(degree = 1:1:9)
    results[:, "LB (l = 1; MC off)"] .= NaN
    results[:, "UB (l = 1; MC off)"] .= NaN
    results[:, "LB (l = 2; MC off)"] .= NaN
    results[:, "UB (l = 2; MC off)"] .= NaN
    results[:, "LB (MC on)"] .= NaN
    results[:, "UB (MC on)"] .= NaN
    assumptions = Dict{Symbol, Any}(:lb => 0, :ub => 1)

    for k in 1:nrow(results)
        bases = [(interacted_bernstein_basis(results[k, :degree], notℓ = 2),
                  interacted_bernstein_basis(results[k, :degree], notℓ = 2)),
                 (interacted_bernstein_basis(results[k, :degree], notℓ = 1),
                  interacted_bernstein_basis(results[k, :degree], notℓ = 1))]

        assumptions[:mutually_consistent] = false
        tp = att(dgp, ℓ = 1)
        r = compute_bounds(tp, bases, assumptions, dgp)
        results[k, 2:3] = [r[:lb], r[:ub]]

        tp = att(dgp, ℓ = 2)
        r = compute_bounds(tp, bases, assumptions, dgp)
        results[k, 4:5] = [r[:lb], r[:ub]]

        assumptions[:mutually_consistent] = true
        r = compute_bounds(tp, bases, assumptions, dgp)
        results[k, 6:7] = [r[:lb], r[:ub]]
    end

    knots = vcat(0, dgp.pscore, 1)
    bases = [(interacted_constantspline_basis(knots, notℓ = 2),
              interacted_constantspline_basis(knots, notℓ = 2)),
             (interacted_constantspline_basis(knots, notℓ = 1),
              interacted_constantspline_basis(knots, notℓ = 1))]

    assumptions[:mutually_consistent] = false
    tp = att(dgp, ℓ = 1)
    r = compute_bounds(tp, bases, assumptions, dgp)
    results[:, "NP LB (l = 1; MC off)"] .= r[:lb]
    results[:, "NP UB (l = 1; MC off)"] .= r[:ub]

    tp = att(dgp, ℓ = 2)
    r = compute_bounds(tp, bases, assumptions, dgp)
    results[:, "NP LB (l = 2; MC off)"] .= r[:lb]
    results[:, "NP UB (l = 2; MC off)"] .= r[:ub]

    assumptions[:mutually_consistent] = true
    r = compute_bounds(tp, bases, assumptions, dgp)
    results[:, "NP LB (MC on)"] .= r[:lb]
    results[:, "NP UB (MC on)"] .= r[:ub]

    return results
end
export simulation_att

"""
    simulation_prte()

This function produces the data to plot Figure 5 of MTW (2021).
It shows that mutual consistency can also help tighten inference on
instrument-dependent parameters.
It can't be used to create other/similar figures.
"""
function simulation_prte()
    dgp = simulation_dgp()
    δ = .2
    tp = prte_plusδpercent(dgp, δ, ℓ = 1)
    knots = vcat(0, dgp.pscore, (1 + δ) .* dgp.pscore, 1)

    results = DataFrame(degree = 1:1:9)
    results[:, "LB (MC off)"] .= NaN
    results[:, "UB (MC off)"] .= NaN
    results[:, "LB (l = 2; NP)"] .= NaN
    results[:, "UB (l = 2; NP)"] .= NaN
    results[:, "LB (l = 2; NP; Decr)"] .= NaN
    results[:, "UB (l = 2; NP; Decr)"] .= NaN
    results[:, "LB (l = 2; NP; Linear)"] .= NaN
    results[:, "UB (l = 2; NP; Linear)"] .= NaN

    for k in 1:nrow(results)
        bases = [(interacted_bernstein_basis(results[k, :degree], notℓ = 2),
                  interacted_bernstein_basis(results[k, :degree], notℓ = 2)),
                 (interacted_constantspline_basis(knots, notℓ = 1),
                  interacted_constantspline_basis(knots, notℓ = 1))]

        assumptions = Dict{Symbol, Any}(:lb => 0, :ub => 1,
                                        :mutually_consistent => false)
        r = compute_bounds(tp, bases, assumptions, dgp)
        results[k, 2:3] = [r[:lb], r[:ub]]

        assumptions[:mutually_consistent] = true
        r = compute_bounds(tp, bases, assumptions, dgp)
        results[k, 4:5] = [r[:lb], r[:ub]]

        assumptions[:decreasing_level] = [(2, 0), (2,1)]
        assumptions[:decreasing_difference] = [2]
        r = compute_bounds(tp, bases, assumptions, dgp)
        results[k, 6:7] = [r[:lb], r[:ub]]

        bases[2] = (interacted_bernstein_basis(1, notℓ = 1),
                    interacted_bernstein_basis(1, notℓ = 1))
        r = compute_bounds(tp, bases, assumptions, dgp)
        results[k, 8:9] = [r[:lb], r[:ub]]
    end

    return results
end
export simulation_prte

"""
    prte_misspecification()

This function produces the data to plot Figure 6 of MTW (2021).
It illustrates the misspecification caused for a PRTE by using the wrong choice
model.
It can't be used to create other/similar figures.
"""
function prte_misspecification(; ℓ_gen = 1)
    dgp = prte_dgp(ℓ_gen = ℓ_gen) # ℓ_gen shouldn't matter here (and doesn't)
    knots = vcat(0, dgp.pscore, 1)

    # Set up PRTE target parameter
    dgp_new = DGP(suppZ = dgp.suppZ,
                  pscore = dgp.pscore,
                  densZ = fill(.25, 4),
                  mtrs = dgp.mtrs)
    tp = prte_newz(dgp, dgp_new, ℓ = 1)

    nrows = 5
    results = DataFrame(name = fill("", nrows),
                        lb = fill(+Inf, nrows), ub = fill(+Inf, nrows))

    # True value
    truth = eval_tp(tp, [dgp.mtrs], dgp)
    results[1, :name] = "True value"
    results[1, :lb], results[1, :ub] = truth, truth

    # IAM model using both instruments
    # --> should get bounds that are non-positive
    bases = [(constantspline_basis(knots), constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(:lb => 0, :ub => 1)
    r = compute_bounds(tp, bases, assumptions, dgp)
    results[2, :name] = "Both instruments with IAM"
    results[2, :lb] = r[:lb]
    results[2, :ub] = r[:ub]

    # Use the first instrument only
    bases = [(interacted_constantspline_basis(knots, notℓ = 2),
              interacted_constantspline_basis(knots, notℓ = 2))]
    r = compute_bounds(tp, bases, assumptions, dgp)
    results[3, :name] = "Instrument 1 with IAM"
    results[3, :lb] = r[:lb]
    results[3, :ub] = r[:ub]

    # Use the second instrument only
    bases = [(interacted_constantspline_basis(knots, notℓ = 1),
              interacted_constantspline_basis(knots, notℓ = 1))]
    r = compute_bounds(tp, bases, assumptions, dgp)
    results[4, :name] = "Instrument 2 with IAM"
    results[4, :lb] = r[:lb]
    results[4, :ub] = r[:ub]

    # Use both instruments together and impose MC
    assumptions = Dict{Symbol, Any}(:lb => 0, :ub => 1,
                                    :mutually_consistent => true)
    bases = [(interacted_constantspline_basis(knots, notℓ = 2),
              interacted_constantspline_basis(knots, notℓ = 2)),
             (interacted_constantspline_basis(knots, notℓ = 1),
              interacted_constantspline_basis(knots, notℓ = 1))]
    r = compute_bounds(tp, bases, assumptions, dgp)
    results[5, :name] = "Both instruments with PM"
    results[5, :lb] = r[:lb]
    results[5, :ub] = r[:ub]
    return results
end
export prte_misspecification
