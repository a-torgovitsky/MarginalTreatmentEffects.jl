# TODO: update documentation
"""
    mtrs_and_weights(
        savedir::String,
        filename::String;
        dgp::DGP = dgp_econometrica(),
        tp::TargetParameter = late(dgp, 0.35, 0.9),
        assumptions::Dict,
        mtroption::String,
        defaults = defaults_econometrica()
    )

This function produces the tex file used to create the figures in MST (2018).
The file path of this tex file will be returned.

# Parameters

- savedir: directory where the tex file will be stored
- filename: name of the tex file
- dgp: data generating process
- tp: target parameter
- basis: vector of tuples, where each tuple is a pair of MTRBasis structs
- assumptions: dictionary of assumptions, including IV-like estimand
- mtroption: either "truth", "max", or "min"
- opts: tuple of four elements:
    1. settings::Dict
    2. colors::Vector{String}
    3. marks::Vector{String}
    4. marksize::Vector{String}
"""
function mtrs_and_weights(
    savedir::String,
    filename::String;
    dgp::DGP,
    tp::TargetParameter,
    basis::Vector{Tuple{MTRBasis, MTRBasis}},
    assumptions::Dict,
    mtroption::String,
    opts::Tuple{Dict, Vector{String}, Vector{String}, Vector{String}, Vector{String}},
    attributes::Dict = Dict("LogLevel" => 0),
    startdf::DataFrame= DataFrame(
        ℓ = Int64[],
        d = Int64[],
        j = Int64[],
        k = Int64[],
        start = Float64[]
    ),
    fixdf::DataFrame= DataFrame(
        ℓ = Int64[],
        d = Int64[],
        j = Int64[],
        k = Int64[],
        fix = Float64[]
    ),
    return_bounds::Bool = false
)
    # initialize
    settings, colors, marks, marksize, linetype = opts # aesthetic information
    m0segments = Vector{Dict}() # keeps track of MTR segments for d = 0
    m1segments = Vector{Dict}() # keeps track of MTR segments for d = 1
    d0weights = Vector{Dict}() # keeps track of weight segments for d = 0
    d1weights = Vector{Dict}() # keeps track of weight segments for d = 1
    legend = Vector{Dict}() # keeps track of legend entries
    aesthetic_counter = 1 # keeps track of colors, marks, and marksize

    # compute bounds
    result = compute_bounds(
        tp,
        basis,
        assumptions,
        dgp,
        attributes,
        startdf,
        fixdf
    )

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
        @error "unsupported" mtroption
    end
    step = 500
    ugrid = (1/step):(1/step):(1 - 1/step)
    ev = DataFrame(z = 1, u = ugrid)
    mtr_results = DataFrame(u = (1/step):(1/step):(1 - 1/step))
    mtr_results[:, "mtr0"] = evaluate_mtr(mtr0, ev)
    mtr_results[:, "mtr1"] = evaluate_mtr(mtr1, ev)

    # Store MTR data in dictionary for Mustache.jl
    # TODO: avoid duplication b/t the for-loops for mtr0 and mtr1
    m0coordinates = df_to_coordinates(mtr_results, :u, :mtr0)
    for coordinate_idx in 1:length(m0coordinates)
        segment = Dict(
            "pathname" => "mtr0" * string(coordinate_idx),
            "coordinates" => m0coordinates[coordinate_idx]
        )
        push!(m0segments, segment)
    end
    m1coordinates = df_to_coordinates(mtr_results, :u, :mtr1)
    for coordinate_idx in 1:length(m1coordinates)
        segment = Dict(
            "pathname" => "mtr1" * string(coordinate_idx),
            "coordinates" => m1coordinates[coordinate_idx]
        )
        push!(m1segments, segment)
    end

    # Collect data for target parameter
    tp_weights = compute_average_weights(tp, dgp)
    tp_d0_coord = df_to_coordinates(tp_weights, :u, 3, steps = 1/step)
    tp_d1_coord = df_to_coordinates(tp_weights, :u, 2, steps = 1/step)
    push!(legend, Dict(
        "color" => colors[aesthetic_counter],
        "mark" => marks[aesthetic_counter],
        "marksize" => marksize[aesthetic_counter],
        "legendtitle" => legendtitle(tp),
        "linetype" => linetype[aesthetic_counter]
    ))
    for coordinate_idx in 1:length(tp_d0_coord)
        segment = Dict(
            "pathname" => "d0" * pathtitle(tp) * string(coordinate_idx),
            "color" => colors[aesthetic_counter],
            "mark" => marks[aesthetic_counter],
            "marksize" => marksize[aesthetic_counter],
            "coordinates" => tp_d0_coord[coordinate_idx],
            "linetype" => linetype[aesthetic_counter]
        )
        push!(d0weights, segment)
    end
    for coordinate_idx in 1:length(tp_d1_coord)
        segment = Dict(
            "pathname" => "d1" * pathtitle(tp) * string(coordinate_idx),
            "color" => colors[aesthetic_counter],
            "mark" => marks[aesthetic_counter],
            "marksize" => marksize[aesthetic_counter],
            "coordinates" => tp_d1_coord[coordinate_idx],
            "linetype" => linetype[aesthetic_counter]
        )
        push!(d1weights, segment)
    end
    aesthetic_counter += 1

    # Collect data for IV-like estimands
    # TODO: within each if-block, the code is quite similar across the IV-like
    # etimands. Can we avoid the duplication?
    ivlike_weights = Vector() # used to compute max and min weights

    if haskey(assumptions, :ivslope) && assumptions[:ivslope]
        ivlike_d0_coord = Vector()
        ivlike_d1_coord = Vector()
        s = ivslope(dgp)
        s_weights = compute_average_weights(s, dgp)
        push!(ivlike_weights, s_weights[:, 2]...)
        push!(ivlike_weights, s_weights[:, 3]...)
        d0_coordinates = df_to_coordinates(s_weights, :u, 3, steps = 1/500)
        push!(ivlike_d0_coord, d0_coordinates...)
        d1_coordinates = df_to_coordinates(s_weights, :u, 2, steps = 1/500)
        push!(ivlike_d1_coord, d1_coordinates...)
        push!(legend, Dict(
            "color" => colors[aesthetic_counter],
            "mark" => marks[aesthetic_counter],
            "marksize" => marksize[aesthetic_counter],
            "legendtitle" => legendtitle(s),
            "linetype" => linetype[aesthetic_counter]
        ))

        # Store data in dictionary for Mustache.jl for d = 0 weights
        for coordinate_idx in 1:length(ivlike_d0_coord)
            segment = Dict(
                "pathname" => "d0" * pathtitle(s) * string(coordinate_idx),
                "color" => colors[aesthetic_counter],
                "mark" => marks[aesthetic_counter],
                "marksize" => marksize[aesthetic_counter],
                "coordinates" => ivlike_d0_coord[coordinate_idx],
                "linetype" => linetype[aesthetic_counter]
            )
            push!(d0weights, segment)
        end

        # Store data in dictionary for Mustache.jl for d = 1 weights
        for coordinate_idx in 1:length(ivlike_d1_coord)
            segment = Dict(
                "pathname" => "d1" * pathtitle(s) * string(coordinate_idx),
                "color" => colors[aesthetic_counter],
                "mark" => marks[aesthetic_counter],
                "marksize" => marksize[aesthetic_counter],
                "coordinates" => ivlike_d1_coord[coordinate_idx],
                "linetype" => linetype[aesthetic_counter]
            )
            push!(d1weights, segment)
        end
        aesthetic_counter += 1
    end

    if haskey(assumptions, :tslsslopeind) && assumptions[:tslsslopeind]
        ivlike_d0_coord = Vector()
        ivlike_d1_coord = Vector()
        s = tslsslope_indicator(dgp)
        s_weights = compute_average_weights(s, dgp)
        push!(ivlike_weights, s_weights[:, 2]...)
        push!(ivlike_weights, s_weights[:, 3]...)
        d0_coordinates = df_to_coordinates(s_weights, :u, 3, steps = 1/500)
        push!(ivlike_d0_coord, d0_coordinates...)
        d1_coordinates = df_to_coordinates(s_weights, :u, 2, steps = 1/500)
        push!(ivlike_d1_coord, d1_coordinates...)
        push!(legend, Dict(
            "color" => colors[aesthetic_counter],
            "mark" => marks[aesthetic_counter],
            "marksize" => marksize[aesthetic_counter],
            "legendtitle" => legendtitle(s),
            "linetype" => linetype[aesthetic_counter]
        ))

        # Store data in dictionary for Mustache.jl for d = 0 weights
        for coordinate_idx in 1:length(ivlike_d0_coord)
            segment = Dict(
                "pathname" => "d0" * pathtitle(s) * string(coordinate_idx),
                "color" => colors[aesthetic_counter],
                "mark" => marks[aesthetic_counter],
                "marksize" => marksize[aesthetic_counter],
                "coordinates" => ivlike_d0_coord[coordinate_idx],
                "linetype" => linetype[aesthetic_counter]
            )
            push!(d0weights, segment)
        end

        # Store data in dictionary for Mustache.jl for d = 1 weights
        for coordinate_idx in 1:length(ivlike_d1_coord)
            segment = Dict(
                "pathname" => "d1" * pathtitle(s) * string(coordinate_idx),
                "color" => colors[aesthetic_counter],
                "mark" => marks[aesthetic_counter],
                "marksize" => marksize[aesthetic_counter],
                "coordinates" => ivlike_d1_coord[coordinate_idx],
                "linetype" => linetype[aesthetic_counter]
            )
            push!(d1weights, segment)
        end
        aesthetic_counter += 1
    end

    if haskey(assumptions, :wald)
        for p in assumptions[:wald]
            ivlike_d0_coord = Vector()
            ivlike_d1_coord = Vector()
            s = wald(dgp; z₀ = p[1], z₁ = p[2])
            s_weights = compute_average_weights(s, dgp)
            push!(ivlike_weights, s_weights[:, 2]...)
            push!(ivlike_weights, s_weights[:, 3]...)
            d0_coordinates = df_to_coordinates(s_weights, :u, 3, steps = 1/500)
            push!(ivlike_d0_coord, d0_coordinates...)
            d1_coordinates = df_to_coordinates(s_weights, :u, 2, steps = 1/500)
            push!(ivlike_d1_coord, d1_coordinates...)
            push!(legend, Dict(
                "color" => colors[aesthetic_counter],
                "mark" => marks[aesthetic_counter],
                "marksize" => marksize[aesthetic_counter],
                "legendtitle" => legendtitle(s),
                "linetype" => linetype[aesthetic_counter]
            ))

            # Store data in dictionary for Mustache.jl for d = 0 weights
            for coordinate_idx in 1:length(ivlike_d0_coord)
                segment = Dict(
                    "pathname" => "d0" * pathtitle(s) * string(coordinate_idx),
                    "color" => colors[aesthetic_counter],
                    "mark" => marks[aesthetic_counter],
                    "marksize" => marksize[aesthetic_counter],
                    "coordinates" => ivlike_d0_coord[coordinate_idx],
                    "linetype" => linetype[aesthetic_counter]
                )
                push!(d0weights, segment)
            end

            # Store data in dictionary for Mustache.jl for d = 1 weights
            for coordinate_idx in 1:length(ivlike_d1_coord)
                segment = Dict(
                    "pathname" => "d1" * pathtitle(s) * string(coordinate_idx),
                    "color" => colors[aesthetic_counter],
                    "mark" => marks[aesthetic_counter],
                    "marksize" => marksize[aesthetic_counter],
                    "coordinates" => ivlike_d1_coord[coordinate_idx],
                    "linetype" => linetype[aesthetic_counter]
                )
                push!(d1weights, segment)
            end
            aesthetic_counter += 1
        end
    end

    if haskey(assumptions, :olsslope) && assumptions[:olsslope]
        ivlike_d0_coord = Vector()
        ivlike_d1_coord = Vector()
        s = olsslope(dgp)
        s_weights = compute_average_weights(s, dgp)
        push!(ivlike_weights, s_weights[:, 2]...)
        push!(ivlike_weights, s_weights[:, 3]...)
        d0_coordinates = df_to_coordinates(s_weights, :u, 3, steps = 1/500)
        push!(ivlike_d0_coord, d0_coordinates...)
        d1_coordinates = df_to_coordinates(s_weights, :u, 2, steps = 1/500)
        push!(ivlike_d1_coord, d1_coordinates...)
        push!(legend, Dict(
            "color" => colors[aesthetic_counter],
            "mark" => marks[aesthetic_counter],
            "marksize" => marksize[aesthetic_counter],
            "legendtitle" => legendtitle(s),
            "linetype" => linetype[aesthetic_counter]
        ))

        # Store data in dictionary for Mustache.jl for d = 0 weights
        for coordinate_idx in 1:length(ivlike_d0_coord)
            segment = Dict(
                "pathname" => "d0" * pathtitle(s) * string(coordinate_idx),
                "color" => colors[aesthetic_counter],
                "mark" => marks[aesthetic_counter],
                "marksize" => marksize[aesthetic_counter],
                "coordinates" => ivlike_d0_coord[coordinate_idx],
                "linetype" => linetype[aesthetic_counter]
            )
            push!(d0weights, segment)
        end

        # Store data in dictionary for Mustache.jl for d = 1 weights
        for coordinate_idx in 1:length(ivlike_d1_coord)
            segment = Dict(
                "pathname" => "d1" * pathtitle(s) * string(coordinate_idx),
                "color" => colors[aesthetic_counter],
                "mark" => marks[aesthetic_counter],
                "marksize" => marksize[aesthetic_counter],
                "coordinates" => ivlike_d1_coord[coordinate_idx],
                "linetype" => linetype[aesthetic_counter]
            )
            push!(d1weights, segment)
        end
        aesthetic_counter += 1
    end

    if haskey(assumptions, :ivslopeind)
        s_ivlike = ivslope_indicator(dgp; support = assumptions[:ivslopeind])
        for support_idx in 1:length(s_ivlike.params[:support])
            zind = s_ivlike.params[:support][support_idx]
            ivlike_d0_coord = Vector()
            ivlike_d1_coord = Vector()
            s = IVLike(
                "IV Slope for 𝟙(Z = " * string(dgp.suppZ[zind]) * ")",
                [s_ivlike.s[support_idx]],
                Dict(:support => [dgp.suppZ[zind]])
            )
            s_weights = compute_average_weights(s, dgp)
            push!(ivlike_weights, s_weights[:, 2]...)
            push!(ivlike_weights, s_weights[:, 3]...)
            d0_coordinates = df_to_coordinates(s_weights, :u, 3, steps = 1/step)
            push!(ivlike_d0_coord, d0_coordinates...)
            d1_coordinates = df_to_coordinates(s_weights, :u, 2, steps = 1/step)
            push!(ivlike_d1_coord, d1_coordinates...)
            push!(legend, Dict(
                "color" => colors[aesthetic_counter],
                "mark" => marks[aesthetic_counter],
                "marksize" => marksize[aesthetic_counter],
                "legendtitle" => legendtitle(s_ivlike)[support_idx],
                "linetype" => linetype[aesthetic_counter]
            ))

            # Store data in dictionary for Mustache.jl for d = 0 weights
            for coordinate_idx in 1:length(ivlike_d0_coord)
                segment = Dict(
                    "pathname" => "d0" * pathtitle(s_ivlike)[support_idx] * string(coordinate_idx),
                    "color" => colors[aesthetic_counter],
                    "mark" => marks[aesthetic_counter],
                    "marksize" => marksize[aesthetic_counter],
                    "coordinates" => ivlike_d0_coord[coordinate_idx],
                    "linetype" => linetype[aesthetic_counter]
                )
                push!(d0weights, segment)
            end

            # Store data in dictionary for Mustache.jl for d = 1 weights
            for coordinate_idx in 1:length(ivlike_d1_coord)
                segment = Dict(
                    "pathname" => "d1" * pathtitle(s_ivlike)[support_idx] * string(coordinate_idx),
                    "color" => colors[aesthetic_counter],
                    "mark" => marks[aesthetic_counter],
                    "marksize" => marksize[aesthetic_counter],
                    "coordinates" => ivlike_d1_coord[coordinate_idx],
                    "linetype" => linetype[aesthetic_counter]
                )
                push!(d1weights, segment)
            end

            aesthetic_counter += 1
        end
    end

    if haskey(assumptions, :saturated) && assumptions[:saturated]
        s_ivlike = make_slist(dgp.suppZ)
        legend_order = [1, 3, 5, 2, 4, 6] # ensures correct legend aesthetics
        for saturated_idx in 1:length(s_ivlike.s)
            ivlike_d0_coord = Vector()
            ivlike_d1_coord = Vector()
            s = IVLike(
                "Saturated; index $saturated_idx",
                [s_ivlike.s[legend_order[saturated_idx]]],
                nothing
            )
            s_weights = compute_average_weights(s, dgp) # BOOKMARK: why are all the weights equal to 0?
            push!(ivlike_weights, s_weights[:, 2]...)
            push!(ivlike_weights, s_weights[:, 3]...)
            d0_coordinates = df_to_coordinates(s_weights, :u, 3, steps = 1/500)
            push!(ivlike_d0_coord, d0_coordinates...)
            d1_coordinates = df_to_coordinates(s_weights, :u, 2, steps = 1/500)
            push!(ivlike_d1_coord, d1_coordinates...)
            push!(legend, Dict(
                "color" => colors[aesthetic_counter],
                "mark" => marks[aesthetic_counter],
                "marksize" => marksize[aesthetic_counter],
                "legendtitle" => legendtitle(s_ivlike)[saturated_idx],
                "linetype" => linetype[aesthetic_counter]
            ))

            # Store data in dictionary for Mustache.jl for d = 0 weights
            for coordinate_idx in 1:length(ivlike_d0_coord)
                segment = Dict(
                    "pathname" => "d0" * pathtitle(s_ivlike)[saturated_idx] * string(coordinate_idx),
                    "color" => colors[aesthetic_counter],
                    "mark" => marks[aesthetic_counter],
                    "marksize" => marksize[aesthetic_counter],
                    "coordinates" => ivlike_d0_coord[coordinate_idx],
                    "linetype" => linetype[aesthetic_counter]
                )
                push!(d0weights, segment)
            end

            # Store data in dictionary for Mustache.jl for d = 1 weights
            for coordinate_idx in 1:length(ivlike_d1_coord)
                segment = Dict(
                    "pathname" => "d1" * pathtitle(s_ivlike)[saturated_idx] * string(coordinate_idx),
                    "color" => colors[aesthetic_counter],
                    "mark" => marks[aesthetic_counter],
                    "marksize" => marksize[aesthetic_counter],
                    "coordinates" => ivlike_d1_coord[coordinate_idx],
                    "linetype" => linetype[aesthetic_counter]
                )
                push!(d1weights, segment)
            end
            aesthetic_counter += 1
        end
    end

    # Update aesthetic information based on weights
    settings[:weightymax] = ceil(max(
        tp_weights[:, 2]...,
        tp_weights[:, 3]...,
        ivlike_weights...
    )) + 1
    settings[:weightymin] = -1 * settings[:weightymax]

    # Create tex file
    # `project` is defined in global scope
    templatefn = joinpath(savedir, project, "tikz-template.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(
        template;
        settings...,
        m0segments = m0segments,
        m1segments = m1segments,
        d0weights = d0weights,
        d1weights = d1weights,
        legend = legend
    )
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
        xlabel = xlabel,
        cycleshift = 0
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
        xlabel = xlabel,
        cycleshift = 2
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
        xlabel = xlabel,
        cycleshift = 0
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
