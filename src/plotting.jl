# Create coordinates for \addplot
# If `steps` is not specified, return a vector of coordinates, where
#   discontinuities separate vectors from one another.
#   This is particularly useful when drawing the MTRs.
# If `steps` is specified, return a vector of coordinates (i.e., endpoints)
#   `steps` controls how far apart the points are.
#   This is particularly useful when drawing tick marks for the weights.
# The MTRs and the weights use different strategies for creating their
# coordinates because we know the weights are piecewise-constant. Hence, we can
# draw them without approximation. On the other hand, we don't know where the
# MTRs are discontinuous. So, we discover where they have discontinuities by
# checking the slopes of the secant lines between two consecutive points. If
# the magnitude of the slopes is larger than `tol`, then we have found a
# discontinuity. These discontinuities determine the segments.
#
# TODO: is there a way to combine these two frameworks? Maybe using the basis
# for the MTRs can help determine whether we should "discover" the cutoffs...?
# TODO: use multiple dispatch to get rid of `isnothing(steps)` conditional.
function df_to_coordinates(df, xindex, yindex; steps = nothing, tol::Number = 7)
    x = round.(df[:, xindex], digits = 4)
    y = round.(df[:, yindex], digits = 4)
    coordinates = Vector{String}() # initialize empty vector of strings
    if isnothing(steps)
        lagdiff = v -> v - vcat(v[2:length(v)], 0)
        slope = lagdiff(y) ./ lagdiff(x)
        disc = findall(abs.(slope) .> tol)
        segment = ""
        for i in 1:nrow(df)
            segment = segment *
                "(" * string(x[i]) * "," * string(y[i]) * ")"
            if i in disc
                push!(coordinates, segment)
                segment = ""
            end
            if i == nrow(df)
                push!(coordinates, segment)
            end
        end
    else
        for segment_idx in 1:nrow(df)
            yval = y[segment_idx]
            if yval ≈ 0
                continue
            end
            endpoints = ""
            lb = x[segment_idx]
            if segment_idx == nrow(df)
                ub = 1
            else
                ub = x[segment_idx + 1]
            end
            grid = round.(range(lb, ub, step = steps), digits = 3)
            unique!(push!(grid, ub)) # ensure that ub is in grid
            for point_idx in 1:length(grid)
                endpoints = endpoints * "(" *
                    string(grid[point_idx]) * "," *
                    string(yval) * ")"
            end
            push!(coordinates, endpoints)
        end
    end
    return coordinates
end

# Generate the title used in the legend
# Q: should this be a property of the TargetParameter struct?
# If we can't specify default values of these properties, then it might break
# existing code.
function legendtitle(tp::TargetParameter)
    if tp.name == "LATE(u₁, u₂)"
        lb = tp.int_limits(1)[1]
        ub = tp.int_limits(1)[2]
        title = "LATE(\$ $(@sprintf("%.2f", lb)), $(@sprintf("%.2f", ub)) \$)"
    else
        @error "WIP" tp.name
    end
    return title
end
function legendtitle(ivlike::IVLike)
    title = ivlike.name
    if occursin("IV Slope for 𝟙(Z == z) for z ∈", ivlike.name)
        title = Vector{String}()
        for z in ivlike.params[:support]
            push!(title, "\$\\mathbb{1}[Z = $z]\$")
        end
    end
    if ivlike.name == "Saturated"
        title = Vector{String}()
        # BUG: the original Figure 5 in MST (2018) doesn't use the support of Z.
        # Instead, they use the indices of Z.
        # It is easier for me to use the indices of Z. If I want to correct
        # this mistake, I somehow need to pass information about the support to
        # this function.
        d_string = ["\$(1 - D)\$", "\$D\$"]
        z_string = "\$\\mathbb{1}[Z = " .* string.(1:(Int(length(ivlike.s) / 2))) .* "]\$"
        # NOTE: without [:], `title` would be a matrix
        title = [d * z for z in z_string, d in d_string][:]
    end
    return title
end

# Generate the title used in path for cross-referencing
function pathtitle(tp::TargetParameter)
    if tp.name == "LATE(u₁, u₂)"
        title = "late"
    else
        @error "WIP" tp.name
    end
    return title
end
function pathtitle(ivlike::IVLike)
    if ivlike.name == "IV Slope"
        title = "ivs"
    elseif ivlike.name == "OLS Slope"
        title = "olss"
    elseif occursin("IV Slope for 𝟙(Z == z) for z ∈", ivlike.name)
        title = "ivnps" .* string.(ivlike.params[:support])
    elseif ivlike.name == "Saturated"
        title = "saturated" .* string(collect(1:length(ivlike.s)))
    end
    return title
end

# Write bounds as an interval
function parse_bounds(result)
    lb = result[:lb]
    ub = result[:ub]
    return ": [\$ $(@sprintf("%.3f", lb)), $(@sprintf("%.3f", ub)) \$]"
end


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
    opts::Tuple{Dict, Vector{String}, Vector{String}, Vector{String}}
)
    # initialize
    settings, colors, marks, marksize = opts # aesthetic information
    m0segments = Vector{Dict}() # keeps track of MTR segments for d = 0
    m1segments = Vector{Dict}() # keeps track of MTR segments for d = 1
    d0weights = Vector{Dict}() # keeps track of weight segments for d = 0
    d1weights = Vector{Dict}() # keeps track of weight segments for d = 1
    legend = Vector{Dict}() # keeps track of legend entries
    aesthetic_counter = 1 # keeps track of colors, marks, and marksize

    # compute bounds
    result = compute_bounds(tp, basis, assumptions, dgp)

    # Collect data for MTRs
    if mtroption == "truth"
        mtr0 = dgp.mtrs[1]
        mtr1 = dgp.mtrs[2]
        settings[:mtrlegendtext] = "DGP MTRs"
    elseif mtroption == "max"
        mtr0 = result[:mtr_ub][1][1]
        mtr1 = result[:mtr_ub][1][2]
        settings[:title] = settings[:title] * parse_bounds(result)
        settings[:mtrlegendtext] = "Maximizing MTRs"
    elseif mtroption == "min"
        mtr0 = result[:mtr_lb][1][1]
        mtr1 = result[:mtr_lb][1][2]
        settings[:title] = settings[:title] * parse_bounds(result)
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
    tp_weights = compute_average_weights(tp)
    tp_d0_coord = df_to_coordinates(tp_weights, :u, 3, steps = 1/step)
    tp_d1_coord = df_to_coordinates(tp_weights, :u, 2, steps = 1/step)
    push!(legend, Dict(
        "color" => colors[aesthetic_counter],
        "mark" => marks[aesthetic_counter],
        "marksize" => marksize[aesthetic_counter],
        "legendtitle" => legendtitle(tp)
    ))
    for coordinate_idx in 1:length(tp_d0_coord)
        segment = Dict(
            "pathname" => "d0" * pathtitle(tp) * string(coordinate_idx),
            "color" => colors[aesthetic_counter],
            "mark" => marks[aesthetic_counter],
            "marksize" => marksize[aesthetic_counter],
            "coordinates" => tp_d0_coord[coordinate_idx]
        )
        push!(d0weights, segment)
    end
    for coordinate_idx in 1:length(tp_d1_coord)
        segment = Dict(
            "pathname" => "d1" * pathtitle(tp) * string(coordinate_idx),
            "color" => colors[aesthetic_counter],
            "mark" => marks[aesthetic_counter],
            "marksize" => marksize[aesthetic_counter],
            "coordinates" => tp_d1_coord[coordinate_idx]
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
            "legendtitle" => legendtitle(s)
        ))

        # Store data in dictionary for Mustache.jl for d = 0 weights
        for coordinate_idx in 1:length(ivlike_d0_coord)
            segment = Dict(
                "pathname" => "d0" * pathtitle(s) * string(coordinate_idx),
                "color" => colors[aesthetic_counter],
                "mark" => marks[aesthetic_counter],
                "marksize" => marksize[aesthetic_counter],
                "coordinates" => ivlike_d0_coord[coordinate_idx]
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
                "coordinates" => ivlike_d1_coord[coordinate_idx]
            )
            push!(d1weights, segment)
        end
        aesthetic_counter += 1
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
            "legendtitle" => legendtitle(s)
        ))

        # Store data in dictionary for Mustache.jl for d = 0 weights
        for coordinate_idx in 1:length(ivlike_d0_coord)
            segment = Dict(
                "pathname" => "d0" * pathtitle(s) * string(coordinate_idx),
                "color" => colors[aesthetic_counter],
                "mark" => marks[aesthetic_counter],
                "marksize" => marksize[aesthetic_counter],
                "coordinates" => ivlike_d0_coord[coordinate_idx]
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
                "coordinates" => ivlike_d1_coord[coordinate_idx]
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
                "legendtitle" => legendtitle(s_ivlike)[support_idx]
            ))

            # Store data in dictionary for Mustache.jl for d = 0 weights
            for coordinate_idx in 1:length(ivlike_d0_coord)
                segment = Dict(
                    "pathname" => "d0" * pathtitle(s_ivlike)[support_idx] * string(coordinate_idx),
                    "color" => colors[aesthetic_counter],
                    "mark" => marks[aesthetic_counter],
                    "marksize" => marksize[aesthetic_counter],
                    "coordinates" => ivlike_d0_coord[coordinate_idx]
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
                    "coordinates" => ivlike_d1_coord[coordinate_idx]
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
                "legendtitle" => legendtitle(s_ivlike)[saturated_idx]
            ))

            # Store data in dictionary for Mustache.jl for d = 0 weights
            for coordinate_idx in 1:length(ivlike_d0_coord)
                segment = Dict(
                    "pathname" => "d0" * pathtitle(s_ivlike)[saturated_idx] * string(coordinate_idx),
                    "color" => colors[aesthetic_counter],
                    "mark" => marks[aesthetic_counter],
                    "marksize" => marksize[aesthetic_counter],
                    "coordinates" => ivlike_d0_coord[coordinate_idx]
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
                    "coordinates" => ivlike_d1_coord[coordinate_idx]
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
    templatefn = joinpath(savedir, "mst2018econometrica", "tikz-template.tex")
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
    return texfn
end
