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
        :ylabelweights => nothing
    )

    colors = ["gray", "darkgray", "lightgray", "darkgray", "lightgray",
              "darkgray", "lightgray"]
    marks = ["*", "x", "diamond*", "square*", "|", "triangle*", "pentagon*"]
    marksize = ["1pt", "2.5pt", "1.5pt", "1.25pt", "3pt", "1.25pt", "1.25pt"]

    return settings, colors, marks, marksize
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
    # TODO: separate mtr data from aesthetic information in settings

    # initialize
    settings, colors, marks, marksize = opts # aesthetic information
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

    println("...store MTR data in dictionary") # DEBUG:

    # Store MTR data in dictionary for Mustache.jl
    settings[:m0coordinates] = df_to_coordinates(mtr_results, :u, :mtr0)
    settings[:m1coordinates] = df_to_coordinates(mtr_results, :u, :mtr1)
    settings[:ylabelweights] = "Weights (where \$\\neq 0\$)"


    # Collect data for target parameter
    tp_weights = compute_average_weights(tp)
    tp_d0_coord = df_to_coordinates(tp_weights, :u, 3, steps = 1/500)
    tp_d1_coord = df_to_coordinates(tp_weights, :u, 2, steps = 1/500)
    println("...collect data for tp") # DEBUG:
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

    println("...collect data for ivlike") # DEBUG:

    # Collect data for IV-like estimands
    ivlike_weights = Vector() # used to compute max and min weights
    if haskey(assumptions, :ivslope) && assumptions[:ivslope]
        println("...collect data for ivslope") # DEBUG:
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

        println("...collect data for ivslope, d = 0") # DEBUG:

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

        println("...collect data for ivslope, d = 1") # DEBUG:

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
        println("...collect data for olsslope") # DEBUG:
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

        println("...collect data for olsslope, d = 0") # DEBUG:

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

        println("...collect data for olsslope, d = 1") # DEBUG:

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
            println("ivslope ind: ", zind, " Z = ", dgp.suppZ[zind])
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

            println("...collect data for ivslope, d = 0, z = " * string(dgp.suppZ[zind])) # DEBUG:

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

            println("...collect data for ivslope, d = 1, z = " * string(dgp.suppZ[zind])) # DEBUG:

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
            println("...collect data for sharp: ", saturated_idx) # DEBUG:
            ivlike_d0_coord = Vector()
            ivlike_d1_coord = Vector()
            s = IVLike(
                "Saturated; index $saturated_idx",
                [s_ivlike.s[legend_order[saturated_idx]]],
                nothing
            )
            println("...compute average weights") # DEBUG:
            s_weights = compute_average_weights(s, dgp) # BOOKMARK: why are all the weights equal to 0?
            println("...finish computing average weights") # DEBUG:
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
            println(colors[aesthetic_counter]) # DEBUG:

            println("...collect data for saturated, d = 0") # DEBUG:

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

            println("...collect data for saturated, d = 1") # DEBUG:

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

    println("...update aesthetic info") # DEBUG:

    # Update aesthetic information based on weights
    settings[:weightymax] = ceil(max(
        tp_weights[:, 2]...,
        tp_weights[:, 3]...,
        ivlike_weights...
    )) + 1
    println("...done updating aesthic info") # DEBUG
    settings[:weightymin] = -1 * settings[:weightymax]

    println("...create tex file") # DEBUG:

    # Create tex file
    templatefn = joinpath(savedir, "mst2018econometrica", "tikz-template.tex")
    template = Mustache.load(templatefn, ("<<", ">>"))
    tex = render(
        template;
        settings...,
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
export mtrs_and_weights
