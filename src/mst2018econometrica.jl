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
- assumptions: dictionary of assumptions, including IV-like estimand
- mtroption: either "truth", "max" or "min"
- defaults: tuple of four elements:
    1. settings::Dict
    2. colors::Vector{String}
    3. marks::Vector{String}
    4. marksize::Vector{String}
"""
function mtrs_and_weights(
    savedir::String,
    filename::String;
    dgp::DGP = dgp_econometrica(),
    tp::TargetParameter = late(dgp, 0.35, 0.9),
    assumptions::Dict,
    mtroption::String,
    defaults = defaults_econometrica()
)
    # TODO: separate mtr data from aesthetic information in settings

    # initialize
    settings, colors, marks, marksize = defaults # aesthetic information
    d0weights = Vector{Dict}() # keeps track of weight segments for d = 0
    d1weights = Vector{Dict}() # keeps track of weight segments for d = 1
    legend = Vector{Dict}() # keeps track of legend entries
    aesthetic_counter = 1 # keeps track of colors, marks, and marksize

    # Collect data for MTRs
    ev = DataFrame(z = 1, u = 0:0.1:1)
    if mtroption == "truth"
        mtr0 = evaluate_mtr(dgp.mtrs[1], ev)
        mtr1 = evaluate_mtr(dgp.mtrs[2], ev)
        settings[:title] = "~"
    else
        @error "unsupported" mtroption
    end
    mtr_results = DataFrame(u = 0:0.1:1)
    mtr_results[:, "mtr0"] = mtr0
    mtr_results[:, "mtr1"] = mtr1

    println("...store MTR data in dictionary") # DEBUG:

    # Store MTR data in dictionary for Mustache.jl
    settings[:m0coordinates] = df_to_coordinates(mtr_results, :u, :mtr0)
    settings[:m1coordinates] = df_to_coordinates(mtr_results, :u, :mtr1)
    settings[:mtrlegendtext] = "DGP MTRs"
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
    println("...finished") # DEBUG:
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
    if haskey(assumptions, :ivslope)
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
                "pathname" => "d0ivslope" * string(coordinate_idx),
                "color" => "darkgray",
                "mark" => "x",
                "marksize" => "2.5pt",
                "coordinates" => ivlike_d0_coord[coordinate_idx]
            )
            push!(d0weights, segment)
        end

        println("...collect data for ivslope, d = 1") # DEBUG:

        # Store data in dictionary for Mustache.jl for d = 1 weights
        for coordinate_idx in 1:length(ivlike_d1_coord)
            segment = Dict(
                "pathname" => "d1ivslope" * string(coordinate_idx),
                "color" => "darkgray",
                "mark" => "x",
                "marksize" => "2.5pt",
                "coordinates" => ivlike_d1_coord[coordinate_idx]
            )
            push!(d1weights, segment)
        end
        aesthetic_counter += 1
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
