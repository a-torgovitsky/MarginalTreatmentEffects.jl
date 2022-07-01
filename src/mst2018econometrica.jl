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
    Dict(
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
end

# This function returns the tex file for producing Figure 1.
function dgp_to_tex(savedir::String)
    dgp = dgp_econometrica()
    tp = late(dgp, 0.35, 0.9)
    ivlike = ivslope(dgp)
    settings = defaults_econometrica()
    d0weights = Vector{Dict}()
    d1weights = Vector{Dict}()
    legend = Vector{Dict}()

    # Collect MTR data
    ev = DataFrame(z = 1, u = 0:0.1:1)
    mtr0 = evaluate_mtr(dgp.mtrs[1], ev)
    mtr1 = evaluate_mtr(dgp.mtrs[2], ev)
    mtr_results = DataFrame(u = 0:0.1:1)
    mtr_results[:, "mtr0"] = mtr0
    mtr_results[:, "mtr1"] = mtr1

    # Store data in dictionary for Mustache.jl
    settings[:m0coordinates] = df_to_coordinates(mtr_results, :u, :mtr0)
    settings[:m1coordinates] = df_to_coordinates(mtr_results, :u, :mtr1)
    settings[:title] = "~"
    settings[:mtrlegendtext] = "DGP MTRs"
    settings[:ylabelweights] = "Weights (where \$\\neq 0\$)"

    # Collect LATE(0.35, 0.90) data
    tp_weights = compute_average_weights(tp)
    tp_d0_coord = df_to_coordinates(tp_weights, :u, 3, steps = 1/500)
    tp_d1_coord = df_to_coordinates(tp_weights, :u, 2, steps = 1/500)

    # Collect IV Slope data
    ivlike_weights = compute_average_weights(ivlike, dgp)
    ivlike_d0_coord = df_to_coordinates(ivlike_weights, :u, 3, steps = 1/500)
    ivlike_d1_coord = df_to_coordinates(ivlike_weights, :u, 2, steps = 1/500)

    # Store data in dictionary for Mustache.jl for d = 0 weights
    d0weights = Vector{Dict}()
    for coordinate_idx in 1:length(tp_d0_coord)
        segment = Dict(
            "pathname" => "d0late" * string(coordinate_idx),
            "color" => "gray",
            "mark" => "*",
            "marksize" => "1pt",
            "coordinates" => tp_d0_coord[coordinate_idx]
        )
        push!(d0weights, segment)
    end
    push!(legend, Dict(
        "color" => "gray",
        "mark" => "*",
        "marksize" => "1pt",
        "legendtitle" => "LATE(\$0.35, 0.90\$)"
    ))
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
    push!(legend, Dict(
        "color" => "darkgray",
        "mark" => "x",
        "marksize" => "2.5pt",
        "legendtitle" => "IV Slope"
    ))

    # Store data in dictionary for Mustache.jl for d = 1 weights
    d1weights = Vector{Dict}()
    for coordinate_idx in 1:length(tp_d1_coord)
        segment = Dict(
            "pathname" => "d1late" * string(coordinate_idx),
            "color" => "gray",
            "mark" => "*",
            "marksize" => "1pt",
            "coordinates" => tp_d1_coord[coordinate_idx]
        )
        push!(d1weights, segment)
    end
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

    # Update aesthetic information based on weights
    settings[:weightymax] = ceil(max(
        tp_weights[:, 2]...,
        tp_weights[:, 3]...,
        ivlike_weights[:, 2]...,
        ivlike_weights[:, 3]...
    )) + 1
    settings[:weightymin] = -1 * settings[:weightymax]

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
    texfn = joinpath(dirname(templatefn), "dgp.tex")
    open(texfn, "w") do file
        write(file, tex)
    end
    return texfn
end
export dgp_to_tex
