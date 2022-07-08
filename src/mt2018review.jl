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

# Figure 1: MTRs and MTE
function run_tikz_mtr(savedir::String, compile::Bool = false)
    dgp = dgp_review()
    texfn = mtr_mte(savedir, "tikz-mtr";
        dgp = dgp
    )
    if compile
        compile_latex(texfn)
    end
end

# Figure 2: weights for conventional target parameters
function run_tikz_weights(savedir::String, compile::Bool = false)
    dgp = dgp_review()
    texfn = conventional_weights(savedir, "tikz-weights";
        dgp = dgp
    )
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
    templatefn = joinpath(savedir, "mt2018review", "tikz-mtr-template.tex")
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
    # setup
    weights = Vector{Dict}()
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
            "opts" => ifelse(segment_idx == length(segment_info), "", "forget plot")
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
            "opts" => ifelse(segment_idx == length(segment_info), "", "forget plot")
        ))
    end
    tp_atu = Dict(
        "legendtitle" => "\$\\text{ATu}\$",
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
