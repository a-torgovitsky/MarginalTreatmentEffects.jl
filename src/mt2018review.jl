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
