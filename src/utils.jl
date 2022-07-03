function make_dirname(savelocation::String;
                      stub::String="unnamed",
                      tag::String="")
    println("Saving in: $savelocation")
    if isempty(tag)
        print("Enter a directory tag "*
              "(or hit Enter to generate from the date): ")
        tag = readline()
        if isempty(tag)
            tag = Libc.strftime("%d-%b-%y-at-%H-%M-%S", time())
        end
    end

    if isempty(stub)
        dirname = tag
    else
        dirname = stub*"-"*tag
    end
    pathstr = joinpath(savelocation, dirname)
    return pathstr
end

function make_dir(pathstr::String)
    if isdir(pathstr)
        @info "Directory already exists, so not copying source files."
    else
        mkpath(pathstr)
        cplist = ["LICENSE", "Project.toml", "README.md",
                  "src", "tex"]
        sourcepathstr = joinpath(@__DIR__, "../")
        for c in cplist
            cp(joinpath(sourcepathstr, c),
               joinpath(pathstr, c),
               follow_symlinks=true)
        end
    end

    # Make a directory for holding output
    resultsdir = joinpath(pathstr, "results")
    already_existed = false
    if isdir(resultsdir)
        @info "Results directory already exists, not modifying."
        already_existed = true
    else
        mkdir(resultsdir)
    end

    # Copy .tex files to results directory
    texdir = joinpath(resultsdir, "../tex")
    cplist = readdir(texdir)
    for c in cplist
        cp(joinpath(texdir, c), joinpath(resultsdir, c))
    end

    return resultsdir, already_existed
end

function setup(savelocation::String;
               stub::String="unnamed",
               tag::String="")
    pathstr = make_dirname(savelocation, stub = stub, tag = tag)
    resultsdir, already_existed = make_dir(pathstr)
    return resultsdir, already_existed
end

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
