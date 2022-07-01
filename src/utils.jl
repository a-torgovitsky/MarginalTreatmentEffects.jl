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
# If steps is not specified, return a string of coordinates to draw curves
# If steps is specified, return a vector of coordinates (i.e., endpoints)
# steps controls how far apart the points are.
# This is particularly useful when drawing tick marks.
function df_to_coordinates(df, xindex, yindex; steps = nothing)
    x = round.(df[:, xindex], digits = 3)
    y = round.(df[:, yindex], digits = 3)
    if isnothing(steps)
        coordinates = ""
        for i in 1:nrow(df)
            coordinates = coordinates *
                "(" * string(x[i]) * "," * string(y[i]) * ")"
        end
    else
        coordinates = Vector{String}() # initialize empty vector of strings
        for segment_idx in 1:(nrow(df) - 1)
            endpoints = ""
            lb = x[segment_idx]
            ub = x[segment_idx + 1]
            yval = y[segment_idx]
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
    return ivlike.name
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
    end
    return title
end

# Write bounds as an interval
function parse_bounds(result)
    lb = result[:lb]
    ub = result[:ub]
    return ": [\$ $(@sprintf("%.3f", lb)), $(@sprintf("%.3f", ub)) \$]"
end
