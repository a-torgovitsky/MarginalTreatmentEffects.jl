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
