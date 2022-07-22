# Generate tag from user input
function generate_tag()
    print("Enter a directory tag "*
          "(or hit Enter to generate from the date): ")
    tag = readline()
    if isempty(tag)
        tag = Libc.strftime("%d-%b-%y-at-%H-%M-%S", time())
    end
    return tag
end

# Create the name of the directory where reproductions will be stored.
function make_dirname(savelocation::String;
                      stub::String = "unnamed",
                      tag::Union{String, Nothing} = nothing)
    println("Saving in: $savelocation")
    if isnothing(tag)
        tag = generate_tag()
    end

    if isempty(stub)
        dirname = tag
    else
        dirname = stub*"-"*tag
    end
    pathstr = joinpath(savelocation, dirname)
    return pathstr
end

# Create self-contained directory where the reproductions will be stored.
# If `copysource` is true, then source code will also be copied into `pathstr`
# and the results will be in a dedicated `results` folder.
function make_dir(pathstr::String; copysource::Bool = false)
    # Create directory that contains output, including source code and results.
    if isdir(pathstr)
        @info "Output directory already exists."
    else
        mkpath(pathstr)
    end

    # Create relevant paths.
    if copysource
        resultsdir = joinpath(pathstr, "results")
    else
        resultsdir = joinpath(pathstr)
    end
    projectdir = joinpath(resultsdir, project) # project-specific results
    sourcepathstr = joinpath(@__DIR__, "../")  # source code

    # Create project-specific directory for holding results.
    already_existed = false
    if isdir(projectdir)
        @info "Results directory for $project already exists, not modifying."
        already_existed = true
    else
        mkpath(projectdir) # creates `resultsdir` if it does not already exist
    end

    # Copy source code and tex templates if `copysource` is true.
    if copysource
        cplist = ["LICENSE", "Project.toml", "README.md",
                  "src", "tex"]
        for c in cplist
            source = joinpath(sourcepathstr, c)
            destination = joinpath(pathstr, c)
            if ispath(destination)
                @info "Source code already exists, so not copying."
                break
            else
                cp(source, destination, follow_symlinks=true)
            end
        end
    end

    # Copy project-specific .tex files to results directory.
    texdir = joinpath(sourcepathstr, "tex", project)
    cplist = readdir(texdir)
    for c in cplist
        source = joinpath(texdir, c)
        destination = joinpath(projectdir, c)
        if ispath(destination)
            @info destination * " already exists, so not copying tex files."
            break
        else
            cp(source, destination)
        end
    end

    return resultsdir, already_existed
end

# Set up directory for reproduction
function setup(savelocation::String;
               stub::String="unnamed",
               tag::Union{String, Nothing}=nothing,
               copysource::Bool = false)
    pathstr = make_dirname(savelocation, stub = stub, tag = tag)
    resultsdir, already_existed = make_dir(pathstr, copysource = copysource)
    return resultsdir, already_existed
end

# Compile tex files
function compile_latex(fn::String)
    oldwd = pwd()
    try
        cd(dirname(fn))
        cstr = `pdflatex -halt-on-error $(basename(fn)) "|" grep -a3 ^!`
        byproducts = ["aux", "log"]
        @suppress begin
            run(cstr)
            run(cstr) # run twice to get references correct
            for extension in byproducts
                remove = "$(splitext(basename(fn))[1]).$(extension)"
                run(`rm $remove`)
            end
        end
        cd(oldwd)
    catch err
        cd(oldwd)
    end
end
