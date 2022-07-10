# Create the name of the directory where reproductions will be stored.
function make_dirname(savelocation::String;
                      stub::String = "unnamed",
                      tag::String = "")
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

# Create self-contained directory where the reproductions will be stored.
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
        destination = joinpath(resultsdir, c)
        if isdir(destination)
            @info destination * " already exists, so not copying tex files."
        else
            cp(joinpath(texdir, c), destination)
        end
    end

    return resultsdir, already_existed
end

# Set up directory for reproduction
function setup(savelocation::String;
               stub::String="unnamed",
               tag::String="")
    pathstr = make_dirname(savelocation, stub = stub, tag = tag)
    resultsdir, already_existed = make_dir(pathstr)
    return resultsdir, already_existed
end

# Compile tex files
function compile_latex(fn::String)
    oldwd = pwd()
    try
        cd(dirname(fn))
        cstr = `pdflatex -halt-on-error $(basename(fn)) "|" grep -a3 ^!`
        @suppress begin
            run(cstr)
            run(cstr) # NOTE: need to run twice to get references correct
            # Q(a-torgovitsky): why not use latexmk to compile pdf?
            # i.e. run(`latexmk -pdf $(basename(fn))`)
            run(`latexmk -c`)
        end
        cd(oldwd)
    catch err
        cd(oldwd)
    end
end
