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
