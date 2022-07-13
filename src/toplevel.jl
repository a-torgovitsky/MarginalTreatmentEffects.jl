function menu(savelocation::String = "."; compile::Bool = false)

    # Initialization
    global texfiles = Vector{String}() # store path of tex files

    # Present project options.
    println("="^80)
    println("Marginal Treatment Effects")
    println("="^80)
    println("What project do you want to replicate?")
    println("\t 0. Everything")
    println("\t 1. Using Instrumental Variables for Inference About \r\t
            Policy Relevant Treatment Parameters")
    println("\t\t Mogstad, Santos, and Torgovitsky (2018, Econometrica)")
    println("\t 2. Identification and Extrapolation of Causal Effects \r\t
            with Instrumental Variables")
    println("\t\t Mogstad and Torgovitsky (2018, Annual Review of Economics)")
    println("\t 3. Policy Evaluation with Multiple Instrumental Variables")
    println("\t\t Mogstad, Torgovitsky, and Walters " *
            "(2021, Journal of Econometrics)")
    print("Enter project choice: ")
    project_choice = readline()
    project_choice = parse(Int64, project_choice)
    @assert project_choice in [0, 1, 2, 3] "Choose 0, 1, 2, or 3."

    # Run everything and short-circuit, or present figure options
    if project_choice == 0
        tag = generate_tag()
        global project = "mst2018econometrica" # tex/project contains templates
        run_mst2018econometrica(0, savelocation; tag = tag)
        project = "mt2018review"
        run_mt2018review(0, savelocation; tag = tag)
        project = "mtw2021econometrics"
        # TODO: save tex files and compile later instead of compiling here
        run_mtw2018econometrics(0, savelocation, compile; tag = tag)
    else
        println("What figure do you want to reproduce?")
        println("\t 0. Everything")
    end

    # Present figure options.
    if project_choice == 1
        global project = "mst2018econometrica"
        println("\t 1. DGP MTRs and Weights for LATE & IV Slope (Figure 1)")
        println("\t 2. LATE Bounds w/ IV Slope (Figure 2)")
        println("\t 3. LATE Bounds w/ IV & OLS Slopes (Figure 3)")
        println("\t 4. LATE Bounds w/ Nonparametric IV Slopes (Figure 4)")
        println("\t 5. Sharp LATE Bounds (Figure 5)")
        println("\t 6. Sharp LATE Bounds w/ Decr. MTRs (Figure 6)")
        println("\t 7. Sharp LATE Bounds w/ Decr., 9th degree MTRs (Figure 7)")
        println("\t 8. Bounds on Family of PRTEs (Figure 8)")
        valid_choices = collect(1:8)
    elseif project_choice == 2
        global project = "mt2018review"
        println("\t 1. DGP MTRs and MTE (Figure 1)")
        println("\t 2. Weights for Conventional Target Parameters (Figure 2)")
        println("\t 3. LATE Extrapolation (Figure 3)")
        println("\t 4. ATT Bounds w/ 4th degree MTRs (Figure 4)")
        println("\t 5. ATT Bounds w/ 9th degree MTRs (Figure 5)")
        println("\t 6. ATT Bounds w/ Different MTR Assumptions (Figure 6)")
        println("\t 7. ATT Bounds w/ Nonparametric MTRs (Figure 7)")
        println("\t 8. ATT Bounds w/ Decr., 9th degree MTRs (Figure 8)")
        println("\t 9. ATT Bounds w/ more IV-like estimands (Figure 9)")
        println("\t 10. LATE Bounds w/ Different Information Sets (Figure 10)")
        println("\t 11. LATE Bounds w/ Different MTR Assumptions (Figure 11)")
        valid_choices = collect(1:11)
    elseif project_choice == 3
        global project = "mtw2021econometrics"
        println("\t 1. Illustration of mutual consistency (Figure 2)")
        println("\t 2. Numerical illustration for ATT (Figure 4)")
        println("\t 3. Numerical illustration for LATE(+20) (Figure 5)")
        println("\t 4. PRTE misspecification example (Figure 6)")
        valid_choices = collect(1:4)
    end

    if project_choice != 0
        print("Enter figure choice: ")
        figure_choice = readline()
        figure_choice = parse(Int64, figure_choice)
        valid_choices = vcat(0, valid_choices)
        @assert figure_choice in valid_choices "Choose one of $valid_choices."
    end

    if project_choice == 1
        run_mst2018econometrica(figure_choice, savelocation)
    elseif project_choice == 2
        run_mt2018review(figure_choice, savelocation)
    elseif project_choice == 3
        run_mtw2018econometrics(figure_choice, savelocation, compile)
    end

    # produce PDFs from tex files
    if compile
        for texfn in texfiles
            compile_latex(texfn)
        end
    end
end
export menu

# Produce figures in MST (2018)
function run_mst2018econometrica(figure_choice::Int64,
                                 savelocation::String;
                                 tag::Union{String, Nothing} = nothing)
    println("Replicating MST (2018)...")
    if figure_choice == 1
        savedir, _ = setup(savelocation, stub = "np-ivs-no-title", tag = tag)
        push!(texfiles, run_np_ivs_notitle(savedir))
    elseif figure_choice == 2
        savedir, _ = setup(savelocation, stub = "np-ivs", tag = tag)
        push!(texfiles, run_np_ivs(savedir))
    elseif figure_choice == 3
        savedir, _ = setup(savelocation, stub = "np-ivs-olss", tag = tag)
        push!(texfiles, run_np_ivs_olss(savedir))
    elseif figure_choice == 4
        # We produce 2 figures here.
        # See the note in the source code for `run_np_ivnps()`, which can
        # be found in `src/mst2018econometrica.jl`.
        savedir, _ = setup(savelocation, stub = "np-ivnps", tag = tag)
        push!(texfiles, run_np_ivnps(savedir)...)
    elseif figure_choice == 5
        savedir, _ = setup(savelocation, stub = "np-sharp", tag = tag)
        push!(texfiles, run_np_sharp(savedir))
    elseif figure_choice == 6
        savedir, _ = setup(savelocation, stub = "np-sharp-decr", tag = tag)
        push!(texfiles, run_np_sharp_decr(savedir))
    elseif figure_choice == 7
        # NOTE: in appendix, K refers to degree, not order
        savedir, _ = setup(savelocation, stub = "np-sharp-decr-k9", tag = tag)
        push!(texfiles, run_np_sharp_decr_k9(savedir))
    elseif figure_choice == 8
        savedir, _ = setup(savelocation, stub = "tikz-extrapolate", tag = tag)
        push!(texfiles, run_tikz_extrapolate(savedir))
    elseif figure_choice == 0
        savedir, _ = setup(savelocation, stub = "everything", tag = tag)
        push!(texfiles, run_np_ivs_notitle(savedir))
        push!(texfiles, run_np_ivs(savedir))
        push!(texfiles, run_np_ivs_olss(savedir))
        push!(texfiles, run_np_ivnps(savedir)...)
        push!(texfiles, run_np_sharp(savedir))
        push!(texfiles, run_np_sharp_decr(savedir))
        push!(texfiles, run_np_sharp_decr_k9(savedir))
        push!(texfiles, run_tikz_extrapolate(savedir))
    else
        @error "ERROR: invalid choice" project_choice figure_choice
    end
end

# Produce figures in MT (2018)
function run_mt2018review(figure_choice::Int64,
                          savelocation::String;
                          tag::Union{String, Nothing} = nothing)
    println("Replicating MST (2018)...")
    if figure_choice == 1
        savedir, _ = setup(savelocation, stub = "tikz-mtr", tag = tag)
        push!(texfiles, run_tikz_mtr(savedir))
    elseif figure_choice == 2
        savedir, _ = setup(savelocation, stub = "tikz-weights", tag = tag)
        push!(texfiles, run_tikz_weights(savedir))
    elseif figure_choice == 3
        savedir, _ = setup(savelocation, stub = "tikz-late-extrap", tag = tag)
        push!(texfiles, run_tikz_late_extrap(savedir))
    elseif figure_choice == 4
        savedir, _ = setup(savelocation, stub = "k4", tag = tag)
        push!(texfiles, run_k4(savedir))
    elseif figure_choice == 5
        savedir, _ = setup(savelocation, stub = "k9", tag = tag)
        push!(texfiles, run_k9(savedir))
    elseif figure_choice == 6
        savedir, _ = setup(savelocation, stub = "kbounds", tag = tag)
        push!(texfiles, run_kbounds(savedir))
    elseif figure_choice == 7
        savedir, _ = setup(savelocation, stub = "np", tag = tag)
        push!(texfiles, run_np(savedir))
    elseif figure_choice == 8
        savedir, _ = setup(savelocation, stub = "k9-decr", tag = tag)
        push!(texfiles, run_k9_decr(savedir))
    elseif figure_choice == 9
        savedir, _ = setup(savelocation, stub = "k9-decr-add-more", tag = tag)
        push!(texfiles, run_k9_decr_add_more(savedir))
    elseif figure_choice == 10
        savedir, _ = setup(savelocation,
                           stub = "late-bounds-information",
                           tag = tag)
        push!(texfiles, run_late_bounds_information(savedir))
    elseif figure_choice == 11
        savedir, _ = setup(savelocation,
                           stub = "late-bounds-assumptions",
                           tag = tag)
        push!(texfiles, run_late_bounds_assumptions(savedir))
    elseif figure_choice == 0
        savedir, _ = setup(savelocation, stub = "everything", tag = tag)
        push!(texfiles, run_tikz_mtr(savedir))
        push!(texfiles, run_tikz_weights(savedir))
        push!(texfiles, run_tikz_late_extrap(savedir))
        push!(texfiles, run_k4(savedir))
        push!(texfiles, run_k9(savedir))
        push!(texfiles, run_kbounds(savedir))
        push!(texfiles, run_np(savedir))
        push!(texfiles, run_k9_decr(savedir))
        push!(texfiles, run_k9_decr_add_more(savedir))
        push!(texfiles, run_late_bounds_information(savedir))
        push!(texfiles, run_late_bounds_assumptions(savedir))
    else
        @error "ERROR: invalid choice" project_choice figure_choice
    end
end

function run_mtw2018econometrics(figure_choice::Int64,
                                 savelocation::String,
                                 compile::Bool;
                                 tag::Union{String, Nothing} = nothing)
    println("Replicating MST (2018)...")
    if figure_choice == 1
        savedir, _ = setup(savelocation, stub = "illustrate-mc", tag = tag)
        run_illustrate_mc(savedir, compile)
    elseif figure_choice == 2
        savedir, _ = setup(savelocation, stub = "simulation-att", tag = tag)
        run_simulation_att(savedir, compile)
    elseif figure_choice == 3
        savedir, _ = setup(savelocation, stub = "simulation-prte", tag = tag)
        run_simulation_prte(savedir, compile)
    elseif figure_choice == 4
        savedir, _ = setup(savelocation,
                           stub = "prte-misspecification",
                           tag = tag)
        run_prte_misspecification(savedir, compile)
        elseif figure_choice == 0
        savedir, _ = setup(savelocation, stub = "everything", tag = tag)
        results_mc = run_illustrate_mc(savedir, compile)
        results_att = run_simulation_att(savedir, compile)
        results_prte = run_simulation_prte(savedir, compile)
        results_prte_misspecification =
            run_prte_misspecification(savedir, compile)
        return Dict(:results_mc => results_mc,
                    :results_att => results_att,
                    :results_prte => results_prte,
                    :results_prte_misspecification =>
                        results_prte_misspecification)
    else
        @error "ERROR: invalid choice" project_choice figure_choice
    end
end
