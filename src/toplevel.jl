function menu(savelocation::String = "."; compile::Bool = false)
    println("="^80)
    println("Marginal Treatment Effects")
    println("="^80)
    println("What project do you want to replicate?")
    println("\t 1. Using Instrumental Variables for Inference About \r\t
            Policy Relevant Treatment Parameters")
    println("\t\t Mogstad, Santos, and Torgovitsky (2018, Econometrica)")
    println("\t 2. Identification and Extrapolation of Causal Effects \r\t
            with Instrumental Variables")
    println("\t\t Mogstad and Torgovitsky (2018, Annual Review of Economics)")
    print("Enter project choice: ")
    project_choice = readline()
    project_choice = parse(Int64, project_choice)

    println("What figure do you want to reproduce?")
    println("\t 0. Everything")
    if project_choice == 1
        global project = "mst2018econometrica" # tex/dirname contains tex templates
        println("\t 1. DGP MTRs and Weights for LATE & IV Slope (Figure 1)")
        println("\t 2. LATE Bounds w/ IV Slope (Figure 2)")
        println("\t 3. LATE Bounds w/ IV & OLS Slopes (Figure 3)")
        println("\t 4. LATE Bounds w/ Nonparametric IV Slopes (Figure 4)")
        println("\t 5. Sharp LATE Bounds (Figure 5)")
        println("\t 6. Sharp LATE Bounds w/ Decr. MTRs (Figure 6)")
        println("\t 7. Sharp LATE Bounds w/ Decr., 9th degree MTRs (Figure 7)")
        println("\t 8. Bounds on Family of PRTEs (Figure 8)")
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
    else
        @error "Invalid answer: only choose 1 or 2"
    end
    print("Enter figure choice: ")
    figure_choice = readline()
    figure_choice = parse(Int64, figure_choice)

    if project_choice == 1
        if figure_choice == 1
            savedir, _ = setup(savelocation, stub = "np-ivs-no-title")
            run_np_ivs_notitle(savedir, compile)
        elseif figure_choice == 2
            savedir, _ = setup(savelocation, stub = "np-ivs")
            run_np_ivs(savedir, compile)
        elseif figure_choice == 3
            savedir, _ = setup(savelocation, stub = "np-ivs-olss")
            run_np_ivs_olss(savedir, compile)
        elseif figure_choice == 4
            savedir, _ = setup(savelocation, stub = "np-ivnps")
            run_np_ivnps(savedir, compile)
        elseif figure_choice == 5
            savedir, _ = setup(savelocation, stub = "np-sharp")
            run_np_sharp(savedir, compile)
        elseif figure_choice == 6
            savedir, _ = setup(savelocation, stub = "np-sharp-decr")
            run_np_sharp_decr(savedir, compile)
        elseif figure_choice == 7
            # NOTE: in appendix, K refers to degree, not order
            savedir, _ = setup(savelocation, stub = "np-sharp-decr-k9")
            run_np_sharp_decr_k9(savedir, compile)
        elseif figure_choice == 8
            savedir, _ = setup(savelocation, stub = "tikz-extrapolate")
            run_tikz_extrapolate(savedir, compile)
        elseif figure_choice == 0
            savedir, _ = setup(savelocation, stub = "everything")
            run_np_ivs_notitle(savedir, compile)
            run_np_ivs(savedir, compile)
            run_np_ivs_olss(savedir, compile)
            run_np_ivnps(savedir, compile)
            run_np_sharp(savedir, compile)
            run_np_sharp_decr(savedir, compile)
            run_np_sharp_decr_k9(savedir, compile)
            run_tikz_extrapolate(savedir, compile)
        else
            @error "ERROR: invalid choice" project_choice figure_choice
        end
    elseif project_choice == 2
        if figure_choice == 1
            savedir, _ = setup(savelocation, stub = "tikz-mtr")
            run_tikz_mtr(savedir, compile)
        elseif figure_choice == 2
            savedir, _ = setup(savelocation, stub = "tikz-weights")
            run_tikz_weights(savedir, compile)
        elseif figure_choice == 3
            savedir, _ = setup(savelocation, stub = "tikz-late-extrap")
            run_tikz_late_extrap(savedir, compile)
        elseif figure_choice == 4
            savedir, _ = setup(savelocation, stub = "k4")
            run_k4(savedir, compile)
        elseif figure_choice == 5
            savedir, _ = setup(savelocation, stub = "k9")
            run_k9(savedir, compile)
        elseif figure_choice == 6
            savedir, _ = setup(savelocation, stub = "kbounds")
            run_kbounds(savedir, compile)
        elseif figure_choice == 7
            savedir, _ = setup(savelocation, stub = "np")
            run_np(savedir, compile)
        elseif figure_choice == 8
            savedir, _ = setup(savelocation, stub = "k9-decr")
            run_k9_decr(savedir, compile)
        elseif figure_choice == 9
            savedir, _ = setup(savelocation, stub = "k9-decr-add-more")
            run_k9_decr_add_more(savedir, compile)
        elseif figure_choice == 10
            savedir, _ = setup(savelocation, stub = "late-bounds-information")
            run_late_bounds_information(savedir, compile)
        elseif figure_choice == 11
            savedir, _ = setup(savelocation, stub = "late-bounds-assumptions")
            run_late_bounds_assumptions(savedir, compile)
        elseif figure_choice == 0
            savedir, _ = setup(savelocation, stub = "everything")
            run_tikz_mtr(savedir, compile)
            run_tikz_weights(savedir, compile)
            run_tikz_late_extrap(savedir, compile)
            run_k4(savedir, compile)
            run_k9(savedir, compile)
            run_kbounds(savedir, compile)
            run_np(savedir, compile)
            run_k9_decr(savedir, compile)
            run_k9_decr_add_more(savedir, compile)
            run_late_bounds_information(savedir, compile)
            run_late_bounds_assumptions(savedir, compile)
        else
            @error "WIP" project_choice figure_choice
        end
    end
end
export menu
