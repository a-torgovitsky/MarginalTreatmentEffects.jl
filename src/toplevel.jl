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
    println("\t 0. All figures")
    if project_choice == 1
        dirname = "mst2018econometrica" # tex/dirname contains tex templates
        println("\t 1. DGP MTRs and Weights for LATE & IV Slope (Figure 1)")
        println("\t 2. LATE Bounds w/ IV Slope (Figure 2)")
        println("\t 3. LATE Bounds w/ IV & OLS Slopes (Figure 3)")
        println("\t 4. LATE Bounds w/ two IV Slopes (Figure 4)")
        println("\t 5. Sharp LATE Bounds (Figure 5)")
        println("\t 6. Sharp LATE Bounds w/ decr. MTRs (Figure 6)")
        println("\t 7. Sharp LATE Bounds w/ decr., 9th degree MTRs (Figure 7)")
        println("\t 8. Bounds on Family of PRTEs (Figure 8)")
    elseif project_choice == 2
        dirname = "mt2018review"
        println("\t 1. DGP MTRs and MTE (Figure 1)")
        println("\t 2. Weights for Conventional Target Parameters (Figure 2)")
        println("\t 3. Extrapolated LATEs (Figure 3)")
        println("\t 4. ATT Bounds w/ 4th degree MTRs (Figure 4)")
        println("\t 5. ATT Bounds w/ 9th degree MTRs (Figure 5)")
        println("\t 6. ATT Bounds w/ different MTR restrictions (Figure 6)")
        println("\t 7. ATT Bounds w/ nonparametric MTRs (Figure 7)")
        println("\t 8. ATT Bounds w/ 9th degree, decr. MTRs (Figure 8)")
        println("\t 9. ATT Bounds w/ more IV-like estimands (Figure 9)")
        println("\t 10. LATE Bounds w/ different IV-like estimands (Figure 10)")
        println("\t 11. LATE Bounds w/ different MTR restrictions (Figure 11)")
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
        else
            @error "WIP" project_choice figure_choice
        end
    elseif project_choice == 2
        @error "WIP" project_choice figure_choice
    end
end
export menu

# Figure 1: DGP MTRs and weights for LATE(0.35, 0.90) and IV Slope
function run_np_ivs_notitle(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    basis = [(constantspline_basis(knots),
              constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslope => true
    )
    opts = defaults_econometrica()
    opts[1][:title] = "~"
    texfn = mtrs_and_weights(savedir, "np-ivs-no-title";
        dgp = dgp,
        tp = late(dgp, 0.35, 0.9),
        basis = basis,
        assumptions = assumptions,
        mtroption = "truth",
        opts = opts
    )
    if compile
        compile_latex(texfn)
    end
end

# Figure 2: maximizing MTRs for LATE(0.35, 0.90) with IV Slope Estimand
function run_np_ivs(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    basis = [(constantspline_basis(knots),
              constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslope => true
    )
    opts = defaults_econometrica()
    opts[1][:title] = "Nonparametric bounds"
    texfn = mtrs_and_weights(savedir, "np-ivs";
        dgp = dgp,
        tp = late(dgp, 0.35, 0.9),
        basis = basis,
        assumptions = assumptions,
        mtroption = "max",
        opts = opts
    )
    if compile
        compile_latex(texfn)
    end
end

# Figure 3: maximizing MTRs for LATE(0.35, 0.90) with IV and OLS Slope Estimand
function run_np_ivs_olss(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    basis = [(constantspline_basis(knots),
              constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslope => true,
        :olsslope => true,
    )
    opts = defaults_econometrica()
    opts[1][:title] = "Nonparametric bounds"
    texfn = mtrs_and_weights(savedir, "np-ivs-olss";
        dgp = dgp,
        tp = late(dgp, 0.35, 0.9),
        basis = basis,
        assumptions = assumptions,
        mtroption = "max",
        opts = opts
    )
    if compile
        compile_latex(texfn)
    end
end

# Figure 4: maximizing MTRs for LATE(0.35, 0.90) with Nonparametric IV Slope
function run_np_ivnps(savedir::String, compile::Bool = false)
    dgp = dgp_econometrica()
    knots = vcat(0, 1, dgp.pscore, 0.35, 0.9)
    basis = [(constantspline_basis(knots),
              constantspline_basis(knots))]
    assumptions = Dict{Symbol, Any}(
        :lb => 0,
        :ub => 1,
        :saturated => false,
        :ivslopeind => [2, 3],
    )
    opts = defaults_econometrica()
    opts[1][:title] = "Nonparametric bounds"
    texfn = mtrs_and_weights(savedir, "np-ivnps";
        dgp = dgp,
        tp = late(dgp, 0.35, 0.9),
        basis = basis,
        assumptions = assumptions,
        mtroption = "max",
        opts = opts
    )
    if compile
        compile_latex(texfn)
    end
end

# Q: move to src/utils.jl?
function compile_latex(fn::String)
    oldwd = pwd()
    try
        cd(dirname(fn))
        cstr = `pdflatex -halt-on-error $(basename(fn)) "|" grep -a3 ^!`
        @suppress begin
            run(cstr)
            run(cstr) # Q: need to run twice to get references correct
            # Q: why not use latexmk to compile pdf?
            # i.e. run(`latexmk -pdf $(basename(fn))`)
            run(`latexmk -c`)
        end
        cd(oldwd)
    catch err
        cd(oldwd)
    end
end
