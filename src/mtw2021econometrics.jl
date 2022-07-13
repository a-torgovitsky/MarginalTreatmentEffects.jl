# Illustration of mutual consistency (Figure 2)
function run_illustrate_mc(savedir::String, compile::Bool = false)
    results = illustrate_mc()
    savedir = joinpath(savedir, project)
    fn_results = joinpath(savedir, "illustrate-mc.csv")
    CSV.write(fn_results, results)
    if compile
        compile_latex(joinpath(savedir, "FigureMTRConsistent.tex"))
        compile_latex(joinpath(savedir, "FigureMTRInconsistent.tex"))
    end
    return results
end
export run_illustrate_mc

# Numerical illustration for ATT (Figure 4)
function run_simulation_att(savedir::String, compile::Bool = false)
    results = simulation_att()
    savedir = joinpath(savedir, project)
    fn_results = joinpath(savedir, "simulation-att.csv")
    CSV.write(fn_results, results)
    if compile
        compile_latex(joinpath(savedir, "FigureATT.tex"))
    end
    return results
end
export run_simulation_att

# Numerical illustration for LATE(+20) (Figure 5)
function run_simulation_prte(savedir::String, compile::Bool = false)
    results = simulation_prte()
    savedir = joinpath(savedir, project)
    fn_results = joinpath(savedir, "simulation-prte.csv")
    CSV.write(fn_results, results)
    if compile
        compile_latex(joinpath(savedir, "FigurePRTE.tex"))
    end
    return results
end
export run_simulation_prte

# PRTE misspecification example (Figure 6)
function run_prte_misspecification(savedir::String, compile::Bool = false)
    results = prte_misspecification()
    savedir = joinpath(savedir, project)
    fn_results = joinpath(savedir, "prte-misspecification.csv")
    CSV.write(fn_results, results)
    if compile
        compile_latex(joinpath(savedir, "FigurePRTEMisspecification.tex"))
    end
    return results
end
export run_prte_misspecification
