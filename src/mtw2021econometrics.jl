# Illustration of mutual consistency (Figure 2)
function run_illustrate_mc(savedir::String)
    results = illustrate_mc()
    savedir = joinpath(savedir, project)
    fn_results = joinpath(savedir, "illustrate-mc.csv")
    CSV.write(fn_results, results)
    texfn1 = joinpath(savedir, "FigureMTRConsistent.tex")
    texfn2 = joinpath(savedir, "FigureMTRInConsistent.tex")
    return results, texfn1, texfn2
end
export run_illustrate_mc

# Numerical illustration for ATT (Figure 4)
function run_simulation_att(savedir::String)
    results = simulation_att()
    savedir = joinpath(savedir, project)
    fn_results = joinpath(savedir, "simulation-att.csv")
    CSV.write(fn_results, results)
    texfn = joinpath(savedir, "FigureATT.tex")
    return results, texfn
end
export run_simulation_att

# Numerical illustration for LATE(+20) (Figure 5)
function run_simulation_prte(savedir::String)
    results = simulation_prte()
    savedir = joinpath(savedir, project)
    fn_results = joinpath(savedir, "simulation-prte.csv")
    CSV.write(fn_results, results)
    texfn = joinpath(savedir, "FigurePRTE.tex")
    return results, texfn
end
export run_simulation_prte

# PRTE misspecification example (Figure 6)
function run_prte_misspecification(savedir::String)
    results = prte_misspecification()
    savedir = joinpath(savedir, project)
    fn_results = joinpath(savedir, "prte-misspecification.csv")
    CSV.write(fn_results, results)
    texfn = joinpath(savedir, "FigurePRTEMisspecification.tex")
    return results, texfn
end
export run_prte_misspecification
