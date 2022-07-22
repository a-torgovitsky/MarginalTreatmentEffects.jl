# For imposing
# E[Y(d) | D = d′, Z = z]
#
# Slightly different than the "s-function" formulation in the paper
# Usually equivalent if fully saturated
#
# Useful here for the PRTE misspecification exercise to impose MC on
# E[Y(d) | D = d', Z = z]
# even when Pr[Z = z] = 0
function compute_Γₘ(bases::Array{Tuple{MTRBasis, MTRBasis}, 1}, dgp::DGP)
    [compute_Γₘ(basis[d + 1], dgp) for basis in bases, d in 0:1]
end

function compute_Γₘ(basis::MTRBasis, dgp::DGP)
    Γₘ = zeros(2, size(dgp.suppZ, 1), length(basis.a), length(basis.b))

    for (i,z) in enumerate(eachrow(dgp.suppZ))
        p = dgp.pscore[i]
        # d′ = 0
        Γₘ[1, i, :, :] = [aj(z) * ibk(p, 1) / (1 - p)
                          for aj in basis.a, ibk in basis.ib]
        # d′ = 1
        Γₘ[2, i, :, :] = [aj(z) * ibk(0, p) / p
                          for aj in basis.a, ibk in basis.ib]
    end
    return Γₘ
end
