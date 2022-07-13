################################################################################
# For tractability, a target parameter must have the form:
#   \sum_{z, \ell, d}
#       int_constant(ℓ, d, z)
#       \times
#       \int_{int_limits(z)[1]}^{int_limits(z)[2]}
#       m(u_{ℓ} | Z_{-ℓ}) du
#
#   So defining a target parameter requires a pair of functions:
#   - int_limits, which takes z and returns a lower and upper limit
#   of integration
#   - int_constant, which takes (ℓ, d, z) and returns a scalar.
#
#   Note that expectation over Z is incorporated into int_constant(ℓ, d, z)
################################################################################
struct TargetParameter
    name::String
    int_limits::Function
    int_constant::Function
    legendtitle::String

    function TargetParameter(name,
                             int_limits,
                             int_constant,
                             legendtitle = nothing)
        if isnothing(legendtitle)
            legendtitle = name
        end
        new(name, int_limits, int_constant, legendtitle)
    end
end
export TargetParameter

function eval_tp(tp::TargetParameter, mtrs::Vector{Tuple{MTR, MTR}}, dgp)
    Γ⭑ = compute_Γ⭑(tp, [(mtr[1].basis, mtr[2].basis) for mtr in mtrs], dgp)
    total = 0
    for ℓ in 1:length(mtrs), d in 0:1
        total += sum(Γ⭑[ℓ, d + 1] .* mtrs[ℓ][d + 1].θ)
    end
    return total
end
export eval_tp

function compute_Γ⭑(tp::TargetParameter,
                    bases::Vector{Tuple{MTRBasis, MTRBasis}}, dgp::DGP)
    [compute_Γ⭑(tp, bases[ℓ][d + 1], d, ℓ, dgp)
     for ℓ in 1:length(bases), d in 0:1]
end

function compute_Γ⭑(tp::TargetParameter, basis::MTRBasis,
                    d::Integer, ℓ::Integer, dgp::DGP)
    Γ⭑ = zeros(length(basis.a), length(basis.b))
    for z in eachrow(dgp.suppZ)
        il = tp.int_limits(z)
        Γ⭑ += [(aj(z) * ibk(il[1], il[2]) * tp.int_constant(ℓ, d, z))
               for aj in basis.a, ibk in basis.ib]
    end
    return Γ⭑
end

function ey1(dgp::DGP; ℓ = 1)
    TargetParameter(
        "EY1",
        z -> (0,1),
        (l,d,z) -> d * (l == ℓ) * find_density(z, dgp),
        "\$\\mathbb{E}[Y(1)]\$"
    )
end
export ey1

function att(dgp::DGP; ℓ = 1)
    prd1 = dot(dgp.pscore, dgp.densZ)
    TargetParameter(
        "ATT",
        z -> (0, find_pscore(z, dgp)),
        (l,d,z) -> ((l == ℓ) * (2*d - 1)/prd1) * find_density(z, dgp)
    )
end
export att

function atu(dgp::DGP; ℓ = 1)
    prd1 = dot(dgp.pscore, dgp.densZ)
    TargetParameter(
        "ATU",
        z -> (find_pscore(z, dgp), 1),
        (l,d,z) -> ((l == ℓ) * (2*d - 1)/(1 - prd1)) * find_density(z, dgp)
    )
end
export atu

function prte_plusδpercent(dgp::DGP, δ; ℓ = 1)
    @assert ℓ == 1 # haven't coded the other case
    limits = z -> (          find_pscore([0, z[2]], dgp),
                   (1 + δ) * find_pscore([1, z[2]], dgp))
    TargetParameter(
        "PRTE_plus$(δ)_ℓ=$(ℓ)",
        limits,
        (l,d,z) -> ((l == ℓ) * (2*d - 1) *
                    find_density(z, dgp) / (limits(z)[2] - limits(z)[1]))
    )
end
export prte_plusδpercent

function ate(dgp::DGP; ℓ = 1)
    TargetParameter(
        "ATE",
        z -> (0, 1),
        (l,d,z) -> ((1 == ℓ) * (2*d - 1))
    )
end
export ate

function late(dgp::DGP, u₁, u₂; ℓ = 1)
    @assert ℓ == 1 # haven't coded the other case
    @assert u₁ <= u₂
    TargetParameter(
        "LATE(u₁, u₂)",
        z -> (u₁, u₂),
        (ℓ, d, z) -> ((1 == ℓ) * (2*d - 1) * find_density(z, dgp) / (u₂ - u₁)),
        "LATE(\$$(@sprintf("%.2f", u₁)), $(@sprintf("%.2f", u₂))\$)",
    )
end
export late

################################################################################
# Special type of PRTE
#   dgpnew and dgpold should have the same "support" points for Z
#   note that you can set their probability to 0 in the old distribution
#   and that's totally fine in the program
#
# All that's shifting here is the marginal distribution of Z
################################################################################
function prte_newz(dgpold::DGP, dgpnew::DGP; ℓ = 1)
    @assert dgpold.suppZ == dgpnew.suppZ
    @assert dgpold.pscore == dgpnew.pscore
    TargetParameter(
        "PRTE New Z",
        z -> (0, find_pscore(z, dgpold)),
        (l, d, z) -> ((l == ℓ) * (2*d - 1) *
                     (find_density(z, dgpnew) - find_density(z, dgpold)))
    )
end
export prte_newz
