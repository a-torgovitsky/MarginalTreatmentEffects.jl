################################################################################
# Goal: compute the average (across instruments) weights on MTRs
#
# The average weights are piecewise-constant functions of u.
# The `compute_average_weights` function will return a data frame whose columns
# are:
#   1. value of u
#   2. average weight for d = 1
#   3. average weight for d = 0
#
# Note that the weight for [a, b] or [a, b) will be associated with u = a, not
# u = b.
# This convention makes it easier to identify the segments.
#
# Here's an example:
#
# |  u   | average weight for d = 1 | average weight for d = 0 |
#  ------ -------------------------- --------------------------
# | 0.0  |         6.66134e-16      |         0.0              |
# | 0.12 |         1.2392           |        -1.2392           |
# | 0.29 |         1.77028          |        -1.77028          |
# | 0.48 |         1.50994          |        -1.50994          |
# | 0.78 |         0.0              |         6.66134e-16      |
#
# This result suggests that:
#   - for u ∈ [0.0, 0.12), the average weight is 0 for d ∈ {0, 1};
#   - for u ∈ [0.12, 0.29), the average weight for d = 1 is 1.2392;
#   - for u ∈ [0.29, 0.48), the average weight for d = 1 is 1.77028;
#   - for u ∈ [0.48, 0.78), the average weight for d = 1 is 1.50994;
#   - for u ∈ [0.78, 1.00], the average weight for d = 1 is 0;
#   - for u ∈ [0.0,  1.00], the average weight for d = 0 is
#         -1 * the average weight for d = 1;
#
# The average weights for IV-like estimands rely on the s(d, z) function.
# The average weights for target parameters are derived by hand.
#
################################################################################

function compute_average_weights(tp::TargetParameter,
                                 dgp::Union{DGP, Nothing} = nothing)
    if tp.name in ["ATT", "ATU"]
        @assert !isnothing(dgp) "Need to specify dgp for $(tp.name)"
    end
    if tp.name == "LATE(u₁, u₂)"
        # Only coded the case of 1 instrument
        u₁ = tp.int_limits(1)[1]
        u₂ = tp.int_limits(1)[2]
        results = DataFrame(u = [0, u₁, u₂, 1])
        results[:, "average weight for d = 1"] .= 0.0
        results[:, "average weight for d = 0"] .= 0.0
        # Only non-trivial weight is 1 / (u₂ - u₁) on the interval [u₁, u₂].
        results[2, 2] =  1 / (u₂ - u₁)
        results[2, 3] = -1 / (u₂ - u₁)
    elseif tp.name == "ATE"
        results = DataFrame(u = [0, 1])
        results[:, "average weight for d = 1"] .= 0.0
        results[:, "average weight for d = 0"] .= 0.0
        # Only non-trivial weight is on the interval [0, 1].
        results[1, 2] =  1.0
        results[1, 3] = -1.0
    elseif tp.name == "ATT"
        results = DataFrame(u = [0, dgp.pscore..., 1])
        results[:, "average weight for d = 1"] .= 0.0
        results[:, "average weight for d = 0"] .= 0.0
        prd1 = dot(dgp.pscore, dgp.densZ)
        weights = reverse(cumsum(dgp.densZ) ./ prd1)
        results[1:length(dgp.pscore), 2] = weights
        results[1:length(dgp.pscore), 3] = -weights
    elseif tp.name == "ATU"
        results = DataFrame(u = [0, dgp.pscore..., 1])
        results[:, "average weight for d = 1"] .= 0.0
        results[:, "average weight for d = 0"] .= 0.0
        prd1 = dot(dgp.pscore, dgp.densZ)
        weights = cumsum(reverse(dgp.densZ)) ./ (1 - prd1)
        results[2:(1 + length(dgp.pscore)), 2] = weights
        results[2:(1 + length(dgp.pscore)), 3] = -weights
    else
        print("Weight Computation is unsupported for this target parameter.")
        return
    end
    return results
end

# The weights can be written in terms of s-like functions:
#     ω₁(u, z) ≡ s(1, z)𝟙[u ≤ p(z)]
#     ω₀(u, z) ≡ s(0, z)𝟙[u > p(z)]
# The average weight across the instruments is therefore given by
#     𝔼[ω₁(u, z)] ≡ ∑ s(1, z)𝟙[u ≤ p(z)] ⋅ ℙ(Z = z)
#     𝔼[ω₀(u, z)] ≡ ∑ s(0, z)𝟙[u > p(z)] ⋅ ℙ(Z = z)
# where the summation is across the support points of Z.
# As u moves from 0 to 1, the indicator controls how many terms enter the sum.
# For 𝔼[ω₁(u, z)], the number of terms decreases; the terms associated with
# smaller values of p(z) leave the summation first.
# For 𝔼[ω₀(u, z)], the number of terms increases; the terms associated with
# smaller values of p(z) enter the summation first.
function compute_average_weights(ivlike::IVLike, dgp::DGP)
    @assert length(ivlike.s) == 1 "Only 1 s-like function is supported."
    results = DataFrame(u = unique(vcat(0, dgp.pscore)))
    order = sortperm(dgp.pscore)
    s = ivlike.s[1] # only 1 s-like function at a time!
    # NOTE: `eachrow()` in the `make_slist()` function yields vectors
    if occursin("Saturated", ivlike.name)
        terms = d -> s.(d, eachrow(dgp.suppZ)) .* dgp.densZ
    else
        terms = d -> s.(d, dgp.suppZ) .* dgp.densZ
    end
    # d = 1 weights
    d1terms = terms(1)
    summands = d1terms[reverse(order)] # order by decreasing pscore
    pushfirst!(summands, 0)
    results[:, "average weight for d = 1"] = reverse(cumsum(summands))
    # d = 0 weights
    d0terms = terms(0)
    summands = d0terms[order] # order by increasing pscore
    pushfirst!(summands, 0)
    results[:, "average weight for d = 0"] = cumsum(summands)
    return results
end
