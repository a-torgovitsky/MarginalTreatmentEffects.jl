################################################################################
# Goal: compute the average (across instruments) weights on MTRs
#
# The average weights are piecewise-constant functions of u.
# The `compute_average_weights` function will return a data frame that stores
# value of the average weight at the endpoints of each nontrivial piece.
#   1. value of u
#   2. average weight for d = 1
#   3. average weight for d = 0
#
# The average weights for IV-like estimands rely on the s(d, z) function.
# The average weights for target parameters are derived manually.
#
################################################################################

function compute_average_weights(tp::TargetParameter)
    if tp.name == "LATE(u₁, u₂)"
        u₁ = tp.int_limits(1)[1] # only coded the case of 1 instrument
        u₂ = tp.int_limits(1)[2] # only coded the case of 1 instrument
        results = DataFrame(u = [0, u₁, u₂, 1])
        results[:, "average weight for d = 1"] .= 0.0
        results[:, "average weight for d = 0"] .= 0.0

        # only non-trivial weight is 1 / (u₂ - u₁) on the interval [u₁, u₂]
        results[2, "average weight for d = 1"] = 1 / (u₂ - u₁)
        results[2, "average weight for d = 0"] = -1 / (u₂ - u₁)
    else
        print("Weight Computation is unsupported for this target parameter.")
    end
    return results
end
export compute_average_weights

function compute_average_weights(ivlike::IVLike, dgp::DGP)
    # only 1 s-list at a time!
    # TODO: allow multiple entries in ivlike.s?
    results = DataFrame(u = unique(vcat(0, dgp.pscore)))
    order = sortperm(dgp.pscore)
    s = ivlike.s[1] # only 1 s-list at a time!
    if occursin("Saturated", ivlike.name)
        # `eachrow()` in the `make_slist()` function yields vectors
        # TODO: resolve this inconsistency
        terms = d -> s.(d, eachrow(dgp.suppZ)) .* dgp.densZ
    else
        terms = d -> s.(d, dgp.suppZ) .* dgp.densZ
    end
    # d = 1
    d1terms = terms(1)
    summands = d1terms[reverse(order)] # order by decreasing pscore
    pushfirst!(summands, 0)
    results[:, "average weight for d = 1"] = reverse(cumsum(summands))
    # d = 0
    d0terms = terms(0)
    summands = d0terms[order] # order by increasing pscore
    pushfirst!(summands, 0)
    results[:, "average weight for d = 0"] = cumsum(summands)
    return results
end
