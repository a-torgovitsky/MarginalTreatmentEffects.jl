function compute_bounds(tp::TargetParameter,
                        bases::Array{Tuple{MTRBasis, MTRBasis}, 1},
                        assumptions::Dict,
                        dgp::DGP,
                        attributes::Dict = Dict("LogLevel" => 0),
                        startdf::DataFrame = DataFrame(
                            ℓ = Int64[],
                            d = Int64[],
                            j = Int64[],
                            k = Int64[],
                            start = Float64[]
                        ),
                        fixdf::DataFrame = DataFrame(
                            ℓ = Int64[],
                            d = Int64[],
                            j = Int64[],
                            k = Int64[],
                            fix = Float64[]
                        ))

    ############################################################################
    # Set up problem
    ############################################################################
    m = Model(() -> Clp.Optimizer())
    for (k, v) in attributes
        set_optimizer_attribute(m, k, v)
    end

    L = length(bases)
    J = [length(basis[d + 1].a) for basis in bases, d in 0:1]
    K = [length(basis[d + 1].b) for basis in bases, d in 0:1]
    @variable(m, θ[ℓ = 1:L, d = 0:1, j = 1:J[ℓ, d + 1], k = 1:K[ℓ, d + 1]])

    # set starting values
    @assert names(startdf) == ["ℓ", "d", "j", "k", "start"]
    for row in 1:nrow(startdf)
        ℓ, d, j, k, start = startdf[row, :]
        set_start_value(θ[ℓ, d, j, k], start)
    end

    # fix values
    @assert names(fixdf) == ["ℓ", "d", "j", "k", "fix"]
    for row in 1:nrow(fixdf)
        ℓ, d, j, k, fix = fixdf[row, :]
        if nrow(fixdf) == 1 # avoid warning about axis element
            row = [row]
        end
        @constraint(m, fixvar[row], θ[ℓ, d, j, k] == fix)
    end

    ############################################################################
    # Observational equivalence for each model
    ############################################################################

    # In MTW (2021), a saturated choice of IV-like specifications were used.
    # In MST (2018) and MT (2018), this assumption was not enforced.
    # To ensure backwards-compatibility, this assumption will be enforced
    # *unless* specified otherwise by the user.
    if (!haskey(assumptions, :saturated) || assumptions[:saturated])
        βₛ = compute_βₛ(dgp)
        Γₛ = compute_Γₛ(bases, dgp)
        @constraint(m, obseq[ℓ = 1:L, s = 1:length(βₛ)],
                    sum(θ[ℓ, d, j, k] * Γₛ[ℓ, d + 1][s, j, k]
                        for d = 0:1,
                            j = 1:J[ℓ, d + 1], k = 1:K[ℓ, d + 1]) == βₛ[s])
    end

    if (haskey(assumptions, :ivslope) && assumptions[:ivslope])
        βₛ = compute_βₛ(dgp, slist = "ivslope")
        Γₛ = compute_Γₛ(bases, dgp, slist = "ivslope")
        @constraint(m, ivslope[ℓ = 1:L, s = 1:length(βₛ)],
                    sum(θ[ℓ, d, j, k] * Γₛ[ℓ, d + 1][s, j, k]
                        for d = 0:1,
                            j = 1:J[ℓ, d + 1], k = 1:K[ℓ, d + 1]) == βₛ[s])
    end

    if (haskey(assumptions, :olsslope) && assumptions[:olsslope])
        βₛ = compute_βₛ(dgp, slist = "olsslope")
        Γₛ = compute_Γₛ(bases, dgp, slist = "olsslope")
        @constraint(m, olsslope[ℓ = 1:L, s = 1:length(βₛ)],
                    sum(θ[ℓ, d, j, k] * Γₛ[ℓ, d + 1][s, j, k]
                        for d = 0:1,
                            j = 1:J[ℓ, d + 1], k = 1:K[ℓ, d + 1]) == βₛ[s])
    end

    if (haskey(assumptions, :ivslopeind))
        βₛ = compute_βₛ(dgp, slist = "ivslopeind",
                        param = assumptions[:ivslopeind])
        Γₛ = compute_Γₛ(bases, dgp, slist = "ivslopeind",
                        param = assumptions[:ivslopeind])
        @constraint(m, ivslopeind[ℓ = 1:L, s = 1:length(βₛ)],
                    sum(θ[ℓ, d, j, k] * Γₛ[ℓ, d + 1][s, j, k]
                        for d = 0:1,
                            j = 1:J[ℓ, d + 1], k = 1:K[ℓ, d + 1]) == βₛ[s])
    end

    if (haskey(assumptions, :tslsslopeind) && assumptions[:tslsslopeind])
        βₛ = compute_βₛ(dgp, slist = "tslsslopeind")
        Γₛ = compute_Γₛ(bases, dgp, slist = "tslsslopeind")
        @constraint(m, tslsslopeind[ℓ = 1:L, s = 1:length(βₛ)],
                    sum(θ[ℓ, d, j, k] * Γₛ[ℓ, d + 1][s, j, k]
                        for d = 0:1,
                            j = 1:J[ℓ, d + 1], k = 1:K[ℓ, d + 1]) == βₛ[s])
    end

    if (haskey(assumptions, :wald))
        βₛ = compute_βₛ(dgp, slist = "wald", param = assumptions[:wald])
        Γₛ = compute_Γₛ(bases, dgp, slist = "wald", param = assumptions[:wald])
        @constraint(m, wald[ℓ = 1:L, s = 1:length(βₛ)],
                    sum(θ[ℓ, d, j, k] * Γₛ[ℓ, d + 1][s, j, k]
                    for d = 0:1, j = 1:J[ℓ, d + 1], k = 1:K[ℓ, d + 1]) == βₛ[s])
    end

    ############################################################################
    # Target parameters
    ############################################################################
    Γ⭑ = compute_Γ⭑(tp, bases, dgp)
    @objective(m, Min, sum(θ[ℓ, d, j, k] * Γ⭑[ℓ, d + 1][j, k]
                           for ℓ in 1:L, d in 0:1,
                               j in 1:J[ℓ, d + 1], k in 1:K[ℓ, d + 1]))

    ############################################################################
    # Add assumptions
    ############################################################################
    @variable(m, θz[ℓ = 1:L, d = 0:1,
                    k = 1:K[ℓ, d + 1], z = 1:size(dgp.suppZ, 1)])
    @constraint(m, θz_definition[ℓ = 1:L, d = 0:1,
                                 k = 1:K[ℓ, d + 1], z = 1:size(dgp.suppZ, 1)],
                sum(θ[ℓ, d, j, k] * bases[ℓ][d + 1].a[j](dgp.suppZ[z, :])
                    for j in 1:J[ℓ, d + 1]) == θz[ℓ, d, k, z])

    # Boundedness
    if (haskey(assumptions, :lb))
        @constraint(m, lb[ℓ = 1:L, d = 0:1,
                          k = 1:K[ℓ, d + 1], z = 1:size(dgp.suppZ, 1)],
                    θz[ℓ, d, k, z] >= assumptions[:lb])
    end
    if (haskey(assumptions, :ub))
        @constraint(m, ub[ℓ = 1:L, d = 0:1,
                          k = 1:K[ℓ, d + 1], z = 1:size(dgp.suppZ, 1)],
                    θz[ℓ, d, k, z] <= assumptions[:ub])
    end
    if (haskey(assumptions, :lb_diff) | haskey(assumptions, :ub_diff))
        for ℓ in 1:L @assert K[ℓ, 1] == K[ℓ, 2] end
    end
    if (haskey(assumptions, :lb_diff))
        @constraint(m, lb_diff[ℓ = 1:L,
                          k = 1:K[ℓ, 1], z = 1:size(dgp.suppZ, 1)],
                    θz[ℓ, 1, k, z] - θz[ℓ, 0, k, z] >= assumptions[:lb_diff])
    end
    if (haskey(assumptions, :ub_diff))
        @constraint(m, ub_diff[ℓ = 1:L,
                          k = 1:K[ℓ, 1], z = 1:size(dgp.suppZ, 1)],
                    θz[ℓ, 1, k, z] - θz[ℓ, 0, k, z] <= assumptions[:ub_diff])
    end

    # Mutual consistency
    if (haskey(assumptions, :mutually_consistent) &&
            assumptions[:mutually_consistent] && (L >= 2))
        Γₘ  = compute_Γₘ(bases, dgp)
        @constraint(m, mutually_consistent[ℓ = 1:L, d = 0:1, d′ = 0:1,
                                           z = 1:size(dgp.suppZ, 1)],
                    sum(θ[ℓ, d, j, k] * Γₘ[ℓ, d + 1][d′ + 1, z, j, k] for
                        j = 1:J[ℓ, d + 1], k = 1:K[ℓ, d + 1]) ==
                    sum(θ[1, d, j, k] * Γₘ[1, d + 1][d′ + 1, z, j, k] for
                        j = 1:J[1, d + 1], k = 1:K[1, d + 1]))
    end

    # Decreasing level
    if (haskey(assumptions, :decreasing_level) &&
            length(assumptions[:decreasing_level]) > 0)
        for (ℓ, d) in assumptions[:decreasing_level]
            con = @constraint(m, [k = 2:K[ℓ, d + 1], z = 1:size(dgp.suppZ, 1)],
                              θz[ℓ, d, k - 1, z] >= θz[ℓ, d, k, z])
        end
    end

    # Decreasing difference
    if (haskey(assumptions, :decreasing_difference) &&
            length(assumptions[:decreasing_difference]) > 0)
        for ℓ in assumptions[:decreasing_difference]
            @assert (K[ℓ, 1] == K[ℓ, 2])
            con = @constraint(m, [k = 2:K[ℓ, 1], z = 1:size(dgp.suppZ, 1)],
                              θz[ℓ, 1, k - 1, z] - θz[ℓ, 0, k - 1, z] >=
                              θz[ℓ, 1, k, z] - θz[ℓ, 0, k, z])
        end
    end

    ############################################################################
    # Solve problem
    ############################################################################
    function recover_optimal() # this dance is due to JuMP's container handling
        θjump = value.(θ)
        θarray = [[NaN for j in 1:J[ℓ, d + 1], k in 1:K[ℓ, d + 1]]
                  for ℓ in 1:L, d in 0:1]
        for key in keys(θjump.data)
            θarray[key[1], key[2] + 1][key[3], key[4]] = θjump.data[key]
        end
        return [(MTR(bases[ℓ][1], θarray[ℓ, 1]),
                 MTR(bases[ℓ][2], θarray[ℓ, 2])) for ℓ in 1:length(bases)]
    end

    results = Dict{Symbol, Any}(:lb => +Inf, :ub => -Inf)
    optimize!(m)
    if !(termination_status(m) in [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED])
        results[:lb] = objective_value(m)
        results[:mtr_lb] = recover_optimal()
        set_objective_sense(m, MOI.MAX_SENSE)
        optimize!(m)
        results[:ub] = objective_value(m)
        results[:mtr_ub] = recover_optimal()
    end

    return results
end
export compute_bounds
