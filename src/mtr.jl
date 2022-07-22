# An MTR is an MTRBasis together with a coefficient vector
#   note this is an MTR for one treatment arm (d = 0 or d = 1)
#   so in practice we are carrying around 2 element tuples of MTR objects
#
# θ[j,k] corresponds to a_{j}(z)b_{k}(u)
@with_kw struct MTR
    basis::MTRBasis
    θ::Matrix{<:Real} # coefficient vector, arranged in a matrix

    function MTR(basis, θ)
        @assert length(basis.a) == size(θ, 1)
        @assert length(basis.b) == size(θ, 2)
        new(basis, θ)
    end
end

function evaluate_mtr(mtr::MTR, ev::DataFrame)
    result = fill(NaN, nrow(ev))
    for (i,r) in enumerate(eachrow(ev))
        result[i] = sum([aj(r[:z]) * bk(r[:u]) for
                            aj in mtr.basis.a, bk in mtr.basis.b]
                        .* mtr.θ)
    end
    return result
end

function evaluate_mtr(mtrs::Tuple{MTR, MTR}, ev::DataFrame)
    result = fill(NaN, nrow(ev))
    for d in 0:1
        result[(ev.d == d)] = evaluate_mtr(mtrs[d+1], ev[(ev.d == d)])
    end
    return result
end
