# An MTRBasis for d = 0 or 1 has functions:
#   \{a_{j}(z)\}_{j=1}^{J}
#   \{b_{k}(u)\}_{k=1}^{K}
#   \{\int_{u}^{v} b_{k}(u')\, du\}_{k=1}^{K}
# The basis is then
#   \{a_{j}(z)b_{k}(u) : j = 1,...,J, k = 1...K\}
#
# The rationale for this formulation is:
#   - Allows me to integrate with respect to u fairly easily
#   which would be difficult with say just a single set of functions b_{k}(u,z)
#   - Allows me to keep track of the coefficients by (d,j,k), which is
#   important for imposing certain types of constraints.
@with_kw struct MTRBasis
    a::Vector{Function}
    b::Vector{Function}
    ib::Vector{Function}

    function MTRBasis(a,b,ib)
        @assert length(b) == length(ib)
        new(a, b, ib)
    end
end

# Take multiple MTRBasis and concatenate them
function MTRBasis(bases::Vector{MTRBasis})
    anew = vcat([basis.a for basis in bases]...)
    bnew = vcat([basis.b for basis in bases]...)
    ibnew = vcat([basis.ib for basis in bases]...)
    MTRBasis(anew, bnew, ibnew)
end

function bernstein_basis(K)
    a = [(z) -> 1]
    b = [u -> bernstein_polynomial(u, k, K) for k in 0:K]
    ib = [(u, v) -> integrate_bernstein_polynomial(u, v, k, K) for k in 0:K]
    MTRBasis(a, b, ib)
end

function interacted_bernstein_basis(K; notℓ = 1)
    @assert notℓ in [1, 2]

    basis = bernstein_basis(K)
    a = [z -> 1, z -> z[notℓ]]
    b = vcat(basis.b)
    ib = vcat(basis.ib)

    MTRBasis(a, b, ib)
end

################################################################################
# The kth Bernstein polynomial of degree K, evaluated at u
#   for k = 0,1,...,K
# The formula is:
#   B_{k,K}(u) = nchoosek(K, k) u^(k)*(u^(K-k))
################################################################################
function bernstein_polynomial(u, k::Integer, K::Integer)
    if ((u < 0) | (u > 1))
        throw(DomainError("Should have u in [0,1]."))
    end
    bernstein_errorcheck(k, K)
    return binomial(K,k) * (u.^k) .* ((1 .- u) .^ (K-k));
end

function bernstein_errorcheck(k::Integer, K::Integer)
    if ((k < 0) | (k > K))
        throw(DomainError("Should have k in 0,1,...,K."))
    end
    if (K < 0)
        throw(DomainError("Should have K >= 0."))
    end
end

################################################################################
# Integrate B_{k,K}(u) over u
#   where B_{k,K} is a Bernstein polynomial
#
# Input:
#   u: An N-vector of lower integration points
#   v: An N-vector of upper integration points
#   k, K: kth Bernstein polynomial of degree K
#       for k = 0,1,...,K
# Output:
#   An N x 1 column vector with ith component:
#       \int_{u(i)}^{v(i)} b_{k,K}(u)\, du
################################################################################
function integrate_bernstein_polynomial(u::Array{<:Real, 1},
                                        v::Array{<:Real, 1},
                                        k::Integer, K::Integer)
    @assert size(u) == size(v)
    @assert all(u .>= 0) & all(u .<= 1)
    @assert all(v .>= 0) * all(v .<= 1)
    @assert all(u .<= v)
    bernstein_errorcheck(k, K)

    BI = zeros(size(u));
    for i = k:1:K
        BI =  BI +
              binomial(K, i) *
              binomial(i, k) *
              (-1)^(i - k) .*
              (v.^(i+1) - u.^(i+1)) ./ (i+1);
    end
    return BI
end

function integrate_bernstein_polynomial(u::Real, v::Real,
                                        k::Integer, K::Integer)
    integrate_bernstein_polynomial([u], [v], k, K)[1]
end

function constantspline_basis(knots::Vector{<:Real})
    a = [(z) -> 1]
    @assert all(0 .<= knots .<= 1)
    @assert (0 in knots) & (1 in knots)
    unique!(knots)
    sort!(knots)
    # The knots partition [0, 1] into half-open partitions of the form
    # [knots[k-1], knots[k]), except the final partition is closed on both ends.
    in_partition = (u, k) -> ifelse(
        k != last(knots),
        Int.(knots[k-1] .<= u .< knots[k]),
        Int.(knots[k-1] .<= u .<= knots[k])
    )
    b = [u -> in_partition(u, k)
         for k in 2:length(knots)]
    ib = [(u, v) -> max(0, (min(v, knots[k]) - max(u, knots[k-1])))
          for k in 2:length(knots)]
    MTRBasis(a, b, ib)
end

function interacted_constantspline_basis(knots::Vector{<:Real}; notℓ = 1)
    @assert notℓ in [1, 2]

    basis = constantspline_basis(knots)
    a = [z -> 1, z -> z[notℓ]]
    b = vcat(basis.b)
    ib = vcat(basis.ib)

    MTRBasis(a, b, ib)
end

function switching_constantspline_basis(knots::Vector{<:Real}; notℓ = 1)
    @assert notℓ in [1, 2]

    basis = interacted_constantspline_basis(knots)
    a = [z -> (1 - z[notℓ]), z -> z[notℓ]]
    MTRBasis(a, basis.b, basis.ib)
end
