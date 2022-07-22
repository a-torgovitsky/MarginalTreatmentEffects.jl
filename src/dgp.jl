@with_kw struct DGP
    suppZ::Array{<:Real, 2}
    densZ::Array{<:Real, 1}
    pscore::Array{<:Real, 1}
    mtrs::Tuple{MTR, MTR}

    function DGP(suppZ, densZ, pscore, mtrs)
        @assert size(suppZ, 1) == length(densZ)
        @assert length(densZ) == length(pscore)
        @assert all(densZ .≥ 0)
        @assert sum(densZ) == 1

        new(suppZ, densZ, pscore, mtrs)
    end
end

function find_pscore(z, dgp::DGP)
    idx = findall([z == dgp.suppZ[i, :] for i in 1:size(dgp.suppZ, 1)])
    @assert length(idx) == 1 # find one and only one match
    return dgp.pscore[idx][1]
end

function find_density(z, dgp::DGP)
    idx = findall([z == dgp.suppZ[i, :] for i in 1:size(dgp.suppZ, 1)])
    @assert length(idx) == 1 # find one and only one match
    return dgp.densZ[idx][1]
end
