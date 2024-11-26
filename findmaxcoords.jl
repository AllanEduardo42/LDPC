# RBP
function
    findmaxcoords!(
        maxresidue::AbstractFloat,
        maxcoords::Vector{<:Integer},
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        ::Nothing
    )

    @fastmath @inbounds for m in eachindex(cn2vn)
        for n in cn2vn[m]
            residue = Residues[m,n]
            if residue > maxresidue
                maxresidue = residue
                maxcoords[1] = m
                maxcoords[2] = n
            end
        end
    end

    return maxresidue
end

# Random-RBP
function
    findmaxcoords!(
        maxresidue::AbstractFloat,
        maxcoords::Vector{<:Integer},
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        samples::Vector{<:Integer},
        rng_sample::AbstractRNG
    )

    M = length(cn2vn)
    rand!(rng_sample,samples,1:M)
    @fastmath @inbounds for m in samples
        for n in cn2vn[m]
            residue = Residues[m,n]
            if residue > maxresidue
                maxresidue = residue
                maxcoords[1] = m
                maxcoords[2] = n
            end
        end
    end

    return maxresidue
end

# Local-RBP and List-RBP
function 
    findmaxcoords!(
        maxresidue::AbstractFloat,
        ::Vector{<:Integer},
        ::Nothing,
        ::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        ::AbstractRNG
    )

    return maxresidue
end