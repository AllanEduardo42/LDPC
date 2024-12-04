# RBP
function
    findmaxcoords!(
        listres1::Vector{<:AbstractFloat},
        listm1::Vector{<:Integer},
        listn1::Vector{<:Integer},
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        ::AbstractRNG
    )

    @fastmath @inbounds for m in eachindex(cn2vn)
        for n in cn2vn[m]
            residue = Residues[m,n]
            if residue > listres1[1]
                listres1[1] = residue
                listm1[1] = m
                listn1[1] = n
            end
        end
    end

end

# Random-RBP
function
    findmaxcoords!(
        listres1::Vector{<:AbstractFloat},
        listm1::Vector{<:Integer},
        listn1::Vector{<:Integer},
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
            if residue > listres1[1]
                listres1[1] = residue
                listm1[1] = m
                listn1[1] = n
            end
        end
    end
end

# Local-RBP and List-RBP
function 
    findmaxcoords!(
        ::Vector{<:AbstractFloat},
        ::Vector{<:Integer},
        ::Vector{<:Integer},
        ::Nothing,
        ::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        ::AbstractRNG
    )

end