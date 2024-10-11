################################################################################
# Allan Eduardo Feitosa
# 24 set 2024
# Specialized methods for the RBP algorithm

include("minsum.jl")

function
    calc_residues!(
        Residues::Union{Matrix{<:AbstractFloat},Nothing},
        maxcoords::Vector{<:Integer},
        maxresidue::AbstractFloat,
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        Lr::Union{Matrix{<:AbstractFloat},Nothing},                      
        Lq::Matrix{<:AbstractFloat},
        signs::Vector{Bool},
        vnmax::Integer,
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        Ms::Matrix{<:AbstractFloat}       
    )
    
    minsum!(Lq,Ms,signs,m,cn2vn)    
    for n in cn2vn[m]
        if n ≠ vnmax
            y = __minsum_RBP!(Ms[m,n],Factors,Lr,m,n)
            maxresidue = _minsum_RBP!(Residues,maxcoords,maxresidue,m,n,y)
        end
    end

    return maxresidue
end

function 
    __minsum_RBP!(
        x::AbstractFloat,
        Factors::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        m::Integer,
        n::Integer
    )

    @fastmath @inbounds y = abs(x - Lr[m,n])*Factors[m,n]

    return y

end

function 
    __minsum_RBP!(
        x::AbstractFloat,
        Factors::Nothing,
        Lr::Nothing,
        m::Integer,
        n::Integer
    )

    @fastmath y = abs(x)

    return y

end

function 
    _minsum_RBP!(
        Residues::Matrix{<:AbstractFloat},
        maxcoords::Vector{Int},
        maxresidue::AbstractFloat,
        m::Integer,
        n::Integer,
        y::AbstractFloat
    )

    @inbounds Residues[m,n] = y

    return maxresidue
end

function 
    _minsum_RBP!(
        Residues::Nothing,
        maxcoords::Vector{Int},
        maxresidue::AbstractFloat,
        m::Integer,
        n::Integer,
        y::AbstractFloat
    )
    if @fastmath y > maxresidue
        maxresidue = y
        @inbounds maxcoords[1] = m
        @inbounds maxcoords[2] = n
    end

    return maxresidue
end

function
    init_residues!(    
        Residues::Union{Matrix{<:AbstractFloat},Nothing},      
        maxcoords::Vector{<:Integer},              
        Lq::Matrix{<:AbstractFloat},
        signs::Vector{Bool},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        Ms::Matrix{<:AbstractFloat}
    )
    
    maxresidue = 0.0
    for m in eachindex(cn2vn)
        maxresidue = minsum_RBP!(
            Residues,
            maxcoords,
            maxresidue,
            nothing,
            nothing,
            Lq,
            signs,
            0,
            m,
            cn2vn,
            Ms
        )
    end
end


function
    find_maxresidue_coords!(
        maxcoords::Vector{<:Integer},
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        samples::Vector{<:Integer},
        rng_sample::AbstractRNG
    )

    maxresidue = 0.0
    M = length(cn2vn)
    rand!(rng_sample,samples,1:M)
    for m in samples
        for n in cn2vn[m]
            @inbounds residue = Residues[m,n]
            if @fastmath residue > maxresidue
                maxresidue = residue
                @inbounds maxcoords[1] = m
                @inbounds maxcoords[2] = n
            end
        end
    end

    return maxresidue
end

function
    find_maxresidue_coords!(
        maxcoords::Vector{<:Integer},
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        samples::Nothing,
        rng_sample::Nothing
    )

    maxresidue = 0.0
    for m in eachindex(cn2vn)
        for n in cn2vn[m]
            @inbounds residue = Residues[m,n]
            if @fastmath residue > maxresidue
                maxresidue = residue
                @inbounds maxcoords[1] = m
                @inbounds maxcoords[2] = n
            end
        end
    end

    return maxresidue
end

