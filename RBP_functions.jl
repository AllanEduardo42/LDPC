################################################################################
# Allan Eduardo Feitosa
# 24 set 2024
# Specialized methods for the RBP algorithm

include("minsum.jl")

function
    minsum_RBP!(
        Residues::Nothing,
        maxcoords::Vector{<:Integer},
        maxresidue::AbstractFloat,
        Factors::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},                           
        Lq::Matrix{<:AbstractFloat},
        signs::Vector{Bool},
        vnmax::Integer,
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer}       
    )
    
    x = 0.0
    y = 0.0
    args = _minsum!(Lq,signs,m,cn2vn)
    for n in cn2vn[m]
        if n ≠ vnmax
            @inbounds x = __minsum!(n,signs[n],args...)
            @fastmath @inbounds y = abs(x - Lr[m,n])*Factors[m,n]
            if @fastmath y > maxresidue
                maxresidue = y
                @inbounds maxcoords[1] = m
                @inbounds maxcoords[2] = n
            end
        end
    end

    return maxresidue
end

function
    minsum_RBP!(
        Residues::Matrix{<:AbstractFloat},
        maxcoords::Vector{<:Integer},
        maxresidue::AbstractFloat,
        Factors::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},                           
        Lq::Matrix{<:AbstractFloat},
        signs::Vector{Bool},
        vnmax::Integer,
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )
    
    x = 0.0
    args = _minsum!(Lq,signs,m,cn2vn)
    for n in cn2vn[m]
        if n ≠ vnmax
            @inbounds x = __minsum!(n,signs[n],args...)
            @fastmath @inbounds Residues[m,n] = abs(x - Lr[m,n])*Factors[m,n]
        end
    end

    return maxresidue

end

function
    minsum_RBP_init!(    
        Residues::Nothing,      
        maxcoords::Vector{<:Integer},              
        Lq::Matrix{<:AbstractFloat},
        signs::Vector{Bool},
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )
    
    x = 0.0
    y = 0.0
    maxresidue = 0.0
    for m in eachindex(cn2vn)
        args = _minsum!(Lq,signs,m,cn2vn)
        for n in cn2vn[m]
            @inbounds x = __minsum!(n,signs[n],args...)
            @fastmath y = abs(x) 
            if @fastmath y > maxresidue
                maxresidue = y
                @inbounds maxcoords[1] = m
                @inbounds maxcoords[2] = n
            end
        end
    end
end

### specialized method for the RBP algorithm


function
    minsum_RBP_init!(     
        Residues::Matrix{<:AbstractFloat}, 
        maxcoords::Vector{<:Integer},                   
        Lq::Matrix{<:AbstractFloat},
        signs::Vector{<:Integer},
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )
    
    x = 0.0

    for m in eachindex(cn2vn)
        args = _minsum!(Lq,signs,m,cn2vn)
        for n in cn2vn[m]
            @inbounds x = __minsum!(n,signs[n],args...)
            @fastmath @inbounds Residues[m,n] = abs(x) 
        end
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

    maxresidue = 0
    M = size(Residues,1)
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
        rng_sample::AbstractRNG
    )

    maxresidue = 0
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