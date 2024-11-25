################################################################################
# Allan Eduardo Feitosa
# 24 set 2024
# Specialized methods for the RBP algorithm

include("minsum.jl")
include("findmaxresidue.jl")

function
    calc_residues!(
        Residues::Union{Matrix{<:AbstractFloat},Nothing},
        maxcoords::Vector{<:Integer},
        maxresidue::AbstractFloat,
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        Ms::Matrix{<:AbstractFloat},
        Lr::Union{Matrix{<:AbstractFloat},Nothing},                      
        Lq::Matrix{<:AbstractFloat},
        signs::Vector{Bool},
        vnmax::Integer,
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        listres1::Union{Vector{<:AbstractFloat},Nothing},
        listadd1::Union{Matrix{<:Integer},Nothing},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listadd2::Union{Matrix{<:Integer},Nothing},
        listaddinv1::Union{Matrix{<:Integer},Nothing},
        listsize1::Integer,
        listsize2::Integer,
        inlist1::Union{Matrix{<:Integer},Nothing}   
    )
    
    minsum!(Lq,Ms,signs,m,cn2vn)
    x = 0.0    
    @fastmath @inbounds for n in cn2vn[m]
        if n ≠ vnmax
            x = calc_residue(Ms,Factors,Lr,m,n)
            if x != 0.0
                maxresidue = findmaxresidue!(Residues,maxcoords,maxresidue,m,n,x,
                    listres1,listadd1,listres2,listadd2,listaddinv1,listsize1,
                    listsize2,inlist1)
            end            
        end
    end

    return maxresidue
end

function 
    calc_residue(
        Ms::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        m::Integer,
        n::Integer
    )

    @fastmath @inbounds return abs(Ms[m,n] - Lr[m,n])*Factors[m,n]

end

# for initialization (Lr[m,n] = 0.0 and Factors[m,n] = 1.0 ∀m,n)
function 
    calc_residue(
        Ms::Matrix{<:AbstractFloat},
        ::Nothing,
        ::Nothing,
        m::Integer,
        n::Integer
    )

    @fastmath @inbounds return abs(Ms[m,n])

end

# RBP and Random-RBP
function 
    findmaxresidue!(
        Residues::Matrix{<:AbstractFloat},
        maxcoords::Vector{<:Integer},
        maxresidue::AbstractFloat,
        m::Integer,
        n::Integer,
        x::AbstractFloat,
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Integer,
        ::Nothing
    )

    @inbounds Residues[m,n] = x

    if @fastmath x > maxresidue
        maxresidue = x
        @inbounds maxcoords[1] = m
        @inbounds maxcoords[2] = n
    end

    return maxresidue
end

# Local-RBP
function 
    findmaxresidue!(
        ::Nothing,
        maxcoords::Vector{<:Integer},
        maxresidue::AbstractFloat,
        m::Integer,
        n::Integer,
        x::AbstractFloat,
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Integer,
        ::Nothing
    )
    if @fastmath x > maxresidue
        maxresidue = x
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
        Ms::Matrix{<:AbstractFloat},
        listres1::Union{Vector{<:AbstractFloat},Nothing},
        listadd1::Union{Matrix{<:Integer},Nothing},
        listaddinv1::Union{Matrix{<:Integer},Nothing},
        listsize1::Integer,
        inlist1::Union{Matrix{<:Integer},Nothing}
    )
    
    maxresidue = 0.0
    for m in eachindex(cn2vn)
        maxresidue = calc_residues!(
            Residues,
            maxcoords,
            maxresidue,
            nothing,
            Ms,
            nothing,
            Lq,
            signs,
            0,
            m,
            cn2vn,
            listres1,
            listadd1,
            nothing,
            nothing,
            listaddinv1,
            listsize1,
            0,
            inlist1
        )
    end

    return maxresidue
end

# RBP
function
    find_maxresidue_coords!(
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
    find_maxresidue_coords!(
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
    find_maxresidue_coords!(
        maxresidue::AbstractFloat,
        ::Vector{<:Integer},
        ::Nothing,
        ::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        ::AbstractRNG
    )

    return maxresidue
end

