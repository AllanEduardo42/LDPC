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

