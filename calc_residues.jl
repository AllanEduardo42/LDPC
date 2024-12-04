################################################################################
# Allan Eduardo Feitosa
# 24 set 2024
# Specialized methods for the RBP algorithm

include("minsum.jl")
include("findmaxresidue.jl")

function
    calc_residues!(
        alpha::AbstractFloat,
        Residues::Union{Matrix{<:AbstractFloat},Nothing},
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        Ms::Matrix{<:AbstractFloat},
        Lr::Union{Matrix{<:AbstractFloat},Nothing},                      
        Lq::Matrix{<:AbstractFloat},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        vnmax::Integer,
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        listres1::Vector{<:AbstractFloat},
        listm1::Vector{<:Integer},
        listn1::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        listsize1::Integer,
        listsize2::Integer,
        inlist::Union{Matrix{<:Integer},Nothing}   
    )
    
    update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi,alpha)
    @inbounds for n in cn2vn[m]
        if n ≠ vnmax
            x = calc_residue(Ms,Factors,Lr,m,n)
            @fastmath if x != 0.0
                findmaxresidue!(Residues,m,n,x,listres1,listm1,listn1,listres2,
                listm2,listn2,listsize1,listsize2,inlist)
            end            
        end
    end

    if listres1[1] == 0 # no update
        listres1[1] = -1
    end

end

function 
    calc_residue(
        Ms::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        m::Integer,
        n::Integer
    )

    @inbounds i = LinearIndices(Ms)[m,n]

    @fastmath @inbounds x = Ms[i] - Lr[i]
    @fastmath if signbit(x)
        x = -x
    end
    @fastmath @inbounds x *= Factors[i]

    return x

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
        alpha::AbstractFloat,  
        Residues::Union{Matrix{<:AbstractFloat},Nothing},              
        Lq::Matrix{<:AbstractFloat},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        Ms::Matrix{<:AbstractFloat},
        listres1::Union{Vector{<:AbstractFloat},Nothing},
        listm1::Union{Vector{<:Integer},Nothing},
        listn1::Union{Vector{<:Integer},Nothing},
        listsize1::Integer,
        inlist::Union{Matrix{<:Integer},Nothing}
    )
    
    for m in eachindex(cn2vn)
        calc_residues!(alpha,Residues,nothing,Ms,nothing,Lq,Lrn,signs,phi,0,m,
            cn2vn,listres1,listm1,listn1,nothing,nothing,nothing,listsize1,0,
            inlist)
    end

end

