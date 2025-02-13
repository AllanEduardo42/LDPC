################################################################################
# Allan Eduardo Feitosa
# 24 set 2024
# Calculate the residues for the RBP algorithm

include("update_residue.jl")

function
    calc_residues!(
        addressinv::Union{Matrix{<:Integer},Nothing},
        residues::Union{Vector{<:AbstractFloat},Nothing},
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
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        listsize::Integer,
        listsize2::Integer,
        count_size::Integer,
        inlist::Union{Matrix{<:Integer},Nothing}   
    )
    # calculate the new check to node messages
    update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)

    # calculate the residues
    @fastmath @inbounds for n in cn2vn[m]
        if n ≠ vnmax
            l = LinearIndices(Ms)[m,n]
            x = calc_residue(Ms,Lr,l)
            count_size = update_residue!(addressinv,residues,m,n,l,x,Factors,listres,listm,listn,
                            listres2,listm2,listn2,listsize,listsize2,count_size,inlist)
        end
    end

    return count_size

end

function 
    calc_residue(
        Ms::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        l::Integer
    )

    @fastmath @inbounds x = Ms[l] - Lr[l]
    @fastmath if signbit(x)
        x = -x
    end

    return x

end

# for initialization (Lr[m,n] = 0.0 and Factors[m,n] = 1.0 ∀m,n)
function 
    calc_residue(
        Ms::Matrix{<:AbstractFloat},
        ::Nothing,
        l::Integer
    )

    @fastmath @inbounds return abs(Ms[l])

end

