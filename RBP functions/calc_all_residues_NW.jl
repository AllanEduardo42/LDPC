################################################################################
# Allan Eduardo Feitosa
# 12 May 2025
# Calculate all alpha for the NW-RBP

include("calc_residue.jl")

function
    calc_all_residues_NW!(
        Lq::Matrix{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        aux::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{<:AbstractFloat},Nothing},
        newLr::Matrix{<:AbstractFloat},
        alpha::Vector{<:AbstractFloat},
    )
    
    @fastmath @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        A, B, C, D = calc_ABCD!(aux,signs,phi,Lq,ci,Nci)
        maxresidue = 0.0
        for vj in Nci
            li = LinearIndices(newLr)[ci,vj]
            newlr = calc_Lr(A,B,C,D,vj,aux,signs,phi)
            newLr[li] = newlr
            residue = calc_residue(newlr,0.0,1.0)
            if residue > maxresidue
                maxresidue = residue
            end
        end
        alpha[ci] = maxresidue
    end
end

# RAW
function
    calc_all_residues_NW!(
        Lq::Matrix{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        ::Nothing,
        ::Nothing,
        newLr::Matrix{<:AbstractFloat},
        alpha::Vector{<:AbstractFloat},
    )
    
    @fastmath @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        maxresidue = 0.0
        for vj in Nci
            li = LinearIndices(newLr)[ci,vj]
            newlr = calc_Lr(Nci,ci,vj,Lq)
            newLr[li] = newlr
            residue = calc_residue(newlr,0.0,1.0)
            if residue > maxresidue
                maxresidue = residue
            end
        end
        alpha[ci] = maxresidue
    end
end