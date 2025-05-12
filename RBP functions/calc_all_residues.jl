################################################################################
# Allan Eduardo Feitosa
# 3 Mar 2025
# Calculate all residues

include("calc_residue.jl")
include("add_residue.jl")

function
    calc_all_residues!(
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        aux::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{<:AbstractFloat},Nothing},
        newLr::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},        
        rbpmatrix::Matrix{<:Integer},
        residues::Vector{<:AbstractFloat},
        coords::Matrix{<:Integer},
        listsizes::Vector{<:Integer},
        relative::Bool
    )
    
    @fastmath @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        A, B, C, D = calc_ABCD!(aux,signs,phi,Lq,ci,Nci)
        for vj in Nci
            li = LinearIndices(Lr)[ci,vj]
            newlr = calc_Lr(A,B,C,D,vj,aux,signs,phi)
            newLr[li] = newlr
            residue = calc_residue(newLr[li],Lr[li],Factors[li],relative,Lq[li],aux)
            add_residue!(rbpmatrix,residues,coords,residue,li,ci,vj,listsizes[1])
        end
    end
end

# TANH
function 
    calc_all_residues!(
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        ::Nothing,
        ::Nothing,
        newLr::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},        
        rbpmatrix::Matrix{<:Integer},
        residues::Vector{<:AbstractFloat},
        coords::Matrix{<:Integer},
        listsizes::Vector{<:Integer},
        relative::Bool
    )

    @fastmath @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        for vj in Nci
            li = LinearIndices(Lr)[ci,vj]
            newlr = calc_Lr(Nci,ci,vj,Lq)
            newLr[li] = newlr
            residue = calc_residue(newLr[li],Lr[li],Factors[li],relative,Lq[li])
            add_residue!(rbpmatrix,residues,coords,residue,li,ci,vj,listsizes[1])
        end
    end
end