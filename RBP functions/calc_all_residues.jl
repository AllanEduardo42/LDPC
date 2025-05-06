################################################################################
# Allan Eduardo Feitosa
# 3 Mar 2025
# Calculate all residues

include("calc_residue.jl")
include("add_residue.jl")

function calc_all_residues!(
    Lq::Matrix{<:AbstractFloat},
    Lr::Matrix{<:AbstractFloat},
    Nc::Vector{Vector{T}} where {T<:Integer},
    Lrj::Union{Vector{<:AbstractFloat},Nothing},
    signs::Union{Vector{Bool},Nothing},
    phi::Union{Vector{<:AbstractFloat},Nothing},
    Ms::Matrix{<:AbstractFloat},
    Factors::Matrix{<:AbstractFloat},        
    rbpmatrix::Matrix{<:Integer},
    residues::Vector{<:AbstractFloat},
    coords::Matrix{<:Integer},
    listsizes::Vector{<:Integer},
    relative::Bool
)
    @fastmath @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        # calculate the new check to node messages
        update_Lr!(Ms,Lq,ci,Nci,Lrj,signs,phi)
        # calculate the residues
        for vj in Nci
            li = LinearIndices(Lr)[ci,vj]
            residue = calc_residue(Ms,Lr,Factors,Lrj,Lq,li,relative)
            add_residue!(rbpmatrix,residues,coords,residue,li,ci,vj,listsizes[1])
        end
    end
end