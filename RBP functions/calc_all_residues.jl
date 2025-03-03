################################################################################
# Allan Eduardo Feitosa
# 3 Mar 2025
# Calculate all residues

include("calc_residue.jl")
include("add_to_list.jl")

function calc_all_residues!(
    Lq::Matrix{<:AbstractFloat},
    Lr::Matrix{<:AbstractFloat},
    cn2vn::Vector{Vector{T}} where {T<:Integer},
    Lrn::Union{Vector{<:AbstractFloat},Nothing},
    signs::Union{Vector{Bool},Nothing},
    phi::Union{Vector{<:AbstractFloat},Nothing},
    Ms::Matrix{<:AbstractFloat},
    Factors::Matrix{<:AbstractFloat},        
    addressinv::Matrix{<:Integer},
    residues::Vector{<:AbstractFloat}
)
    for m in 1:M 
        # calculate the new check to node messages
        update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)
        # calculate the residues
        @fastmath @inbounds for n in cn2vn[m]
            residue, index = calc_residue(Ms,Lr,Factors,Lrn,Lq,m,n)
            residues[addressinv[index]] = residue
        end
    end
end

function calc_all_residues_list!(
    Lq::Matrix{<:AbstractFloat},
    Lr::Matrix{<:AbstractFloat},
    cn2vn::Vector{Vector{T}} where {T<:Integer},
    Lrn::Union{Vector{<:AbstractFloat},Nothing},
    signs::Union{Vector{Bool},Nothing},
    phi::Union{Vector{<:AbstractFloat},Nothing},
    Ms::Matrix{<:AbstractFloat},
    Factors::Matrix{<:AbstractFloat},        
    listsizes::Vector{<:Integer},
    listres1::Vector{<:AbstractFloat},
    indices_res1::Vector{<:Integer},
    inlist::Union{Matrix{<:Integer},Nothing}    
)
    for m in 1:M 
        # calculate the new check to node messages
        update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)
        # calculate the residues
        for n in cn2vn[m]
            residue, li = calc_residue(Ms,Lr,Factors,Lrn,Lq,m,n)
            add_to_list!(inlist,listres1,indices_res1,residue,li,
                listsizes[1])
        end
    end
end

function calc_all_residues_local!(
    Lq::Matrix{<:AbstractFloat},
    Lr::Matrix{<:AbstractFloat},
    cn2vn::Vector{Vector{T}} where {T<:Integer},
    Lrn::Union{Vector{<:AbstractFloat},Nothing},
    signs::Union{Vector{Bool},Nothing},
    phi::Union{Vector{<:AbstractFloat},Nothing},
    Ms::Matrix{<:AbstractFloat},
    Factors::Matrix{<:AbstractFloat},        
    max_residue::Vector{<:AbstractFloat},
    max_indices::Vector{<:Integer}
)
    for m in 1:M
        # calculate the new check to node messages
        update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)
        # calculate the residues
        for n in cn2vn[m]
            residue, li = calc_residue(Ms,Lr,Factors,Lrn,Lq,m,n)
            if residue > max_residue[1]
                max_residue[2] = max_residue[1]
                max_residue[1] = residue 
                max_indices[2] = max_indices[1]
                max_indices[1] = li
            elseif residue > max_residue[2]
                max_residue[2] = residue
                max_indices[2] = li
            end   
        end
    end
end