################################################################################
# Allan Eduardo Feitosa
# 3 Mar 2025
# Calculate local residues

include("calc_residue.jl")
include("update_lists.jl")

function calc_local_residues!(
    Lq::Matrix{<:AbstractFloat},
    Lr::Matrix{<:AbstractFloat},
    cn2vn::Vector{Vector{T}} where {T<:Integer},
    vn2cn::Vector{Vector{T}} where {T<:Integer},
    Lrn::Union{Vector{<:AbstractFloat},Nothing},
    signs::Union{Vector{Bool},Nothing},
    phi::Union{Vector{<:AbstractFloat},Nothing},
    Ms::Matrix{<:AbstractFloat},
    Factors::Matrix{<:AbstractFloat},        
    addressinv::Matrix{<:Integer},
    residues::Vector{<:AbstractFloat},
    cnmax::Integer,
    vnmax::Integer
)
    for m in vn2cn[vnmax]
        if m ≠ cnmax    
            # calculate the new check to node messages
            update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)
            # calculate the residues
            @fastmath @inbounds for n in cn2vn[m]
                if n ≠ vnmax
                    residue, li = calc_residue(Ms,Lr,Factors,Lrn,Lq,m,n)
                    residues[addressinv[li]] = residue
                end
            end
        end
    end
end

function calc_local_residues_list!(
    Lq::Matrix{<:AbstractFloat},
    Lr::Matrix{<:AbstractFloat},
    cn2vn::Vector{Vector{T}} where {T<:Integer},
    vn2cn::Vector{Vector{T}} where {T<:Integer},
    Lrn::Union{Vector{<:AbstractFloat},Nothing},
    signs::Union{Vector{Bool},Nothing},
    phi::Union{Vector{<:AbstractFloat},Nothing},
    Ms::Matrix{<:AbstractFloat},
    Factors::Matrix{<:AbstractFloat},        
    listsizes::Vector{<:Integer},
    listres1::Vector{<:AbstractFloat},
    indices_res1::Vector{<:Integer},
    listres2::Union{Vector{<:AbstractFloat},Nothing},
    indices_res2::Union{Vector{<:Integer},Nothing},
    inlist::Union{Matrix{<:Integer},Nothing},
    new_listsize2::Integer,
    cnmax::Integer,
    vnmax::Integer
)
    for m in vn2cn[vnmax]
        if m ≠ cnmax    
            # calculate the new check to node messages
            update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)
            # calculate the residues
            for n in cn2vn[m]
                if n ≠ vnmax
                    residue, li = calc_residue(Ms,Lr,Factors,Lrn,Lq,m,n)
                    new_listsize2 = update_list2!(listres1,indices_res1,
                                        listres2,indices_res2,listsizes,
                                        new_listsize2,inlist,li,residue)
                end
            end
        end
    end

    return new_listsize2
end

function calc_local_residues_local!(
    Lq::Matrix{<:AbstractFloat},
    Lr::Matrix{<:AbstractFloat},
    cn2vn::Vector{Vector{T}} where {T<:Integer},
    vn2cn::Vector{Vector{T}} where {T<:Integer},
    Lrn::Union{Vector{<:AbstractFloat},Nothing},
    signs::Union{Vector{Bool},Nothing},
    phi::Union{Vector{<:AbstractFloat},Nothing},
    Ms::Matrix{<:AbstractFloat},
    Factors::Matrix{<:AbstractFloat},
    max_residue::Vector{<:AbstractFloat},
    max_indices::Vector{<:Integer},
    cnmax::Integer,
    vnmax::Integer
)
    for m in vn2cn[vnmax]
        if m ≠ cnmax    
            # calculate the new check to node messages
            update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)
            # calculate the residues
            for n in cn2vn[m]
                if n ≠ vnmax
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
    end
end

