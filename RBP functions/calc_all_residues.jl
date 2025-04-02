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
    residues::Vector{<:AbstractFloat},
    M::Integer
)
    @fastmath @inbounds for m in 1:M 
        vns = cn2vn[m]
        # calculate the new check to node messages
        update_Lr!(Ms,Lq,m,vns,Lrn,signs,phi)
        # calculate the residues
        for n in vns
            residues[addressinv[m,n]] = calc_residue(Ms,Lr,Factors,Lrn,Lq,m,n)
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
    coords::Matrix{<:Integer},
    inlist::Union{Matrix{<:Integer},Nothing},
    M::Integer   
)
    @fastmath @inbounds for m in 1:M 
        vns = cn2vn[m]
        # calculate the new check to node messages
        update_Lr!(Ms,Lq,m,vns,Lrn,signs,phi)
        # calculate the residues
        for n in vns
            residue = calc_residue(Ms,Lr,Factors,Lrn,Lq,m,n)
            add_to_list!(inlist,listres1,coords,residue,m,n,listsizes[1])
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
    maxcoords::Vector{<:Integer},
    M::Integer
)
    @fastmath @inbounds for m in 1:M
        vns = cn2vn[m]
        # calculate the new check to node messages
        update_Lr!(Ms,Lq,m,vns,Lrn,signs,phi)
        # calculate the residues
        for n in vns
            residue = calc_residue(Ms,Lr,Factors,Lrn,Lq,m,n)
            if residue > max_residue[1]
                max_residue[2] = max_residue[1]
                max_residue[1] = residue 
                maxcoords[3] = maxcoords[1]
                maxcoords[4] = maxcoords[2]
                maxcoords[1] = m
                maxcoords[2] = n
            elseif residue > max_residue[2]
                max_residue[2] = residue
                maxcoords[3] = m
                maxcoords[4] = n
            end   
        end
    end
end