################################################################################
# Allan Eduardo Feitosa
# 3 Mar 2025
# Calculate all residues

include("calc_residue.jl")
include("add_residue.jl")

function calc_all_residues!(
    Lq::Matrix{<:AbstractFloat},
    Lr::Matrix{<:AbstractFloat},
    cn2vn::Vector{Vector{T}} where {T<:Integer},
    Lrn::Union{Vector{<:AbstractFloat},Nothing},
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
    @fastmath @inbounds for m in eachindex(cn2vn)
        vns = cn2vn[m]
        # calculate the new check to node messages
        update_Lr!(Ms,Lq,m,vns,Lrn,signs,phi)
        # calculate the residues
        for n in vns
            li = LinearIndices(Lr)[m,n]
            residue = calc_residue(Ms,Lr,Factors,Lrn,Lq,li,relative)
            add_residue!(rbpmatrix,residues,coords,residue,li,m,n,listsizes[1])
        end
    end
end

# function calc_all_residues_list!(
#     Lq::Matrix{<:AbstractFloat},
#     Lr::Matrix{<:AbstractFloat},
#     cn2vn::Vector{Vector{T}} where {T<:Integer},
#     Lrn::Union{Vector{<:AbstractFloat},Nothing},
#     signs::Union{Vector{Bool},Nothing},
#     phi::Union{Vector{<:AbstractFloat},Nothing},
#     Ms::Matrix{<:AbstractFloat},
#     Factors::Matrix{<:AbstractFloat},        
#     listsizes::Vector{<:Integer},
#     listres1::Vector{<:AbstractFloat},
#     coords::Matrix{<:Integer},
#     inlist::Union{Matrix{<:Integer},Nothing},
#     M::Integer   
# )
#     @fastmath @inbounds for m in 1:M 
#         vns = cn2vn[m]
#         # calculate the new check to node messages
#         update_Lr!(Ms,Lq,m,vns,Lrn,signs,phi)
#         # calculate the residues
#         for n in vns
#             li = LinearIndices(Lr)[m,n]
#             residue = calc_residue(Ms,Lr,Factors,Lrn,Lq,li)
            
#         end
#     end
# end

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
    relative::Bool
)
    @fastmath @inbounds for m in eachindex(cn2vn)
        vns = cn2vn[m]
        # calculate the new check to node messages
        update_Lr!(Ms,Lq,m,vns,Lrn,signs,phi)
        # calculate the residues
        for n in vns
            li = LinearIndices(Lr)[m,n]
            residue = calc_residue(Ms,Lr,Factors,Lrn,Lq,li,relative)
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