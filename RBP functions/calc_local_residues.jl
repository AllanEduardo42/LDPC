################################################################################
# Allan Eduardo Feitosa
# 3 Mar 2025
# Calculate local residues

include("calc_residue.jl")
# include("update_lists.jl")

# function calc_local_residues!(
#     Lq::Matrix{<:AbstractFloat},
#     Lr::Matrix{<:AbstractFloat},
#     cn2vn::Vector{Vector{T}} where {T<:Integer},
#     cns::Vector{<:Integer},
#     Lrn::Union{Vector{<:AbstractFloat},Nothing},
#     signs::Union{Vector{Bool},Nothing},
#     phi::Union{Vector{<:AbstractFloat},Nothing},
#     Ms::Matrix{<:AbstractFloat},
#     Factors::Matrix{<:AbstractFloat},        
#     addressinv::Matrix{<:Integer},
#     residues::Vector{<:AbstractFloat},
#     cnmax::Integer,
#     vnmax::Integer
# )
#     @fastmath @inbounds for m in cns
#         if m ≠ cnmax
#             vns = cn2vn[m]    
#             # calculate the new check to node messages
#             update_Lr!(Ms,Lq,m,vns,Lrn,signs,phi)
#             # calculate the residues
#             for n in vns
#                 if n ≠ vnmax
#                     li = LinearIndices(Lr)[m,n]
#                     residues[addressinv[li]] = calc_residue(Ms,Lr,Factors,Lrn,Lq,li)
#                 end
#             end
#         end
#     end
# end

# function calc_local_residues_list!(
#     Lq::Matrix{<:AbstractFloat},
#     Lr::Matrix{<:AbstractFloat},
#     cn2vn::Vector{Vector{T}} where {T<:Integer},
#     cns::Vector{<:Integer},
#     Lrn::Union{Vector{<:AbstractFloat},Nothing},
#     signs::Union{Vector{Bool},Nothing},
#     phi::Union{Vector{<:AbstractFloat},Nothing},
#     Ms::Matrix{<:AbstractFloat},
#     Factors::Matrix{<:AbstractFloat},        
#     listsizes::Vector{<:Integer},
#     listres1::Vector{<:AbstractFloat},
#     coords1::Matrix{<:Integer},
#     listres2::Union{Vector{<:AbstractFloat},Nothing},
#     coords2::Union{Matrix{<:Integer},Nothing},
#     inlist::Union{Matrix{<:Integer},Nothing},
#     new_listsize2::Integer,
#     cnmax::Integer,
#     vnmax::Integer
# )
#     @fastmath @inbounds for m in cns
#         if m ≠ cnmax
#             vns = cn2vn[m] 
#             # calculate the new check to node messages
#             update_Lr!(Ms,Lq,m,vns,Lrn,signs,phi)
#             # calculate the residues
#             for n in vns
#                 if n ≠ vnmax
#                     li = LinearIndices(Lr)[m,n]
#                     residue = calc_residue(Ms,Lr,Factors,Lrn,Lq,li)
#                     new_listsize2 = update_list2!(listres1,coords1,
#                                         listres2,coords2,listsizes,
#                                         new_listsize2,inlist,li,m,n,residue)
#                 end
#             end
#         end
#     end

#     return new_listsize2
# end

function calc_local_residues_local!(
    Lq::Matrix{<:AbstractFloat},
    Lr::Matrix{<:AbstractFloat},
    cn2vn::Vector{Vector{T}} where {T<:Integer},
    cns::Vector{<:Integer},
    Lrn::Union{Vector{<:AbstractFloat},Nothing},
    signs::Union{Vector{Bool},Nothing},
    phi::Union{Vector{<:AbstractFloat},Nothing},
    Ms::Matrix{<:AbstractFloat},
    Factors::Matrix{<:AbstractFloat},
    max_residue::Vector{<:AbstractFloat},
    maxcoords::Vector{<:Integer},
    cnmax::Integer,
    vnmax::Integer
)
    @fastmath @inbounds for m in cns
        if m ≠ cnmax    
            vns = cn2vn[m] 
            # calculate the new check to node messages
            update_Lr!(Ms,Lq,m,vns,Lrn,signs,phi)
            # calculate the residues
            for n in vns
                if n ≠ vnmax
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
    end
end

