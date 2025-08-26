################################################################################
# Allan Eduardo Feitosa
# 05 Ago 2025
# RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/findmaxnode.jl")
include("./RBP functions/calc_residue.jl")
include("./List functions/find_list_pos.jl")

function
    List_VN_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        decayfactor::Float64,
        num_reps::Int,
        newLr::Matrix{Float64},
        Residues::Matrix{Float64},
        list::Vector{Float64},
        Factors::Matrix{Float64},
        coords::Matrix{Int},
        inlist::Matrix{Bool},
        listsize::Int,
        rbp_not_converged::Bool
    )
    
    @fastmath @inbounds for e in 1:num_reps

        # # display("e = $e")

        # # display(sum(inlist))

        # # display([residues coords'])

        # if max residue is equal to zero, refill the list
        if list[1] == 0.0
            for ci in eachindex(Nc)
                for vj in Nc[ci]
                    li = LinearIndices(Residues)[ci,vj]
                    residue = Residues[li]
                    # add residue to main list
                    add_to_list!(inlist,list,coords,residue,li,ci,vj,listsize)
                end
            end
            # if max residue is still zero, List-RBP has converged
            if list[1] == 0.0
                rbp_not_converged = false
                break
            end
        end
        
        vjmax = coords[2,1]
              
        Nvjmax = Nv[vjmax]
        for ci in Nvjmax
            li = LinearIndices(Lr)[ci,vjmax]
            # 2) Decay the RBP factor corresponding to the maximum residue
            Factors[li] *= decayfactor
            # 3) update check to node message Lr[cnmax,vnmax]
            Lr[li] = newLr[li]
            # 4) remove from list
            Residues[li] = 0.0
            if inlist[li]
                pos = find_list_pos(li,listsize,coords)
                remove_from_list!(li,listsize,list,coords,inlist,pos)
            end
        end

        Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
        bitvector[vjmax] = signbit(Ld)

        for ci in Nvjmax
            # 5) update Nv messages Lq[ci,vnmax]
            li = LinearIndices(Lq)[ci,vjmax]
            Lq[li] = tanh(0.5*(Ld - Lr[li]))
            # 6) calculate residues
            Nci = Nc[ci]
            for vj in Nci
                if vj ≠ vjmax
                    li = LinearIndices(Lr)[ci,vj]
                    newlr = calc_Lr(Nci,ci,vj,Lq)
                    li = LinearIndices(Lr)[ci,vj]
                    newLr[li] = newlr
                    residue = abs(newlr - Lr[li])*Factors[li]
                    Residues[li] = residue
                    if inlist[li]
                        pos = find_list_pos(li,listsize,coords)
                        remove_from_list!(li,listsize,list,coords,inlist,pos)
                    end
                    add_to_list!(inlist,list,coords,residue,li,ci,vj,listsize)
                end
            end
        end
    end

    return rbp_not_converged

end

# function
#     List_VN_RBP!(
#         bitvector::Vector{Bool},
#         Lq::Matrix{Float64},
#         Lr::Matrix{Float64},
#         Lf::Vector{Float64},
#         Nc::Vector{Vector{Int}},
#         Nv::Vector{Vector{Int}},
#         signs::Union{Vector{Bool},Nothing},
#         phi::Union{Vector{Float64},Nothing},
#         decayfactor::Float64,
#         num_reps::Int,
#         newLr::Matrix{Float64},
#         Factors::Matrix{Float64},
#         coords::Matrix{Int},
#         inlist::Matrix{Bool},
#         Residues::Matrix{Float64},
#         list::Vector{Float64},
#         local_list::Union{Vector{Float64},Nothing},
#         local_coords::Union{Matrix{Int},Nothing},
#         listsize::Int,
#         listsize2::Int,
#         rbp_not_converged::Bool
#     )
    
#     @fastmath @inbounds for e in 1:num_reps

#         # display("e = $e")

#         # display(sum(inlist))

#         # display([list coords'])

#         # if max residue is equal to zero, refill the list
#         if list[1] == 0.0
#             for ci in eachindex(Nc)
#                 for vj in Nc[ci]
#                     li = LinearIndices(Residues)[ci,vj]
#                     residue = Residues[li]
#                     # add residue to main list
#                     add_to_list!(inlist,list,coords,residue,li,ci,vj,listsize)
#                 end
#             end
#             # if max residue is still zero, List-RBP has converged
#             if list[1] == 0.0
#                 rbp_not_converged = false
#                 break
#             end
#         end
        
#         vjmax = coords[2,1]
              
#         Nvjmax = Nv[vjmax]
#         for ci in Nvjmax
#             li = LinearIndices(Lr)[ci,vjmax]
#             # 2) Decay the RBP factor corresponding to the maximum residue
#             Factors[li] *= decayfactor
#             # 3) update check to node message Lr[cnmax,vnmax]
#             Lr[li] = newLr[li]
#             # 4) remove from list
#             Residues[li] = 0.0
#             if inlist[li]
#                 # display("$li is on the list")
#                 pos = find_list_pos(li,listsize,coords)
#                 remove_from_list!(li,listsize,list,coords,inlist,pos)
#             end
#         end

#         # display([list coords'])

#         Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
#         bitvector[vjmax] = signbit(Ld)

#         # display([local_list local_coords'])

#         for ci in Nvjmax
#             # 5) update Nv messages Lq[ci,vnmax]
#             li = LinearIndices(Lq)[ci,vjmax]
#             Lq[li] = tanh(0.5*(Ld - Lr[li]))
#             # 6) calculate residues
#             Nci = Nc[ci]
#             for vj in Nci
#                 if vj ≠ vjmax
#                     li = LinearIndices(Lr)[ci,vj]
#                     newlr = calc_Lr(Nci,ci,vj,Lq)
#                     li = LinearIndices(Lr)[ci,vj]
#                     newLr[li] = newlr
#                     residue = abs(newlr - Lr[li])*Factors[li]
#                     Residues[li] = residue
#                     if inlist[li]
#                         # display("$li is on the local list")
#                         pos = find_list_pos(li,listsize,coords)
#                         remove_from_list!(li,listsize,list,coords,inlist,pos)
#                     end
#                     # add residue to local list
#                     add_to_list!(nothing,local_list,local_coords,residue,li,ci,vj,
#                                                                         listsize2)
#                 end
#             end
#         end
#         # display([local_list local_coords'])
#         # display([list coords'])
#         update_main_list!(list,coords,local_list,local_coords,listsize,listsize2,
#                                                                     inlist)

#         # display([list coords'])
#     end

#     return rbp_not_converged

# end
