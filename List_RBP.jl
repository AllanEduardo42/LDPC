################################################################################
# Allan Eduardo Feitosa
# 07 Apr 2025
# RBP and List-RBP Algorithm with residual decaying factor

include("./RBP functions/calc_residue.jl")
include("./List functions/update_main_list.jl")
include("./List functions/add_to_list.jl")
include("./List functions/remove_from_list.jl")
include("./List functions/find_list_pos.jl")

function
    List_RBP!(
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
        Factors::Matrix{Float64},
        coords::Matrix{Int},
        inlist::Matrix{Bool},
        Residues::Matrix{Float64},
        list::Vector{Float64},
        local_list::Union{Vector{Float64},Nothing},
        local_coords::Union{Matrix{Int},Nothing},
        listsize::Int,
        listsize2::Int,
        rbp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_reps

        # display("e = $e")

       # 1) Find largest residue and coordenates

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
        cimax = coords[1,1]
        vjmax = coords[2,1]
        limax = coords[3,1]

        # 2) Decay the RBP factor corresponding to the maximum residue
        Factors[limax] *= decayfactor
        # 3) update check to node message Lr[cnmax,vnmax]
        Lr[limax] = newLr[limax]
        # 4) set maximum residue to zero
        Residues[limax] = 0.0
        remove_from_list!(limax,listsize,list,coords,inlist,1)

        Nvjmax = Nv[vjmax]
        Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
        bitvector[vjmax] = signbit(Ld)

        for ci in Nvjmax
            if ci ≠ cimax
                # 5) update Nv messages Lq[ci,vnmax]
                li = LinearIndices(Lq)[ci,vjmax]
                Lq[li] = tanh(0.5*(Ld - Lr[li]))
                # 6) calculate residues
                Nci = Nc[ci]
                for vj in Nci
                    if vj ≠ vjmax
                        newlr = calc_Lr(Nci,ci,vj,Lq)
                        li = LinearIndices(Lr)[ci,vj]
                        newLr[li] = newlr
                        residue = abs(newlr - Lr[li])*Factors[li]
                        Residues[li] = residue
                        # if egde "li" is on the main list, find position and remove
                        if inlist[li]
                            pos = find_list_pos(li,listsize,coords)
                            remove_from_list!(li,listsize,list,coords,inlist,pos)
                        end
                        # add residue to local list
                        add_to_list!(nothing,local_list,local_coords,residue,li,
                                                                ci,vj,listsize2)
                    end
                end
            end
        end
        update_main_list!(list,coords,local_list,local_coords,listsize,listsize2,
                                                                         inlist)
    end

    return rbp_not_converged
end

