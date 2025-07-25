################################################################################
# Allan Eduardo Feitosa
# 07 Apr 2025
# RBP and List-RBP Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/update_local_list.jl")
include("./RBP functions/update_global_list.jl")
include("./RBP functions/remove_residue!.jl")
include("./RBP functions/calc_residue.jl")

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
        residues::Vector{Float64},
        local_residues::Union{Vector{Float64},Nothing},
        local_coords::Union{Matrix{Int},Nothing},
        listsizes::Vector{Int},
        rbp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_reps

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        if residues[1] == 0.0
            init_list_RBP!(Lq,Lr,Nc,nothing,nothing,newLr,Factors,inlist,
                                            residues,coords,listsizes)
            if residues[1] == 0.0
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
        remove_residue!(limax,listsizes[1],residues,coords,inlist,1)

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
                        update_local_list!(residues,coords,local_residues,
                            local_coords,listsizes,inlist,li,ci,vj,residue)
                    end
                end
            end
        end
        update_global_list!(residues,coords,local_residues,local_coords,listsizes,
                                                                    inlist)
    end

    return rbp_not_converged
end

