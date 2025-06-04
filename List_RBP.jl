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
        aux::Vector{Float64},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        decayfactor::Float64,
        num_edges::Int,
        newLr::Matrix{Float64},
        Factors::Matrix{Float64},
        coords::Matrix{Int},
        inlist::Matrix{Bool},
        residues::Vector{Float64},
        local_residues::Union{Vector{Float64},Nothing},
        local_coords::Union{Matrix{Int},Nothing},
        listsizes::Vector{Int},
        relative::Bool,
        rbp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        if residues[1] == 0.0
            init_list_RBP!(Lq,Lr,Nc,aux,signs,phi,newLr,Factors,inlist,
                                            residues,coords,listsizes,relative)
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

        # 3) update check to node message Lr[cimax,vjmax]
        RBP_update_Lr!(limax,Lr,newLr,cimax,vjmax,Nc[cimax],Lq,aux,signs,phi)

        # 4) set maximum residue to zero
        remove_residue!(limax,listsizes[1],residues,coords,inlist,1)

        # 5) Calculate Ld of vjmax and bitvector[vjmax]
        Nvjmax = Nv[vjmax]
        Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
        bitvector[vjmax] = signbit(Ld)

        for ci in Nvjmax
            if ci ≠ cimax
                # 6) update Nv messages Lq[ci,vjmax]
                li = LinearIndices(Lq)[ci,vjmax]
                Lq[li] = Ld - Lr[li]
                # 7) calculate residues
                Nci = Nc[ci]    
                A, B, C, D = calc_ABCD!(aux,signs,phi,Lq,ci,Nci)
                for vj in Nci
                    if vj ≠ vjmax
                        newlr = calc_Lr(A,B,C,D,vj,aux,signs,phi)
                        li = LinearIndices(Lr)[ci,vj]
                        newLr[li] = newlr                                              
                        residue = calc_residue(newlr,Lr[li],Factors[li],
                                                        relative,Lq[li])
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

# RAW
function
    List_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        ::Nothing,
        ::Nothing,
        ::Nothing,
        decayfactor::Float64,
        num_edges::Int,
        newLr::Matrix{Float64},
        Factors::Matrix{Float64},
        coords::Matrix{Int},
        inlist::Matrix{Bool},
        residues::Vector{Float64},
        local_residues::Union{Vector{Float64},Nothing},
        local_coords::Union{Matrix{Int},Nothing},
        listsizes::Vector{Int},
        relative::Bool,
        rbp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        if residues[1] == 0.0
            init_list_RBP!(Lq,Lr,Nc,nothing,nothing,nothing,newLr,Factors,inlist,
                                            residues,coords,listsizes,relative)
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
        for ci in Nvjmax
            if ci ≠ cimax
                # 5) update Nv messages Lq[ci,vnmax]
                Lq[ci,vjmax] = calc_Lq(Nvjmax,ci,vjmax,Lr,Lf)
                # 6) calculate residues
                Nci = Nc[ci]
                for vj in Nci
                    if vj ≠ vjmax
                        newlr = calc_Lr(Nci,ci,vj,Lq)
                        li = LinearIndices(Lr)[ci,vj]
                        newLr[li] = newlr
                        residue = calc_residue_raw(newlr,Lr[li],Factors[li],
                                                        relative,Lq[li])
                        update_local_list!(residues,coords,local_residues,
                            local_coords,listsizes,inlist,li,ci,vj,residue)
                    end
                end
            end
        end
    end
    update_global_list!(residues,coords,local_residues,local_coords,listsizes,
                                                                    inlist)

    # 7) update bitvector
    for vj in eachindex(Nv)
        ci = Nv[vj][1]
        li = LinearIndices(Lr)[ci,vj]
        bitvector[vj] = signbit(Lr[li] + Lq[li])
    end

    return rbp_not_converged
end

