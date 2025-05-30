################################################################################
# Allan Eduardo Feitosa
# 07 Apr 2025
# RBP and List-RBP Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/findmaxedge.jl")
include("./RBP functions/update_local_list.jl")
include("./RBP functions/update_global_list.jl")
include("./RBP functions/remove_residue!.jl")
include("./RBP functions/calc_residue.jl")

function
    RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        aux::Union{Vector{Float64},Nothing},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        decayfactor::Float64,
        num_edges::Int,
        newLr::Matrix{Float64},
        Factors::Matrix{Float64},
        coords::Matrix{Int},
        Edges::Matrix{<:Integer},
        residues::Vector{Float64},
        local_residues::Union{Vector{Float64},Nothing},
        local_coords::Union{Matrix{Int},Nothing},
        listsizes::Vector{Int},
        relative::Bool,
        bp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        max_edge, max_residue = findmaxedge(residues,local_residues)
        if max_residue == 0.0
            if max_edge == 0
                bp_not_converged = false
                break # i.e., BP has converged
            elseif max_edge == 1 # if list-RBP
                calc_all_residues!(Lq,Lr,Nc,aux,signs,phi,newLr,Factors,Edges,
                residues,coords,listsizes,relative)
                if residues[1] == 0.0
                    bp_not_converged = false
                    break
                end
                cimax = coords[1,1]
                vjmax = coords[2,1]
                limax = coords[3,1]
            end
        else
            cimax = coords[1,max_edge]
            vjmax = coords[2,max_edge]
            limax = coords[3,max_edge]
        end

        # 2) Decay the RBP factor corresponding to the maximum residue
        Factors[limax] *= decayfactor

        # 3) update check to node message Lr[cimax,vjmax]
        RBP_update_Lr!(limax,Lr,newLr,cimax,vjmax,Nc[cimax],Lq,aux,signs,phi)

        # 4) set maximum residue to zero
        remove_residue!(limax,listsizes[1],residues,coords,Edges,max_edge)

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
                            local_coords,listsizes,Edges,li,ci,vj,residue)
                    end
                end
            end
        end
        # if List-RBP: update list 1 
        update_global_list!(residues,coords,local_residues,local_coords,listsizes,
            Edges)
    end

    return bp_not_converged
end

# RAW
function
    RBP!(
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
        Edges::Matrix{<:Integer},
        residues::Vector{Float64},
        ::Nothing,
        ::Nothing,
        ::Vector{Int},
        relative::Bool,
        bp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        max_edge, _ = findmaxedge(residues,nothing)
        if max_edge == 0.0
            bp_not_converged = false
            break # i.e., BP has converged
        else
            cimax = coords[1,max_edge]
            vjmax = coords[2,max_edge]
            limax = coords[3,max_edge]
        end

        # 2) Decay the RBP factor corresponding to the maximum residue
        Factors[limax] *= decayfactor
        # 3) update check to node message Lr[cnmax,vnmax]
        Lr[limax] = newLr[limax]
        # 4) set maximum residue to zero
        residues[max_edge] = 0.0

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
                        residues[Edges[li]] = calc_residue_raw(newlr,Lr[li],
                                                Factors[li],relative,Lq[li])
                    end
                end
            end
        end
    end

    # 7) update bitvector
    for vj in eachindex(Nv)
        ci = Nv[vj][1]
        li = LinearIndices(Lr)[ci,vj]
        bitvector[vj] = signbit(Lr[li] + Lq[li])
    end

    return bp_not_converged
end

