################################################################################
# Allan Eduardo Feitosa
# 07 Apr 2025
# RBP and List-RBP Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/findmaxedge.jl")
include("./RBP functions/update_lists.jl")
include("./RBP functions/remove_residue!.jl")
include("./RBP functions/calc_residue.jl")

function
    RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        Nv::Vector{Vector{T}} where {T<:Integer},
        Lrj::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{<:AbstractFloat},Nothing},
        decayfactor::AbstractFloat,
        num_edges::Integer,
        newLr::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},
        coords::Matrix{<:Integer},
        rbpmatrix::Matrix{<:Integer},
        residues::Vector{<:AbstractFloat},
        local_residues::Union{Vector{<:AbstractFloat},Nothing},
        local_coords::Union{Matrix{<:Integer},Nothing},
        listsizes::Vector{<:Integer},
        relative::Bool,
        bp_not_converged::Bool
    )

    # count = 0 

    @fastmath @inbounds for e in 1:num_edges
    # @fastmath @inbounds while count < 200000

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        max_edge, max_residue = findmaxedge(residues,local_residues)
        if max_residue == 0.0
            if max_edge == 0
                bp_not_converged = false
                break # i.e., BP has converged
            elseif max_edge == -1 # if list-RBP
                calc_all_residues!(Lq,Lr,Nc,Lrj,signs,phi,newLr,Factors,rbpmatrix,
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
        RBP_update_Lr!(limax,Lr,newLr,cimax,vjmax,Nc[cimax],Lq,Lrj,signs,phi)

        # 4) set maximum residue to zero
        remove_residue!(limax,listsizes[1],residues,coords,rbpmatrix,max_edge)

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
                pLr, countzeros, vj0 = calc_pLr(Lq,ci,Nci,Lrj) 
                for vj in Nci
                    if vj ≠ vjmax
                        li = LinearIndices(Lr)[ci,vj]
                        # count += 1
                        newlr = fast_Lr(Lrj,pLr,countzeros,vj0,vj)
                        newLr[li]= newlr                                                    
                        residue = calc_residue(newlr,Lr[li],Factors[li],
                                                        relative,Lq[li],Lrj)
                        update_local_list!(residues,coords,local_residues,
                            local_coords,listsizes,rbpmatrix,li,ci,vj,residue)
                    end
                end
            end
        end
        # if List-RBP: update list 1 
        update_global_list!(residues,coords,local_residues,local_coords,listsizes,
            rbpmatrix)
    end

    # display(count)

    return bp_not_converged
end

