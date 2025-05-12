################################################################################
# Allan Eduardo Feitosa
# 08 Apr 2025
# Genius RBP Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/findmaxedge.jl")
include("./RBP functions/update_lists.jl")
include("./RBP functions/remove_residue!.jl")

function
    genius_RBP!(
        bitvector::Vector{Bool},
        cword::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{<:AbstractFloat},Nothing},
        decayfactor::AbstractFloat,
        num_edges::Integer,
        Ms::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},
        coords::Matrix{<:Integer},
        rbpmatrix::Matrix{<:Integer},
        residues::Vector{<:AbstractFloat},
        local_residues::Union{Vector{<:AbstractFloat},Nothing},
        local_coords::Union{Matrix{<:Integer},Nothing},
        listsizes::Vector{<:Integer},
        relative::Bool,
        rbp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        max_edge, max_residue = findmaxedge(residues,local_residues)
        if max_residue == 0.0
            if max_edge == 0
                rbp_not_converged = true
                break # i.e., RBP has converged
            else
                calc_all_residues!(Lq,Lr,cn2vn,Lrn,signs,phi,Ms,Factors,rbpmatrix,
                residues,coords,listsizes,relative)
                if residues[1] == 0.0
                    rbp_not_converged = true
                    break
                end
                cnmax = coords[1,1]
                vnmax = coords[2,1]
                limax = coords[3,1]
            end
        else
            cnmax = coords[1,max_edge]
            vnmax = coords[2,max_edge]
            limax = coords[3,max_edge]
        end

        # 2) Decay the RBP factor corresponding to the maximum residue
        Factors[limax] *= decayfactor

        # 3) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(limax,Lr,Ms,cnmax,vnmax,cn2vn[cnmax],Lq,Lrn,signs,phi)

        # 4) set maximum residue to zero
        remove_residue!(limax,listsizes[1],residues,coords,rbpmatrix,max_edge)

        # 5) update vn2cn messages Lq[vnmax,m] and bitvector[vnmax]
        bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn[vnmax],Lrn)

        TOTALBITERROR[e] = sum(bitvector .≠ cword)

        # 6) calculate residues
        for m in vn2cn[vnmax]
            if m ≠ cnmax
                vns = cn2vn[m]    
                # calculate the new check to node messages
                update_Lr!(Ms,Lq,m,vns,Lrn,signs,phi)
                # calculate the residues
                for n in vns
                    if n ≠ vnmax
                        li = LinearIndices(Lr)[m,n]
                        residue = calc_residue(Ms,Lr,Factors,Lrn,Lq,li,relative)
                        update_local_list!(residues,coords,local_residues,local_coords,
                            listsizes,rbpmatrix,li,m,n,residue)
                    end
                end
            end
        end
        # update list 1 
        update_global_list!(residues,coords,local_residues,local_coords,listsizes,
            rbpmatrix)

    end

    return rbp_not_converged
end

