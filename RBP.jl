################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/findmaxedge.jl")
include("./RBP functions/update_lists.jl")
include("./RBP functions/remove_maxresidue!.jl")

function
    RBP!(
        bitvector::Vector{Bool},
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
        localresidues::Union{Vector{<:AbstractFloat},Nothing},
        localcoords::Union{Matrix{<:Integer},Nothing},
        listsizes::Vector{<:Integer},
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        max_edge , max_residue = findmaxedge(residues, localresidues)
        if max_residue == 0.0
            if max_edge == 0
                break # i.e., RBP has converged
            else
                calc_all_residues!(Lq,Lr,cn2vn,Lrn,signs,phi,Ms,Factors,rbpmatrix,
                residues,coords,listsizes)
                if residues[1] == 0.0
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
        remove_maxresidue!(limax,listsizes[1],residues,coords,rbpmatrix,max_edge)

        # 5) update vn2cn messages Lq[vnmax,m] and bitvector[vnmax]
        bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn[vnmax],Lrn)

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
                        residue = calc_residue(Ms,Lr,Factors,Lrn,Lq,li)
                        update_list2!(residues,coords,localresidues,localcoords,
                            listsizes,rbpmatrix,li,m,n,residue)
                    end
                end
            end
        end
        # update list 1 
        update_list1!(residues,coords,localresidues,localcoords,listsizes,
        rbpmatrix)

    end
end

