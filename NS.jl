################################################################################
# Allan Eduardo Feitosa
# 07 Apr 2025
# RBP and List-RBP Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/findmaxedge.jl")
include("./RBP functions/update_lists.jl")
include("./RBP functions/remove_residue!.jl")

function
    NS!(
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
        Factors::Vector{<:AbstractFloat},
        residues::Vector{<:AbstractFloat},
        relative::Bool,
        rbp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        cnmax = findmaxnode(residues)
        if cnmax == 0
            rbp_not_converged  = false
            break # i.e., RBP has converged
        end
        residues[cnmax] = 0.0
        Factors[cnmax] *= decayfactor

        for n in cn2vn[cnmax]
            li = LinearIndices(Ms)[cnmax,n]
            Lr[li] = Ms[li]
            bitvector[n] = update_Lq!(Lq,Lr,Lf[n],n,vn2cn[n],Lrn)
            for m in vn2cn[n]
                maxresidue = 0.0
                if m ≠ cnmax
                    vns = cn2vn[m]    
                    # calculate the new check to node messages
                    update_Lr!(Ms,Lq,m,vns,Lrn,signs,phi)
                    # calculate the residues
                    for n2 in vns
                        li = LinearIndices(Lr)[m,n2]
                        residue = abs(Ms[li] - Lr[li])*Factors[m]
                        if residue > maxresidue
                            maxresidue = residue
                        end
                    end
                end
                residues[m] = maxresidue
            end
        end
    end

    return rbp_not_converged
end