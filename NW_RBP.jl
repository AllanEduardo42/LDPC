################################################################################
# Allan Eduardo Feitosa
# 06 mai 2025
# Nodewise RBP Algorithm with residual decaying factor

include("./RBP functions/findmaxedge.jl")
include("./RBP functions/calc_residue.jl")

function
    NW_RBP!(
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
        M::Integer,
        Ms::Matrix{<:AbstractFloat},
        Factors::Vector{<:AbstractFloat},
        alpha::Vector{<:AbstractFloat},
        bp_not_converged::Bool
    )

    @fastmath @inbounds for m in 1:M

        # display("m = $m")

        # 1) Find largest alpha
        cimax = findmaxnode(alpha)
        # display(findmax(alpha))
        if cimax == 0
            bp_not_converged = false
            break # i.e., BP has converged
        end

        # 2) Set maximum alpha to zero
        alpha[cimax] = 0.0

        # 3) Decay the maximum alpha
        Factors[cimax] *= decayfactor

        for vj in Nc[cimax]
            # 4) update check to node messages Lr[cimax,vnmax]
            li = LinearIndices(Lr)[cimax,vj]
            Lr[li] = Ms[li]
            # 5) update Nv messages Lq[vj,ci] and bitvector[vj]
            Nvj = Nv[vj]
            bitvector[vj] = update_Lq!(Lq,Lr,Lf,vj,Nvj,Lrj)
            # 6) calculate residues
            for ci in Nvj
                if ci ≠ cimax
                    Nci = Nc[ci]    
                    # calculate the new check to node messages
                    update_Lr!(Ms,Lq,ci,Nci,Lrj,signs,phi)
                    # calculate alpha
                    maxresidue = 0.0
                    for vk in Nci
                        # if vk ≠ vj
                            li = LinearIndices(Lr)[ci,vk]
                            residue = calc_residue(Ms[li],Lr[li],Factors[ci])
                            if residue > maxresidue
                                maxresidue = residue
                            end
                        # end
                    end
                    alpha[ci] = maxresidue
                end
            end
        end
    end

    return bp_not_converged
end

