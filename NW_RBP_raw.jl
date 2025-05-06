################################################################################
# Allan Eduardo Feitosa
# 06 mai 2025
# Nodewise RBP Algorithm with residual decaying factor (RAW)

include("./RBP functions/findmaxedge.jl")
include("./RBP functions/calc_residue.jl")

function
    NW_RBP_raw!(
        bitvector::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        Nv::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        ::Nothing,
        ::Nothing,
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
            Nvj = Nv[vj]
            for ci in Nvj
                if ci ≠ cimax
                    # 5) update Nv messages Lq[ci,vnmax]
                    Lq[ci,vj] = calc_Lq(Nvj,ci,vj,Lr,Lf)   
                    # 6) calculate alpha
                    Nci = Nc[ci]
                    maxresidue = 0.0
                    for vk in Nci
                        # if vk ≠ vj
                            li = LinearIndices(Lr)[ci,vk]
                            oldLr = Lr[li]
                            newLr = calc_Lr(Nci,ci,vk,Lq)
                            Ms[li] = newLr
                            residue = calc_residue(newLr,oldLr,Factors[ci])
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

    # 7) update bitvector
    for vj in eachindex(Nv)
        ci = Nv[vj][1]
        li = LinearIndices(Lr)[ci,vj]
        bitvector[vj] = signbit(Lr[li] + Lq[li])
    end

    return bp_not_converged
end