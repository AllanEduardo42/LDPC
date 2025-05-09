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
        Lrj::Vector{<:AbstractFloat},
        decayfactor::AbstractFloat,
        M::Integer,
        newLr::Matrix{<:AbstractFloat},
        Factors::Vector{<:AbstractFloat},
        alpha::Vector{<:AbstractFloat},
        bp_not_converged::Bool
    )

    # count = 0

    @fastmath @inbounds for m in 1:M
    # @fastmath @inbounds while count < 200000

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
            Lr[li] = newLr[li]
            # 5) update Nv messages Lq[vj,ci] and bitvector[vj]
            Nvj = Nv[vj]
            Ld = calc_Ld(vj,Nvj,Lf,Lr)
            bitvector[vj] = signbit(Ld)
            # 6) calculate residues
            for ci in Nvj
                if ci ≠ cimax
                    # 6) update Nv messages Lq[ci,vj]
                    li = LinearIndices(Lq)[ci,vj]
                    Lq[li] = Ld - Lr[li]
                    Nci = Nc[ci]    
                    # calculate the new check to node messages
                    pLr, countzeros, vj0 = calc_pLr(Lq,ci,Nci,Lrj) 
                    # calculate alpha
                    maxresidue = 0.0
                    for vk in Nci
                        # if vk ≠ vj
                        # count += 1
                            li = LinearIndices(Lr)[ci,vk]
                            newlr = fast_Lr(Lrj,pLr,countzeros,vj0,vk)
                            newLr[li]= newlr
                            residue = calc_residue(newlr,Lr[li],Factors[ci])
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

    # display(count)

    return bp_not_converged
end

# RAW
function
    NW_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        Nv::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        decayfactor::AbstractFloat,
        M::Integer,
        newLr::Matrix{<:AbstractFloat},
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
            Lr[li] = newLr[li]
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
                            newlr = calc_Lr(Nci,ci,vk,Lq)
                            newLr[li] = newlr
                            residue = calc_residue(newlr,oldLr,Factors[ci])
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

