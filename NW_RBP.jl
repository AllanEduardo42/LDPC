################################################################################
# Allan Eduardo Feitosa
# 06 mai 2025
# Nodewise RBP Algorithm with residual decaying factor

include("./RBP functions/findmaxnode.jl")
include("./RBP functions/calc_residue.jl")
include("./RBP functions/RBP_update_Lr.jl")

function
    NW_RBP!(
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
        Factors::Vector{Float64},
        alpha::Vector{Float64},
        bp_not_converged::Bool,
        raw::Bool
    )

    @fastmath @inbounds if raw

        for m in 1:num_reps

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
        
            for vk in Nc[cimax]
                # 4) update check to node messages Lr[cimax,vnmax]
                li = LinearIndices(Lr)[cimax,vk]
                Lr[li] = newLr[li]
                Nvk = Nv[vk]
                for ci in Nvk
                    if ci ≠ cimax
                        # 5) update Nv messages Lq[ci,vnmax]
                        Lq[ci,vk] = calc_Lq(Nvk,ci,vk,Lr,Lf)   
                        # 6) calculate alpha
                        Nci = Nc[ci]
                        maxresidue = 0.0
                        for vj in Nci
                            newlr = calc_Lr(Nci,ci,vj,Lq)
                            li = LinearIndices(Lr)[ci,vj]
                            newLr[li] = newlr
                            residue = calc_residue_VN_NW_raw(newlr,Lr[li],Factors[ci])
                            if residue > maxresidue
                                maxresidue = residue
                            end
                        end
                        alpha[ci] = maxresidue
                    end
                end
            end
        end

        # 7) update bitvector
        for vk in eachindex(Nv)
            Ld = calc_Ld(vk,Nv[vk],Lf,Lr)
            bitvector[vk] = signbit(Ld)
        end
    else
        for m in 1:num_reps

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

            Ncimax = Nc[cimax]
            for vk in Ncimax
                # 4) update check to node messages Lr[cimax,vnmax]
                li = LinearIndices(Lr)[cimax,vk]
                RBP_update_Lr!(li,Lr,newLr,cimax,vk,Ncimax,Lq,signs,phi)
                # 5) update Nv messages Lq[vk,ci] and bitvector[vk]
                Nvk = Nv[vk]
                Ld = calc_Ld(vk,Nvk,Lf,Lr)
                bitvector[vk] = signbit(Ld)
                # 6) calculate residues
                for ci in Nvk
                    if ci ≠ cimax
                        # 6) update Nv messages Lq[ci,vk]
                        li = LinearIndices(Lq)[ci,vk]
                        Lq[li] = tanhLq(Ld - Lr[li],signs)
                        Nci = Nc[ci]    
                        # calculate the new check to node messages
                        A,B,C,D = calc_ABCD!(Lq,ci,Nci,signs,phi)
                        # calculate alpha
                        maxresidue = 0.0
                        for vj in Nci
                            li = LinearIndices(Lr)[ci,vj]
                            newlr = calc_Lr(A,B,C,D,vj,Lq[li],signs,phi)
                            newLr[li] = newlr
                            residue = abs(newlr - Lr[li])*Factors[ci]
                            if residue > maxresidue
                                maxresidue = residue
                            end
                        end
                        alpha[ci] = maxresidue
                    end
                end
            end
        end
    end

    return bp_not_converged
end

# RAW
function
    NW_RBP!(
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
        num_reps::Int,
        newLr::Matrix{Float64},
        Factors::Vector{Float64},
        alpha::Vector{Float64},
        bp_not_converged::Bool
    )

    @fastmath @inbounds for m in 1:num_reps

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
       
        for vk in Nc[cimax]
            # 4) update check to node messages Lr[cimax,vnmax]
            li = LinearIndices(Lr)[cimax,vk]
            Lr[li] = newLr[li]
            Nvk = Nv[vk]
            for ci in Nvk
                if ci ≠ cimax
                    # 5) update Nv messages Lq[ci,vnmax]
                    Lq[ci,vk] = calc_Lq(Nvk,ci,vk,Lr,Lf)   
                    # 6) calculate alpha
                    Nci = Nc[ci]
                    maxresidue = 0.0
                    for vj in Nci
                        newlr = calc_Lr(Nci,ci,vj,Lq)
                        li = LinearIndices(Lr)[ci,vj]
                        newLr[li] = newlr
                        residue = calc_residue_VN_NW_raw(newlr,Lr[li],Factors[ci])
                        if residue > maxresidue
                            maxresidue = residue
                        end
                    end
                    alpha[ci] = maxresidue
                end
            end
        end
    end

    # 7) update bitvector
    for vk in eachindex(Nv)
        ci = Nv[vk][1]
        li = LinearIndices(Lr)[ci,vk]
        bitvector[vk] = signbit(Lr[li] + Lq[li])
    end

    return bp_not_converged
end

