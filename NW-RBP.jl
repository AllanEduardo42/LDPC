################################################################################
# Allan Eduardo Feitosa
# 06 mai 2025
# Nodewise RBP Algorithm with residual decaying factor

include("./RBP functions/findmaxnode.jl")
include("./RBP functions/calc_residue.jl")

function
    NW_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        phi::Union{Vector{Float64},Nothing},
        msum_factor::Union{Float64,Nothing},
        msum2::Bool,
        num_reps::Int,
        newLr::Matrix{Float64},
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
    
        for vk in Nc[cimax]
            # 4) update check to node messages Lr[cimax,vnmax]
            li = LinearIndices(Lr)[cimax,vk]
            Lr[li] = newLr[li]
            Nvk = Nv[vk]
            Ld = calc_Ld(vk,Nvk,Lf,Lr)
            bitvector[vk] = signbit(Ld)
            for ci in Nvk
                alp = 0.0
                if ci ≠ cimax
                    # 5) update Nv messages Lq[ci,vk]
                    li = LinearIndices(Lq)[ci,vk]
                    Lq[li] = tanh(0.5*(Ld - Lr[li]))   
                    # 6) calculate alpha
                    Nci = Nc[ci]
                    for vj in Nci
                        li = LinearIndices(Lr)[ci,vj]
                        alp, _ = calc_residue!(Lq,Lr,newLr,li,ci,vj,Nci,1.0,alp,msum_factor)
                    end
                    alpha[ci] = alp
                end
            end
        end
    end

    return bp_not_converged
end
