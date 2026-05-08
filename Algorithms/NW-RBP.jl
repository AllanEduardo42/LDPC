################################################################################
# Allan Eduardo Feitosa
# 06 mai 2025
# Nodewise RBP Algorithm with residual decaying factor

include("Auxiliary functions/RBP_functions.jl")

function NW_RBP!(
    bitvector::Vector{Bool},
    V2C::Matrix{Float64},
    C2V::Matrix{Float64},
    prior_LLRs::Vector{Float64},
    Nc::Vector{Vector{Int}},
    Nv::Vector{Vector{Int}},
    msum_factor::Union{Float64,Nothing},
    msum2::Bool,
    num_reps::Int,
    newC2V::Matrix{Float64},
    alpha::Vector{Float64}      
)

    rbp_not_converged = true

    @fastmath @inbounds for m in 1:num_reps

        # display("m = $m")

        # 1) Find largest alpha
        cimax = findmaxnode(alpha)
        # display(findmax(alpha))
        if cimax == 0
            rbp_not_converged = false
            break # i.e., BP has converged
        end

        # 2) Set maximum alpha to zero
        alpha[cimax] = 0.0
    
        for vk in Nc[cimax]
            # 4) update C2V messages C2V[cimax,vnmax]
            li = LinearIndices(C2V)[cimax,vk]
            C2V[li] = newC2V[li]
            Nvk = Nv[vk]
            post_LLR = calc_post_LLR(vk,Nvk,prior_LLRs,C2V)
            bitvector[vk] = signbit(post_LLR)
            for ci in Nvk
                alp = 0.0
                if ci ≠ cimax
                    # 5) update V2C messages V2C[ci,vk]
                    li = LinearIndices(V2C)[ci,vk]
                    V2C[li] = tanh(0.5*(post_LLR - C2V[li]))   
                    # 6) calculate alpha
                    Nci = Nc[ci]
                    for vj in Nci
                        li = LinearIndices(C2V)[ci,vj]
                        newc2v = calc_C2V(Nci,ci,vj,V2C,msum_factor)
                        newC2V[li] = newc2v
                        residual = abs(newc2v - C2V[li])
                        if residual > alp
                            alp = residual
                        end
                    end
                    alpha[ci] = alp
                end
            end
        end
    end

    return rbp_not_converged
end
