################################################################################
# Allan Eduardo Feitosa
# 25 Fev 2026
# UBP-RBP

include("./RBP functions/findmaxedge.jl")

function
    UBP_RBP!(
        bitvector::Vector{Bool},
        V2C::Matrix{Float64},
        C2V::Matrix{Float64},
        prior_LLRs::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        phi::Union{Vector{Float64},Nothing},
        msum_factor::Union{Float64,Nothing},
        msum2::Bool,
        num_reps::Int,
        newC2V::Matrix{Float64},
        Residuals::Matrix{Float64},
        alpha::Vector{Float64},
        UBP::Vector{Bool}
    )

    rbp_not_converged = true
    
    # for e in 1:num_reps
    @fastmath @inbounds for e in 1:num_reps

        # display("e = $e")

        # 1) Find largest residual and coordenates
        cimax, vjmax = findmaxedge_UBP(Residuals,alpha,Nc,UBP)
        if cimax == 0.0
            rbp_not_converged = false
            break # i.e., BP has converged
        end

        UBP[cimax] = false
        
        limax = LinearIndices(C2V)[cimax,vjmax]
        # 2) update check to node message C2V[ci,vjmax]
        if msum2
            C2V[limax] = calc_Lr_no_opt(Nc[cimax],cimax,vjmax,V2C)
        else
            C2V[limax] = newC2V[limax]
        end
        # 3) set maximum residual to zero
        Residuals[limax] = 0.0

        # 4) update LLR[vjmax]
        Nvjmax = Nv[vjmax]
        post_LLR = calc_post_LLR(vjmax,Nvjmax,prior_LLRs,C2V)
        bitvector[vjmax] = signbit(post_LLR)

        for ci in Nvjmax
            if ci ≠ cimax
                UBP[ci] = true
                # 5) update Nv messages V2C[ci,vnmax]
                li = LinearIndices(V2C)[ci,vjmax]
                alp = Residuals[li]
                V2C[li] = tanh_V2C(post_LLR,C2V[li],msum_factor)
                # 6) calculate Residuals
                Nci = Nc[ci]
                for vj in Nci
                    if vj ≠ vjmax
                        li = LinearIndices(C2V)[ci,vj]
                        newc2v = calc_C2V(Nci,ci,vj,V2C,msum_factor)
                        newC2V[li] = newc2v
                        residual = abs(newc2v - C2V[li])
                        if residual > alp
                            alp = residual
                        end
                        Residuals[li] = residual
                    end
                end
                alpha[ci] = alp
            end
        end
    end

    return rbp_not_converged
end

