################################################################################
# Allan Eduardo Feitosa
# 16 Jan 2026
# CI-RBP Algorithm

include("Auxiliary functions/CI_functions.jl")

function CI_RBP!(
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
    Dn::Vector{Float64},
    Prob::Vector{Float64},
    gamma::Float64
)

    rbp_not_converged = true
    
    # for e in 1:num_reps
    @inbounds @fastmath for e in 1:num_reps

        #display("e = $e")

        # 1) Find node vjmax with largest D
        vjmax = find_vjmax(Dn)
        if vjmax == 0
            rbp_not_converged = false
            break # i.e., BP has converged
        end        
            
        if Dn[vjmax] < gamma            
            cimax, vjmax = findmaxedge(Residuals,alpha,Nc)
        else
            cimax = find_cimax(Residuals,Nv,vjmax)
        end
        if cimax == 0.0
            rbp_not_converged = false
            break # i.e., BP has converged
        end

        limax = LinearIndices(C2V)[cimax,vjmax]
        # 2) update check to node message C2V[ci,vjmax]
        if msum2
            C2V[limax] = calc_C2V_no_opt(Nc[cimax],cimax,vjmax,V2C)
        else
            C2V[limax] = newC2V[limax]
        end
        # 3) set maximum residual to zero
        Residuals[limax] = 0.0

        # 4) update LLR[vjmax]
        Nvjmax = Nv[vjmax]
        post_LLR = calc_post_LLR(vjmax,Nvjmax,prior_LLRs,C2V)
        bitvector[vjmax] = signbit(post_LLR)

        Prob[vjmax] = calc_prob(post_LLR)
        Dn[vjmax] = 0.0

        # update alpha[cimax]
        maxalp = 0.0
        for vj in Nc[cimax]
            residual = Residuals[cimax,vj]
            if residual > maxalp
                maxalp = residual
            end
        end
        alpha[cimax] = maxalp

        for ci in Nvjmax
            if ci ≠ cimax
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
                        calc_Dn!(Dn,Prob,newC2V,prior_LLRs,vj,Nv)                    
                    end
                end
                alpha[ci] = alp
            end
        end
    end

    return rbp_not_converged
end


