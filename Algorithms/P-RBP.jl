################################################################################
# Allan Eduardo Feitosa
# 13 Mar 2026
# RBP Algorithm

include("Auxiliary functions/RBP_functions.jl")

function parallel_RBP!(
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
    Residuals::Matrix{Float64},
    alpha::Vector{Float64},
    par_par::Int,
    list::Matrix{Int},
    post_LLR::Vector{Float64}
)

    rbp_not_converged = true
    
    # for e in 1:num_reps
    e = 0
    @fastmath @inbounds while e ≤ num_reps

        # display("e = $e")

        for i in 1:par_par

            # 1) Find largest residual and coordenates
            cimax, vjmax = findmaxedge(Residuals,alpha,Nc)
            # display((cimax,vjmax))
            list[i,1] = cimax
            list[i,2] = vjmax
            if cimax == 0
                rbp_not_converged = false
                break # i.e., BP has converged
            end        

            # 2) update C2V message C2V[cimax,vjmax]
            limax = LinearIndices(C2V)[cimax,vjmax]
            if msum2
                C2V[limax] = calc_C2V_no_opt(Nc[cimax],cimax,vjmax,V2C)
            else
                C2V[limax] = newC2V[limax]
            end
            e += 1
            # 3) set maximum residual to zero
            Residuals[limax] = 0.0

            # 4) update LLR[vjmax]
            Nvjmax = Nv[vjmax]
            post_LLR[i] = calc_post_LLR(vjmax,Nvjmax,prior_LLRs,C2V)
            bitvector[vjmax] = signbit(post_LLR[i])
        end
        if rbp_not_converged    
            for i in 1:par_par
                cimax = list[i,1]
                vjmax = list[i,2]
                Nvjmax = Nv[vjmax]
                for ci in Nvjmax
                    if ci ≠ cimax
                        # 5) update V2C messages V2C[ci,vnmax]
                        li = LinearIndices(V2C)[ci,vjmax]
                        V2C[li] = tanh_V2C(post_LLR[i],C2V[li],msum_factor)
                    end
                end
            end

            for i in 1:par_par
                cimax = list[i,1]
                vjmax = list[i,2]
                Nvjmax = Nv[vjmax]
                for ci in Nvjmax
                    if ci ≠ cimax
                        li = LinearIndices(V2C)[ci,vjmax]
                        alp = Residuals[li]
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
        end
    end

    return rbp_not_converged
end

