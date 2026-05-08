################################################################################
# Allan Eduardo Feitosa
# 13 Mar 2026
# C-RBP, C&R-RBP and C&DR-RBP Algorithms

include("Auxiliary functions/RBP_functions.jl")

function C_and_R_RBP!(
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
    decayfactor::Float64,
    Factors::Matrix{Float64},
    switch_R::Bool
)
    
    rbp_not_converged = true

    # for e in 1:num_reps
    @fastmath @inbounds for e in 1:num_reps

        # display("e = $e")

        # 1) Find largest residual and coordenates
        cimax, vjmax = findmaxedge(Residuals,alpha,Nc)
        if cimax == 0.0
            rbp_not_converged = false
            break # i.e., BP has converged
        end 
        
        ### Consensus

        Nvjmax = Nv[vjmax]

        for ci in Nvjmax
            # 2) update C2V message C2V[cimax,vjmax]
            li = LinearIndices(C2V)[ci,vjmax]

            if msum2
                C2V[li] = calc_C2V_no_opt(Nc[ci],ci,vjmax,V2C)
            else
                C2V[li] = newC2V[li]    # NO EXTRA CALCULATIONS!!!
            end
            # 3) set maximum residual to zero
            Residuals[li] = 0.0

            # 4) Decay the RBP factor corresponding to the maximum residual
            Factors[li] *= decayfactor
        end
        
        # 5) update LLR[vjmax]
        Nvjmax = Nv[vjmax]
        post_LLR = calc_post_LLR(vjmax,Nvjmax,prior_LLRs,C2V)
        bitvector[vjmax] = signbit(post_LLR)

        for ci in Nvjmax
            if switch_R || ci ≠ cimax
                # 6) update V2C messages V2C[ci,vnmax]
                li = LinearIndices(V2C)[ci,vjmax]
                alp = Residuals[li]
                V2C[li] = tanh_V2C(post_LLR,C2V[li],msum_factor)
                # 7) calculate Residuals
                Nci = Nc[ci]
                for vj in Nci
                    if vj ≠ vjmax
                        li = LinearIndices(C2V)[ci,vj]
                        newc2v = calc_C2V(Nci,ci,vj,V2C,msum_factor)
                        newC2V[li] = newc2v
                        residual = abs(newc2v - C2V[li])
                        residual *= Factors[li]
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


