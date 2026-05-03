################################################################################
# Allan Eduardo Feitosa
# 10 Mar 2026
# RBP with D1VN processing Scheme
# Tofar C.-Y. Chang et al. “Enhanced Informed Dynamic BP Decoding Scheduling Strategies for 5G NR LDPC Codes”. In: 2022 IEEE 96th Vehicular Technology Conference (VTC2022-Fall). 2022, pp. 1–6.


include("./RBP functions/findmaxedge.jl")

function
    RBP_D1VN!(
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
        F::Vector{Bool},
        degree_vn::Vector{Int}
    )

    rbp_not_converged = true

    # for e in 1:num_reps + 1
    @fastmath @inbounds for e in 1:num_reps + 1

        # display("e = $e") 
        
        if e ≤ num_reps

            # 1) Find largest residual and coordenates
            cimax, vjmax = findmaxedge_D1VN(Residuals,alpha,Nc,F)
            if cimax == 0.0
                rbp_not_converged = false
                break # i.e., BP has converged
            end

            if degree_vn[vjmax] == 1
                F[vjmax] = true
                Ncimax = Nc[cimax]
                maxresidue = 0.0
                vjmax2 = Ncimax[1]
                for vj in Ncimax                    
                    if vj != vjmax
                        residual = Residuals[cimax,vj]
                        if residual > maxresidue
                            vjmax2 = vj
                            maxresidue = residual
                        end
                    end
                end
                vjmax = vjmax2
            end

            # 2) C2V update
            limax = LinearIndices(V2C)[cimax,vjmax]
            C2V[limax] = newC2V[limax]

            # 3) set maximum residual to zero
            Residuals[limax] = 0.0 

            # 4) Update LLR and estimate bit message
            Nvjmax = Nv[vjmax]
            post_LLR = calc_post_LLR(vjmax,Nvjmax,prior_LLRs,C2V)
            bitvector[vjmax] = signbit(post_LLR)

            for ci in Nvjmax
                if ci ≠ cimax
                    # 5) update V2C messages
                    li = LinearIndices(V2C)[ci,vjmax]
                    # alp = Residuals[li]
                    V2C[li] = tanh_V2C(post_LLR,C2V[li],msum_factor)
                    # 6) calculate Residuals
                    Nci = Nc[ci]
                    for vj in Nci
                        if vj ≠ vjmax && !F[vj]
                            li = LinearIndices(C2V)[ci,vj]
                            newc2v = calc_C2V(Nci,ci,vj,V2C,msum_factor)
                            newC2V[li] = newc2v
                            residual = abs(newc2v - C2V[li])
                            # if residual > alp
                            #     alp = residual
                            # end
                            Residuals[li] = residual
                        end
                    end
                    # alpha[ci] = alp
                end
            end
        else
            F .= false
            for vjmax in eachindex(degree_vn)
                if degree_vn[vjmax] == 1
                    cimax = Nv[vjmax][1]
                    # 2) C2V update
                    li = LinearIndices(V2C)[cimax,vjmax]
                    C2V[li] = newC2V[li]

                    # 3) set maximum residual to zero
                    Residuals[li] = 0.0 

                    # 4) Update LLR and estimate bit message
                    Nvjmax = Nv[vjmax]
                    post_LLR = calc_post_LLR(vjmax,Nvjmax,prior_LLRs,C2V)
                    bitvector[vjmax] = signbit(post_LLR)            
                end
            end
        end
    end

    return rbp_not_converged
end


