################################################################################
# Allan Eduardo Feitosa
# 13 Mar 2026
# OV-RBP Algorithm

include("Auxiliary functions/init_OV.jl")

function OV_RBP!(
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
    Residuals::Vector{Float64},
    LLRs::Vector{Float64},
    newLLRs::Vector{Float64},
    C::Vector{Bool},
    C1::Vector{Bool},
    upc::Vector{Int},
    max_upc::Int
)

    rbp_not_converged = true

    for e in 1:num_reps

        # display(e)

        if sum(C1) > 0
            max_residue = 0.0
            vjmax = 0
            for vj in eachindex(Nv)
                if C1[vj]
                    residual = Residuals[vj]
                    if residual > max_residue
                        max_residue = residual
                        vjmax = vj
                    end
                end
            end
            if vjmax == 0
                throw(error("vjmax == 0"))
            end
        elseif sum(C) > 0
            max_residue = 0.0
            vjmax = 0
            for vj in eachindex(Nv)
                if C[vj]
                    residual = Residuals[vj]
                    if residual > max_residue
                        max_residue = residual
                        vjmax = vj
                    end
                end
            end
            if vjmax == 0
                throw(error("vjmax == 0"))
            end
            # Find new maximum number of unsatisfied parity check equations
            max_upc = 1
            for p in upc
                if p > max_upc
                    max_upc = p
                end
            end
        else
            max_residue = 0.0
            vjmax = 0
            for vj in eachindex(Nv)
                residual = Residuals[vj]
                if residual > max_residue
                    max_residue = residual
                    vjmax = vj
                end
            end
        end

        if vjmax == 0
            rbp_not_converged = false
            break # i.e., BP has converged
        end

        Nvjmax = Nv[vjmax]

        # 12 - 14
        for ci in Nvjmax
            li = LinearIndices(C2V)[ci,vjmax]
            C2V[li] = newC2V[li]
        end

        # 15
        newLLR_vjmax = newLLRs[vjmax]

        #16 - 19
        if C[vjmax]
            LLRs[vjmax] += newLLR_vjmax
            # LLRs[vjmax] = newLvjmax
            C[vjmax] = false
            C1[vjmax] = false
        else
            LLRs[vjmax] = newLLR_vjmax
        end


        # 20
        Residuals[vjmax] = 0        

        #21 - 23
        post_LLR = calc_post_LLR(vjmax,Nvjmax,prior_LLRs,C2V)
        bitvector[vjmax] = signbit(LLRs[vjmax])
        for ci in Nvjmax
            li = LinearIndices(C2V)[ci,vjmax]
            V2C[li] = tanh_V2C(post_LLR,C2V[li],msum_factor)
            # 3 - 7
            Nci = Nc[ci]
            for vj in Nci
                if vj != vjmax
                    newC2V[ci,vj] = calc_C2V(Nci,ci,vj,V2C,msum_factor) 
                    oldllr = LLRs[vj]
                    newllr = calc_post_LLR(vj,Nv[vj],prior_LLRs,newC2V)
                    newLLRs[vj] = newllr
                    Residuals[vj] = abs(newllr - oldllr)        
                    if sign(oldllr)*sign(newllr) < 0
                        C[vj] = true
                        count_upc = 0
                        for ca in Nv[vj]
                            if _calc_syndrome(bitvector,Nc[ca])
                                count_upc += 1
                            end
                        end
                        upc[vj] = count_upc
                        if count_upc == max_upc  # every VN in C1 is in C
                            C1[vj] = true
                        else
                            C1[vj] = false
                        end
                    else
                        C[vj] = false
                        C1[vj] = false
                    end
                end
            end
        end   
    end

    return rbp_not_converged

end