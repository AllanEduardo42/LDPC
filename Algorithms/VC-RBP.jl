################################################################################
# Allan Eduardo Feitosa
# 01 Dec 2025
# VC-RBP

include("Auxiliary functions/VC_functions.jl")

function VC_RBP!(
    bitvector::Vector{Bool},
    V2C::Matrix{Float64},
    C2V::Matrix{Float64},
    prior_LLRs::Vector{Float64},
    Nc::Vector{Vector{Int}},
    Nv::Vector{Vector{Int}},
    msum_factor::Union{Float64,Nothing},
    msum2::Bool,
    num_reps::Int,
    Residuals::Matrix{Float64},
    alpha::Vector{Float64}
)

    rbp_not_converged = true
    
    @fastmath for e in 1:num_reps

        # display("e = $e")

        # 1) Find largest residual and coordenates
        cimax, vjmax = findmaxedge_VC(Residuals,alpha,Nv)
        if vjmax == 0.0
            rbp_not_converged = false
            break # i.e., BP has converged
        end

        # 2) Set to zero max residual
        Residuals[cimax,vjmax] = 0.0

        Ncimax = Nc[cimax]        
        for vj in Ncimax
            if vj ≠ vjmax
                # 5) update messages C2V[cimax, vj]
                li = LinearIndices(C2V)[cimax,vj]
                alp = Residuals[li]
                C2V[li] = calc_C2V(Ncimax,cimax,vj,V2C,msum_factor)
                # 6) calculate Residuals
                Nvj = Nv[vj]
                post_LLR = calc_post_LLR(vj,Nvj,prior_LLRs,C2V)
                bitvector[vj] = signbit(post_LLR)
                for ci in Nvj
                    if ci ≠ cimax
                        li = LinearIndices(V2C)[ci,vj]
                        newv2c = post_LLR - C2V[li]
                        oldv2c = 2*atanh(V2C[li])
                        residual = abs(newv2c - oldv2c)
                        if residual > alp
                            alp = residual
                        end
                        Residuals[li] = residual
                        V2C[li] = tanh_V2C(newv2c,0.0,msum_factor)
                    end
                end
                alpha[vj] = alp
            end
        end
    end

    return rbp_not_converged
end

