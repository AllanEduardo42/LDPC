################################################################################
# Allan Eduardo Feitosa
# 27 Jun 2025
# SVNF Algorithm

include("Auxiliary functions/RBP_functions.jl")

function SVNF!( 
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
    twoLS::Int,
    N::Int      
)

    rbp_not_converged = true

    count = 0
    count_zeros = 0

    @inbounds while count < num_reps

        # display("count = $count")
    
        for vn in 0:N-1

            vj = rem(vn + twoLS,N) + 1   # jump to non_punctured nodes in 5GNR

            # display("vj = $vj")

            # 1) Find largest residual and coordenates
            Nvj = Nv[vj]
            cimax, vjmax = findmaxedge_SVNF(Residuals,vj,Nvj,Nc)
            if cimax ≠ 0
                count_zeros = 0

                # 2) update check to node message C2V[cimax,vjmax]
                Ncimax = Nc[cimax]
                limax = LinearIndices(C2V)[cimax,vjmax]
                if msum2
                    C2V[limax] = calc_C2V_no_opt(Ncimax,cimax,vjmax,V2C)
                else
                    C2V[limax] = newC2V[limax]
                end
                
                count += 1
                
                # 3) set maximum residual to zero
                Residuals[limax] = 0.0

                # 4) Calculate post_LLR of vjmax and bitvector[vjmax]
                Nvjmax = Nv[vjmax]
                post_LLR = calc_post_LLR(vjmax,Nvjmax,prior_LLRs,C2V)
                bitvector[vjmax] = signbit(post_LLR)

                for ci in Nvjmax
                    if ci ≠ cimax
                        # 6) update Nv messages V2C[ci,vjmax]
                        li = LinearIndices(V2C)[ci,vjmax]
                        V2C[li] = tanh_V2C(post_LLR,C2V[li],msum_factor)
                        # 7) calculate Residuals
                        Nci = Nc[ci]    
                        for vj in Nci
                            if vj ≠ vjmax
                                li = LinearIndices(C2V)[ci,vj]
                                newc2v = calc_C2V(Nci,ci,vj,V2C,msum_factor)
                                newC2V[li] = newc2v
                                Residuals[li] = abs(newc2v - C2V[li])
                            end
                        end
                    end
                end
            else
                count_zeros += 1        
            end
            if count > num_reps
                break
            end
        end
        if count_zeros ≥ N
            rbp_not_converged = false
            break
        end
    end

    return rbp_not_converged

end