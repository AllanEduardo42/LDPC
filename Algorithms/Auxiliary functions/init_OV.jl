################################################################################
# Allan Eduardo Feitosa
# 22 Feb 2024
# auxiliary functions for OV-RBP

function init_OV!(
    V2C::Matrix{Float64},
    Nc::Vector{Vector{Int}},
    Nv::Vector{Vector{Int}},
    newC2V::Matrix{Float64}, 
    Residuals::Vector{Float64},
    msum_factor::Union{Float64,Nothing},
    prior_LLRs::Vector{Float64},
    LLRs::Vector{Float64},
    newLLRs::Vector{Float64},
    C::Vector{Bool},
    C1::Vector{Bool},
    upc::Vector{Int}
)

    @inbounds @fastmath begin
        LLRs .= copy(prior_LLRs)
        # 3 - 7
        for ci in eachindex(Nc)
            Nci = Nc[ci]
            for vj in Nci
                newC2V[ci,vj] = calc_C2V(Nci,ci,vj,V2C,msum_factor) 
            end
        end
        #8 - 10                    
        C .= false                      # set C of VNs whose sign of the LLR changes
        upc .= 0                        # unsatisfied parity check equations
        max_upc = 0                     # maximum number of unsatified parity check equations 
        for vj in eachindex(Nv)
            oldllr = LLRs[vj]
            newllr = calc_post_LLR(vj,Nv[vj],prior_LLRs,newC2V)
            newLLRs[vj] = newllr
            Residuals[vj] = abs(newllr - oldllr)        
            if sign(oldllr)*sign(newllr) < 0
                C[vj] = true
                count_upc = 0
                for ci in Nv[vj]
                    if _calc_syndrome(bitvector,Nc[ci])
                        count_upc += 1
                    end
                end
                upc[vj] = count_upc
                if count_upc > max_upc
                    max_upc = count_upc
                end
            end
        end

        # 11
        C1 .= false                     # set of VNs with upc = max_upc
        for vj in eachindex(Nv)
            p = upc[vj]
            if C[vj] && (p == max_upc)  # every VN in C1 is in C
                C1[vj] = true
            end
        end 
    end

    return max_upc
end