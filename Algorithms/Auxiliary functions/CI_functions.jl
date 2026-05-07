################################################################################
# Allan Eduardo Feitosa
# 07 May 2026
# CI-RBP auxiliary functions

function calc_Dn!(
    Dn::Vector{Float64},
    Prob::Vector{Float64},
    newC2V::Matrix{Float64},
    prior_LLRs::Vector{Float64},
    vj::Int,
    Nv::Vector{Vector{Int}}
)

    @inbounds @fastmath begin
        llr = calc_post_LLR(vj,Nv[vj],prior_LLRs,newC2V)
        new_prob = calc_prob(llr)
        Dn[vj] = abs(Prob[vj] - new_prob)
    end

end

function calc_prob(llr::Float64)
    @fastmath begin
        exp_llr = exp(llr)
        return exp_llr/(1 + exp_llr)
    end
end

function find_vjmax(Dn::Vector{Float64})

    max_dn = 0.0
    vjmax = 0
    @inbounds @fastmath for vj in eachindex(Dn)
        dn = Dn[vj]
        if dn > max_dn
            max_dn = dn
            vjmax = vj
        end
    end

    return vjmax

end

function find_cimax(
    Residuals::Matrix{Float64},
    Nv::Vector{Vector{Int}},
    vjmax::Int
)

    maxresidue = 0.0
    cimax = 0
    @inbounds @fastmath for ci in Nv[vjmax]
        residual = Residuals[ci,vjmax]
        if residual > maxresidue
            maxresidue = residual
            cimax = ci
        end
    end

    return cimax

end