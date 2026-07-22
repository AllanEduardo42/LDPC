################################################################################
# Allan Eduardo Feitosa
# 17 set 2024
# LBP Sum-Product Algorithm

function LBP!(
    bitvector::Vector{Bool},
    V2C::Matrix{Float64},
    C2V::Matrix{Float64},        
    prior_LLRs::Vector{Float64},
    Nc::Vector{Vector{Int}},
    Nv::Vector{Vector{Int}},
    msum_factor::Union{Float64,Nothing}
)

    @fastmath @inbounds for ci in eachindex(Nc)
        # V2C updates 
        Nci = Nc[ci]     
        for vj in Nci # for every vj in Neighborhood(ci)
            V2C[ci,vj] = calc_V2C(Nv[vj],ci,vj,C2V,prior_LLRs)
        end
        # C2V updates
        for vj in Nci
            C2V[ci,vj] = calc_C2V(Nci,ci,vj,V2C,msum_factor)
        end
    end
    for vj in eachindex(Nv)
        post_LLR = calc_post_LLR(vj,Nv[vj],prior_LLRs,C2V)
        bitvector[vj] = signbit(post_LLR)
    end
end

function calc_V2C(
    Nvj::Vector{Int},
    ci::Int,
    vj::Int,
    C2V::Matrix{Float64},
    prior_LLRs::Vector{Float64}
)

    @fastmath @inbounds begin
        lq = prior_LLRs[vj]
        for ca in Nvj
            if ca ≠ ci
                lq += C2V[ca,vj]
            end
        end
        return tanh(0.5*lq)
    end
end