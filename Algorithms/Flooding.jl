################################################################################
# Allan Eduardo Feitosa
# 27 set 2024
# Flooding Sum-Product Algorithm

function flooding!(
    bitvector::Vector{Bool},
    V2C::Matrix{Float64},
    C2V::Matrix{Float64},
    prior_LLRs::Vector{Float64},
    Nc::Vector{Vector{Int}},
    Nv::Vector{Vector{Int}},
    msum_factor::Union{Float64,Nothing}
)

    @inbounds @fastmath begin

        # C2V update
        for ci in eachindex(Nc)
            Nci = Nc[ci]
            for vj in Nci
                C2V[ci,vj] = calc_C2V(Nci,ci,vj,V2C,msum_factor)
            end
        end

        # V2C update
        for vj in eachindex(Nv)
            Nvj = Nv[vj]
            post_LLR = calc_post_LLR(vj,Nvj,prior_LLRs,C2V)
            bitvector[vj] = signbit(post_LLR)
            for ci in Nvj
                li = LinearIndices(V2C)[ci,vj]
                V2C[li] = tanh_V2C(post_LLR,C2V[li],msum_factor)
            end
        end 
    end
end