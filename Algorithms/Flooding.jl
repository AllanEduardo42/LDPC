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
    phi::Union{Vector{Float64},Nothing},
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

### if mode == "MKAY"

function flooding!(
    bitvector::Vector{Bool},
    q::Array{Float64,3},
    r::Array{Float64,3},
    f::Matrix{Float64},
    Nc::Vector{Vector{Int}},
    Nv::Vector{Vector{Int}},
    δq::Vector{Float64},
    ::Nothing,
    ::Nothing
)

    @fastmath @inbounds begin

        # horizontal update

        for ci in eachindex(Nc)
            Nci = Nc[ci]
            # S = length(Nci)-1
            for vj in Nci
                δq[vj] = q[ci,vj,1] - q[ci,vj,2]
            end
            for vj in Nci
                δr = calc_δr(Nci,vj,δq)
                r[ci,vj,1] = 0.5*(1+δr)
                r[ci,vj,2] = 0.5*(1-δr)
                # r[ci,vj,1], r[ci,vj,2] = calc_r(q,ci,vj,Nci,S)
            end

        end
        
        # vertical update        
        for vj in eachindex(Nv)
            Nvj = Nv[vj]    
            for ci in Nvj
                Ld1, Ld2 = calc_Ld(r,f,ci,vj,Nvj)
                a = Ld1 + Ld2
                q[ci,vj,1] = Ld1/a
                q[ci,vj,2] = Ld2/a
            end
        end

        for vj in eachindex(Nv) 
            d0 = f[vj,1]
            d1 = f[vj,2]
            for ci in Nv[vj]
                d0 *= r[ci,vj,1]
                d1 *= r[ci,vj,2]
            end
            bitvector[vj] = d1 > d0
        end
    end
end