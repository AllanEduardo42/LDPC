################################################################################
# Allan Eduardo Feitosa
# 17 set 2024
# LBP Sum-Product Algorithm

function LBP!(
    bitvector::Vector{Bool},
    Lq::Matrix{Float64},
    Lr::Matrix{Float64},        
    Lf::Vector{Float64},
    Nc::Vector{Vector{Int}},
    Nv::Vector{Vector{Int}},
    phi::Union{Vector{Float64},Nothing},
    msum_factor::Union{Float64,Nothing}
)

    @fastmath @inbounds for ci in eachindex(Nc)
        # Lq updates 
        Nci = Nc[ci]     
        for vj in Nci # for every vj in Neighborhood(ci)
            Lq[ci,vj] = calc_V2C(Nv[vj],ci,vj,Lr,Lf)
        end
        # Lr updates
        for vj in Nci
            Lr[ci,vj] = calc_C2V(Nci,ci,vj,Lq,msum_factor)
        end
    end
    for vj in eachindex(Nv)
        Ld = calc_post_LLR(vj,Nv[vj],Lf,Lr)
        bitvector[vj] = signbit(Ld)
    end
end

function calc_V2C(
    Nvj::Vector{Int},
    ci::Int,
    vj::Int,
    C2V::Matrix{Float64},
    Lf::Vector{Float64}
)

    @fastmath @inbounds begin
        lq = Lf[vj]
        for ca in Nvj
            if ca ≠ ci
                lq += C2V[ca,vj]
            end
        end
        return tanh(0.5*lq)
    end
end