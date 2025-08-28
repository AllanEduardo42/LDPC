################################################################################
# Allan Eduardo Feitosa
# 17 set 2024
# LBP Sum-Product Algorithm

include("update_Lq.jl")
include("update_Lr.jl")
include("calc_syndrome.jl")

function
    VN_LBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},        
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}}
    )

    @fastmath @inbounds for vj in eachindex(Nv)
        Nvj = Nv[vj]
        # Lr updates
        for ci in Nvj
            Nci = Nc[ci]
            Lr[ci,vj] = calc_Lr(Nci,ci,vj,Lq)
        end
        # Lq updates   
        Ld = calc_Ld(vj,Nv[vj],Lf,Lr)
        bitvector[vj] = signbit(Ld)
        for ci in Nvj # for every vj in Neighborhood(ci)
            li = LinearIndices(Lr)[ci,vj]
            Lq[li] = tanh(0.5*(Ld - Lr[li]))
        end        
    end
end