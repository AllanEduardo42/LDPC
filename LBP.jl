################################################################################
# Allan Eduardo Feitosa
# 17 set 2024
# LBP Sum-Product Algorithm

include("update_Lq.jl")
include("update_Lr.jl")
include("calc_syndrome.jl")

function
    LBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},        
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing}
    )

    @fastmath @inbounds for ci in eachindex(Nc)
        # Lq updates 
        Nci = Nc[ci]     
        for vj in Nci # for every vj in Neighborhood(ci)
            Lq[ci,vj] = calc_Lq(Nv[vj],ci,vj,Lr,Lf)
        end
        # Lr updates
        for vj in Nci
            Lr[ci,vj] = calc_Lr(Nci,ci,vj,Lq)
        end
    end
    for vj in eachindex(Nv)
        Ld = calc_Ld(vj,Nv[vj],Lf,Lr)
        bitvector[vj] = signbit(Ld)
    end
end