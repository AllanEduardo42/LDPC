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
        aux::Vector{Float64},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing}
    )

    @fastmath @inbounds for ci in eachindex(Nc)
        # Lq updates 
        Nci = Nc[ci]     
        for vj in Nci # for every vj in Neighborhood(ci)
            Ld = calc_Ld(vj,Nv[vj],Lf,Lr)
            Lq[ci,vj] = Ld - Lr[ci,vj]
        end
        # Lr updates
        A, B, C, D = calc_ABCD!(aux,signs,phi,Lq,ci,Nci)
        for vj in Nci
            li = LinearIndices(Lr)[ci,vj]
            Lr[li] = calc_Lr(A,B,C,D,vj,aux,signs,phi)
            bitvector[vj] = signbit(Lr[li] + Lq[li])
        end
    end
end