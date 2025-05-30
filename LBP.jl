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
        Lrj::Vector{Float64},
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
        update_Lr!(Lr,Lq,ci,Nci,Lrj,signs,phi)
        for vj in Nci
            bitvector[vj] = signbit(Lq[ci,vj] + Lr[ci,vj])
        end
    end
end

### old instantaneous LBP
# function
#     LBP!(
#         bitvector::Vector{Bool},
#         Lq::Matrix{Float64},
#         Lr::Matrix{Float64},
#         Lf::Vector{Float64},
#         Ldn::Vector{Float64},
#         Nc::Vector{Vector{T}} where {TInteger},
#         Nv::Vector{Vector{T}} where {TInteger},
#         Lrj::Vector{Float64},
#         syndrome::Vector{Bool},
#         visited_vns::Vector{Bool},
#         i::Integer
#     )

#     visited_vns .= false
#     for ci in eachindex(Nc)
#         # Lq updates       
#         @fastmath @inbounds for vj in Nc[ci] # for every vj in Neighborhood(ci)
#             if visited_vns[vj]
#                 Lq[vj,ci] = Ldn[vj] - Lr[ci,vj]
#             else
#             Ldn[vj], bitvector[vj] = update_Lq!(Lq,Lr,Lf[vj],vj,Nv,Lrj)
#                 visited_vns[vj] = true
#             end
#         end
#         # Lr updates
#         pLr = 1.0
#         @fastmath @inbounds for vj in Nc[ci]
#             Lrj[vj] = tanh(0.5*Lq[vj,ci])
#             pLr *= Lrj[vj]
#         end
#         @fastmath @inbounds for vj in Nc[ci]
#             Ldn[vj] -= Lr[ci,vj]
#             x = pLr/Lrj[vj]
#             if abs(x) < 1 # controls divergent values of Lr                
#                 Lr[ci,vj] = 2*atanh(x)                
#             else
#                 Lr[ci,vj] = x*INFFLOAT
#             end
#             Ldn[vj] += Lr[ci,vj]
#             bitvector[vj] = signbit(Ldn[vj])
#         end
#         # calc syndrome
#         # if i > 1
#         #     @inbounds syndrome[ci] = _calc_syndrome(bitvector,Nc[ci])
#         #     if iszero(syndrome)
#         #         break
#         #     end
#         # end
#     end
# end