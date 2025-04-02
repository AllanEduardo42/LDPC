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
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},        
        Lf::Vector{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat}
    )

    @fastmath @inbounds for m in eachindex(cn2vn)
        # Lq updates       
        for n in cn2vn[m] # for every n in Neighborhood(m)
            _ = update_Lq!(Lq,Lr,Lf[n],n,vn2cn[n],Lrn)
        end
        # Lr updates
        update_Lr!(Lr,Lq,m,cn2vn[m],Lrn,nothing,nothing)
    end

    @fastmath @inbounds for n in eachindex(vn2cn)
        m = vn2cn[n][1]
        bitvector[n] = signbit(Lq[n,m] + Lr[m,n])
    end

end

### old instantaneous LBP
# function
#     LBP!(
#         bitvector::Vector{Bool},
#         Lq::Matrix{<:AbstractFloat},
#         Lr::Matrix{<:AbstractFloat},
#         Lf::Vector{<:AbstractFloat},
#         Ldn::Vector{<:AbstractFloat},
#         cn2vn::Vector{Vector{T}} where {T<:Integer},
#         vn2cn::Vector{Vector{T}} where {T<:Integer},
#         Lrn::Vector{<:AbstractFloat},
#         syndrome::Vector{Bool},
#         visited_vns::Vector{Bool},
#         i::Integer
#     )

#     visited_vns .= false
#     for m in eachindex(cn2vn)
#         # Lq updates       
#         @fastmath @inbounds for n in cn2vn[m] # for every n in Neighborhood(m)
#             if visited_vns[n]
#                 Lq[n,m] = Ldn[n] - Lr[m,n]
#             else
#             Ldn[n], bitvector[n] = update_Lq!(Lq,Lr,Lf[n],n,vn2cn,Lrn)
#                 visited_vns[n] = true
#             end
#         end
#         # Lr updates
#         pLr = 1.0
#         @fastmath @inbounds for n in cn2vn[m]
#             Lrn[n] = tanh(0.5*Lq[n,m])
#             pLr *= Lrn[n]
#         end
#         @fastmath @inbounds for n in cn2vn[m]
#             Ldn[n] -= Lr[m,n]
#             x = pLr/Lrn[n]
#             if abs(x) < 1 # controls divergent values of Lr                
#                 Lr[m,n] = 2*atanh(x)                
#             else
#                 Lr[m,n] = x*INFFLOAT
#             end
#             Ldn[n] += Lr[m,n]
#             bitvector[n] = signbit(Ldn[n])
#         end
#         # calc syndrome
#         # if i > 1
#         #     @inbounds syndrome[m] = _calc_syndrome(bitvector,cn2vn[m])
#         #     if iszero(syndrome)
#         #         break
#         #     end
#         # end
#     end
# end