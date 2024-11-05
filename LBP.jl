################################################################################
# Allan Eduardo Feitosa
# 17 set 2024
# LBP Sum-Product Algorithm

include("update_Lq.jl")
include("calc_syndrome.jl")

function
    LBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
        syndrome::Vector{Bool},
        Ldn::Vector{<:AbstractFloat},
        visited_vns::Vector{Bool},
        ilbp::Bool
    )

    visited_vns .= false
    for m in eachindex(cn2vn)
        # Lq updates       
        @fastmath @inbounds for n in cn2vn[m] # for every n in Neighborhood(m)
            if visited_vns[n]
                Lq[n,m] = Ldn[n] - Lr[m,n]
            else
            Ldn[n], bitvector[n] = update_Lq!(Lq,Lr,Lf[n],n,vn2cn,Lrn)
                visited_vns[n] = true
            end
        end
        # Lr updates
        pLr = 1.0
        @fastmath @inbounds for n in cn2vn[m]
            Lrn[n] = tanh(0.5*Lq[n,m])
            pLr *= Lrn[n]
        end
        @fastmath @inbounds for n in cn2vn[m]
            Ldn[n] -= Lr[m,n]
            x = pLr/Lrn[n]
            if abs(x) < 1 # controls divergent values of Lr                
                Lr[m,n] = 2*atanh(x)                
            else
                Lr[m,n] = x*INFFLOAT
            end
            Ldn[n] += Lr[m,n]
            bitvector[n] = signbit(Ldn[n])
        end
        if ilbp
            # calc syndrome
            @inbounds syndrome[m] = _calc_syndrome(bitvector,cn2vn[m])
            if iszero(syndrome)
                break
            end
        end
    end
end