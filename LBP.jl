################################################################################
# Allan Eduardo Feitosa
# 17 set 2024
# LBP Sum-Product Algorithm

include("update_Lq.jl")
include("calc_syndrome.jl")

function
    LBP!(
        d::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
        syndrome::Vector{Bool},
        Ldn::Vector{<:AbstractFloat},
        visited_vns::Vector{Bool},
        ilbp::Bool
    )

    visited_vns .*= false
    for m in eachindex(cn2vn)
        # Lq updates       
        for n in cn2vn[m] # every n in Neighborhood(m)
            if @inbounds visited_vns[n]
                @fastmath Lq[n,m] = Ldn[n] - Lr[m,n]
            else
                @inbounds Ldn[n], d[n] = update_Lq!(Lq,Lr,Lf[n],n,vn2cn,Lrn)
                @inbounds visited_vns[n] = true
            end
        end
        # Lr updates
        pLr = 1.0
        for n in cn2vn[m]
            @fastmath @inbounds Lrn[n] = tanh(0.5*Lq[n,m])
            @fastmath @inbounds pLr *= Lrn[n]
        end
        for n in cn2vn[m]
            @fastmath @inbounds x = pLr/Lrn[n]
            if @fastmath abs(x) < 1 # controls divergent values of Lr
                @fastmath @inbounds Ldn[n] -= Lr[m,n]
                @fastmath @inbounds Lr[m,n] = 2*atanh(x)
                @fastmath @inbounds Ldn[n] += Lr[m,n]
                @inbounds d[n] = signbit(Ldn[n])
            end
        end
        if ilbp
            # calc syndrome
            @inbounds syndrome[m] = _calc_syndrome(d,cn2vn[m])
            if iszero(syndrome)
                break
            end
        end
    end
end