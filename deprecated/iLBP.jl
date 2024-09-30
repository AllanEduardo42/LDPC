################################################################################
# Allan Eduardo Feitosa
# 26 set 2024
# instantaneos LBP Sum-Product Algorithm

include("update_Lq.jl")

function
    iLBP!(
        d::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
        syndrome::Vector{Bool},
        Ldn::Vector{<:AbstractFloat},
        visited_vns::Vector{Bool}
    )

    visited_vns .*= false
    for m in eachindex(cn2vn)
        # Lq updates   
        for n in cn2vn[m] # every n in Neighborhood(m)
            if visited_vns[n]
                Lq[n,m] = Ldn[n] - Lr[m,n]
            else
                Ldn[n], _ = update_Lq!(Lq,Lr,Lf[n],n,vn2cn,Lrn)
                visited_vns[n] = true
            end
        end          
        # Lr updates
        pLr = 1.0
        for n in cn2vn[m]
            @inbounds @fastmath Lrn[n] = tanh(0.5*Lq[n,m])
            @inbounds @fastmath pLr *= Lrn[n]
        end
        for n in cn2vn[m]
            @inbounds @fastmath x = pLr/Lrn[n]
            if abs(x) < 1 # controls divergent values of Lr
                Ldn[n] -= Lr[m,n]
                @inbounds @fastmath Lr[m,n] = 2*atanh(x)
                Ldn[n] += Lr[m,n]
                d[n] = signbit(Ldn[n])
            end
        end
        # calc syndrome
        @inbounds syndrome[m] = _calc_syndrome(d,cn2vn[m])
        if iszero(syndrome)
            # println("iLBP: stopped at m=$m ")
            # println()
            break
        end 
    end

end