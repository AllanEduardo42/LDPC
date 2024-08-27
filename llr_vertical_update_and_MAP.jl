################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Vertical update and MAP estimate of the LLR based Sum-Product Algorithm

function llr_vertical_update_and_MAP(Lq::Matrix{Float64}, Lr::Matrix{Float64},
                                     d::Vector{Int64}, N::Int64, ΔLf::Vector{Float64},
                                     indices_m::Vector{Vector{Int64}})

    for n = 1:N
        Ld = ΔLf[n]
        for m in indices_m[n]
            Ld += Lr[m,n]
        end
        if Ld < 0
            d[n] = 1
        else
            d[n] = 0
        end
        for m in indices_m[n]
            Lq[n,m] = Ld - Lr[m,n]
        end
    end  

    return Lq, d

end

function llr_init_q(Lq::Matrix{Float64},N,ΔLf::Vector{Float64},indices_m)

    for n=1:N
        for m in indices_m[n]
            Lq[n,m] = ΔLf[n]
        end
    end

    return Lq

end