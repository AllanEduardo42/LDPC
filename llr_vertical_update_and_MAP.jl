################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Vertical update and MAP estimate of the LLR based Sum-Product Algorithm

function llr_vertical_update_and_MAP!(Lq::Matrix{Float64},
                                      d::Vector{Bool},
                                      Lr::Matrix{Float64},
                                      ΔLf::Vector{Float64},
                                      indices_m::Vector{Vector{Int64}})

    d .*= false
    n = 0
    for indices in indices_m
        n += 1
        @inbounds Ld = ΔLf[n]
        for m in indices
            @fastmath @inbounds Ld += Lr[m,n]
        end
        if Ld < 0
            d[n] = 1
        end
        for m in indices
            @fastmath @inbounds Lq[n,m] = Ld - Lr[m,n]
        end
    end
end


function llr_init_q!(Lq::Matrix{Float64},
                     ΔLf::Vector{Float64},
                     indices_m::Vector{Vector{Int64}})

    n = 0
    for indices in indices_m
        n += 1
        for m in indices
            @inbounds Lq[n,m] = ΔLf[n]
        end
    end
end