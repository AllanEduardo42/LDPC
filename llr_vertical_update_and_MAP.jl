################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Vertical update and MAP estimate of the LLR based Sum-Product Algorithm (with
# "Inf" restriction)

include("MAP.jl")

function 
    llr_vertical_update_and_MAP!(
        Lq::Matrix{<:AbstractFloat},
        d::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        ΔLf::Vector{<:AbstractFloat},
        nodes2checks::Vector{Vector{T}} where {T<:Integer}
    )

    node = 0
    for checks in nodes2checks
        node += 1
        @inbounds d[node], Ld = _MAP!(
            checks,
            ΔLf[node],
            view(Lr,:,node)
        )
        for check in checks
            @inbounds @fastmath Lq[check,node] = Ld - Lr[check,node]
        end
    end
end
