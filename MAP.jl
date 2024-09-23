################################################################################
# Allan Eduardo Feitosa
# 19 sep 2024
# Maximum A Posteriori estimate for LDPC code based on the SPA algorithm

function 
    MAP!(
        d::Vector{Bool},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        ΔLf::Vector{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat}
    )
    
    d .*= false
    node = 0
    for checks in nodes2checks
        node += 1
        @inbounds Ld = ΔLf[node]
        for check in checks
            @inbounds @fastmath Ld += Lr[check,node]
        end
        if Ld < 0
            @inbounds d[node] = 1
        end
    end
end