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
    
    node = 0
    for checks in nodes2checks
        node += 1
        @inbounds d[node], _ = _MAP!(
            checks,
            ΔLf[node],
            view(Lr,:,node)
        )
    end

end

function 
    _MAP!(
        checks::Vector{<:Integer},
        ΔLf::AbstractFloat,
        Lr::AbstractVector{<:AbstractFloat}
    )
    Ld = ΔLf
    for check in checks
        @inbounds @fastmath Ld += Lr[check]
    end
    
    return signbit(Ld), Ld

end