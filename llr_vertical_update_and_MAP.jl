################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Vertical update and MAP estimate of the LLR based Sum-Product Algorithm (with
# "Inf" restriction)

function
    _llr_vertical_update_and_MAP!(
        v_Lq::AbstractVector{<:AbstractFloat},
        v_Lr::AbstractVector{<:AbstractFloat},
        ΔLf::AbstractFloat,
        checks::Vector{<:Integer};
        MAP = true
    )

    Ld = calc_Ld(
        checks,
        ΔLf,
        v_Lr
    )
    for check in checks
        @inbounds @fastmath v_Lq[check] = Ld - v_Lr[check]
    end

    return MAP ? signbit(Ld) : nothing
end

function 
    calc_Ld(
        checks::Vector{<:Integer},
        ΔLf::AbstractFloat,
        Lr::AbstractVector{<:AbstractFloat}
    )
    Ld = ΔLf
    for check in checks
        @inbounds @fastmath Ld += Lr[check]
    end
    
    return Ld

end

function 
    MAP!(
        d::Vector{Bool},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        ΔLf::Vector{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat}
    )
    
    node = 0
    Ld = 0.0
    for checks in nodes2checks
        node += 1
        Ld = calc_Ld(
                    checks,
                    ΔLf[node],
                    view(Lr,:,node)
                )
        @inbounds d[node] = signbit(Ld)
    end

end