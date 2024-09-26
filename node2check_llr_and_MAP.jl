################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Vertical update and MAP estimate of the LLR based Sum-Product Algorithm (with
# "Inf" restriction)

function
    node2check_llr_and_MAP!(
        vLq::AbstractVector{<:AbstractFloat},
        vLr::AbstractVector{<:AbstractFloat},
        ΔLf::AbstractFloat,
        checks::Vector{<:Integer};
        MAP = true
    )

    Ld = calc_Ld(
        checks,
        ΔLf,
        vLr
    )
    for check in checks
        @inbounds @fastmath vLq[check] = Ld - vLr[check]
    end

    return MAP ? signbit(Ld) : nothing
end

function 
    calc_Ld(
        checks::Vector{<:Integer},
        ΔLf::AbstractFloat,
        vLr::AbstractVector{<:AbstractFloat}
    )
    Ld = ΔLf
    for check in checks
        @inbounds @fastmath Ld += vLr[check]
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