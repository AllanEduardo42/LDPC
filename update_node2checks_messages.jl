################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Vertical update and MAP estimate of the LLR based Sum-Product Algorithm (with
# "Inf" restriction)

function
    update_node2checks_messages!(
        vq::AbstractMatrix{<:AbstractFloat},
        vr::AbstractMatrix{<:AbstractFloat},
        Ld::Vector{<:AbstractFloat},
        checks::Vector{<:Integer}
    )

    
    for check in checks
        @inbounds @fastmath Ld[1] *= vr[check,1]
        @inbounds @fastmath Ld[2] *= vr[check,2]
    end
    for check in checks
        @inbounds @fastmath q1 = Ld[1] / (vr[check,1] + eps())  
        @inbounds @fastmath q2 = Ld[2] / (vr[check,2] + eps())  
        @fastmath a = q1 + q2
        @inbounds @fastmath vq[check,1] = q1/a
        @inbounds @fastmath vq[check,2] = q2/a
    end

    return signbit(Ld[1]-Ld[2])

end

function
    update_node2checks_messages!(
        vLq::AbstractVector{<:AbstractFloat},
        vLr::AbstractVector{<:AbstractFloat},
        Lf::AbstractFloat,
        checks::Vector{<:Integer}
    )

    Ld = calc_Ld(
        checks,
        Lf,
        vLr
    )
    for check in checks
        @inbounds @fastmath vLq[check] = Ld - vLr[check]
    end

    return Ld, signbit(Ld)
end

function
    update_node2check_message(
        nodes2checks_nmax::Vector{<:Integer},
        check::Integer,
        Lq::AbstractFloat,
        Lr::AbstractVector{<:AbstractFloat}
    )
    
    for c in nodes2checks_nmax
        if c != check
            @inbounds @fastmath Lq += Lr[c]
        end
    end

    return Lq
end

function 
    calc_Ld(
        checks::Vector{<:Integer},
        Lf::AbstractFloat,
        vLr::AbstractVector{<:AbstractFloat}
    )
    Ld = Lf
    for check in checks
        @inbounds @fastmath Ld += vLr[check]
    end
    
    return Ld

end

function 
    MAP!(
        d::Vector{Bool},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lf::Vector{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat}
    )
    
    node = 0
    Ld = 0.0
    for checks in nodes2checks
        node += 1
        Ld = calc_Ld(
                    checks,
                    Lf[node],
                    view(Lr,:,node)
                )
        @inbounds d[node] = signbit(Ld)
    end

end

