################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Vertical update and MAP estimate of the Sum-Product Algorithm

function 
    vertical_update_and_MAP!(
        q::Array{<:AbstractFloat,3},
        r::Array{<:AbstractFloat,3},
        f::Matrix{<:AbstractFloat},
        nodes2checks::Vector{Vector{T}} where {T<:Integer}
    )

    d = zeros(Bool,N)
    node = 0
    for checks in nodes2checks
        node += 1
        @inbounds d0 = f[node,1]
        @inbounds d1 = f[node,2]
        for check in checks
            @inbounds d0 *= r[check,node,1]
            @inbounds d1 *= r[check,node,2]
        end
        if d1 > d0
            @inbounds d[node] = 1
        end
        for check in checks
            @inbounds q0 = d0 / (r[check,node,1]+eps())
            @inbounds q1 = d1 / (r[check,node,2]+eps())
            α = q0 + q1
            @inbounds q[check,node,1] = q0/α
            @inbounds q[check,node,2] = q1/α
        end
    end

    return d
end

