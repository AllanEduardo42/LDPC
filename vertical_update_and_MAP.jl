################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Vertical update and MAP estimate of the Sum-Product Algorithm

function 
    vertical_update_and_MAP!(
        q::Array{Float64,3},
        d::Vector{Bool},
        r::Array{Float64,3},
        f::Matrix{Float64},
        indices_col::Vector{Vector{Int64}}
    )

    d .*= 0
    n = 0
    for indices in indices_col
        n += 1
        @inbounds d0 = f[n,1]
        @inbounds d1 = f[n,2]
        for m in indices
            @inbounds d0 *= r[m,n,1]
            @inbounds d1 *= r[m,n,2]
        end
        if d1 > d0
            @inbounds d[n] = 1
        end
        for m in indices
            @inbounds q0 = d0 / r[m,n,1]
            @inbounds q1 = d1 / r[m,n,2]
            α = q0 + q1
            @inbounds q[n,m,1] = q0/α
            @inbounds q[n,m,2] = q1/α
        end
    end
end

function
    init_q!(
        q::Array{Float64, 3},
        f::Matrix{Float64},
        indices_col::Vector{Vector{Int64}}
    )

    n = 0
    for indices in indices_col
        n += 1
        for m in indices
            @inbounds q[n,m,1] = f[n,1]
            @inbounds q[n,m,2] = f[n,2]
        end
    end
end