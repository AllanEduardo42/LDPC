function vertical_update_and_MAP!(q::Array{Float64,3},
                                  d::Vector{Int64},
                                  r::Array{Float64,3},
                                  f::Matrix{Float64},
                                  indices_m::Vector{Vector{Int64}})
    d .*= 0
    N = length(indices_m)
    for n = 1:N
        @inbounds d0 = f[n,1]
        @inbounds d1 = f[n,2]
        for m in indices_m[n]
            @inbounds d0 *= r[m,n,1]
            @inbounds d1 *= r[m,n,2]
        end
        if d1 > d0
            @inbounds d[n] = 1
        end
        for m in indices_m[n]
            @inbounds q0 = d0 / r[m,n,1]
            @inbounds q1 = d1 / r[m,n,2]
            α = q0 + q1
            @inbounds q[n,m,1] = q0/α
            @inbounds q[n,m,2] = q1/α
        end
    end
end

function init_q!(q::Array{Float64, 3},
                 N::Int64,
                 f::Matrix{Float64},
                 indices_m::Vector{Vector{Int64}})

    for n=1:N
        for m in indices_m[n]
            @inbounds q[n,m,1] = f[n,1]
            @inbounds q[n,m,2] = f[n,2]
        end
    end
end