################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Vertical update and MAP estimate of the Sum-Product Algorithm

function 
    vertical_update_and_MAP!(
        q::Array{<:AbstractFloat,3},
        d::Vector{Bool},
        r::Array{<:AbstractFloat,3},
        f::Matrix{<:AbstractFloat},
        nodes2checks::Vector{Vector{T}} where {T<:Integer}
    )

    d .*= 0
    n = 0
    for indices in nodes2checks
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
            @inbounds q0 = d0 / (r[m,n,1]+eps())
            @inbounds q1 = d1 / (r[m,n,2]+eps())
            α = q0 + q1
            @inbounds q[m,n,1] = q0/α
            @inbounds q[m,n,2] = q1/α
        end
    end
end

# function 
#     vertical_update_and_MAP!(
#         q::Array{AbstractFloat,3},
#         d::Vector{Bool},
#         r::Array{AbstractFloat,3},
#         f::Matrix{<:AbstractFloat},
#         nodes2checks::Vector{Vector{T}} where {T<:Integer}
#     )

#     d .*= 0
#     n = 0
#     for indices in nodes2checks
#         n += 1
#         for m in indices
#             @inbounds q[m,n,1] = f[n,1]
#             @inbounds q[m,n,2] = f[n,2]
#             for mm in indices
#                 if mm ≠ m
#                     @inbounds q[m,n,1] *= r[mm,n,1]
#                     @inbounds q[m,n,2] *= r[mm,n,2]
#                 end
#             end
#         end
#         m = indices[1]
#         @inbounds d0 = q[m,n,1] * r[m,n,1]
#         @inbounds d1 = q[m,n,2] * r[m,n,2]
#         if d1 > d0
#             @inbounds d[n] = 1
#         end
#         for m in indices
#             α = q[m,n,1] + q[m,n,2] 
#             @inbounds q[m,n,1] /= α
#             @inbounds q[m,n,2] /= α
#         end
#     end
# end

function
    init_q!(
        q::Array{<:AbstractFloat, 3},
        f::Matrix{<:AbstractFloat},
        nodes2checks::Vector{Vector{T}} where {T<:Integer}
    )

    n = 0
    for indices in nodes2checks
        n += 1
        for m in indices
            @inbounds q[m,n,1] = f[n,1]
            @inbounds q[m,n,2] = f[n,2]
        end
    end
end