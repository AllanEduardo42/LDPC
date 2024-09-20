################################################################################
# Allan Eduardo Feitosa
# 19 sep 2024
# Flooding SPA algorithm

include("simple_horizontal_update.jl")
include("vertical_update_and_MAP.jl")

function
    flooding!(
        d::Vector{Bool},
        r::Array{<:AbstractFloat,3}, 
        q::Array{<:AbstractFloat,3},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        f::Matrix{<:AbstractFloat},
    )
     ### Conventional simplified SPA

    δQ = q[:,:,1]-q[:,:,2]    

    simple_horizontal_update!(
        r,
        δQ,
        checks2nodes
    )
    vertical_update_and_MAP!(
        q,
        d,
        r,
        f,
        nodes2checks
    )
end