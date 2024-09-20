################################################################################
# Allan Eduardo Feitosa
# 19 sep 2024
# Flooding LLR-SPA algorithm

include("llr_horizontal_update.jl")
include("llr_vertical_update_and_MAP.jl")

function 
    llr_flooding!(
        d::Vector{Bool},
        Lr::Matrix{<:AbstractFloat}, 
        Lq::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        ΔLf::Vector{<:AbstractFloat},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        sn::Union{Vector{Int8},Nothing},
        phi::Union{Vector{<:AbstractFloat},Nothing},
    )

    llr_horizontal_update!(
        Lr,
        Lq,
        checks2nodes,
        Lrn,
        sn,
        phi
    )
    llr_vertical_update_and_MAP!(
        Lq,
        d,
        Lr,
        ΔLf,
        nodes2checks
    )
end