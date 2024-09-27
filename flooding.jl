################################################################################
# Allan Eduardo Feitosa
# 26 set 2024
# Flooding Sum-Product Algorithm

include("check2node_llr.jl")
include("node2check_llr_and_MAP.jl")

function
    flooding!(
        d::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        sn::Union{Vector{Bool},Nothing},
        phi::Union{Vector{<:AbstractFloat},Nothing}
    )

    # horizontal update
    check = 0
    for nodes in checks2nodes
        check += 1
        check2node_llr!(
            view(Lr,check,:),
            view(Lq,check,:),
            nodes,
            Lrn,
            sn,
            phi
        )
    end

    # vertical update
    node = 0
    for checks in nodes2checks
        node += 1
        _, d[node] = node2check_llr_and_MAP!(
                        view(Lq,:,node),
                        view(Lr,:,node),
                        Lf[node],
                        checks
                    )
    end
end