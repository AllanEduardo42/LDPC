################################################################################
# Allan Eduardo Feitosa
# 17 set 2024
# LBP Sum-Product Algorithm

include("check2node_llr.jl")
include("node2check_llr_and_MAP.jl")

function
    LBP!(
        d::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        ΔLf::Vector{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
    )

    check = 0
    for nodes in checks2nodes
        check += 1
        # vertical update        
        for node in nodes
            node2check_llr_and_MAP!(
                view(Lq,:,node),
                view(Lr,:,node),
                ΔLf[node],
                nodes2checks[node];
                MAP = false
            )
        end
        # horizontal update
        check2node_llr!(
            view(Lr,check,:),
            view(Lq,check,:),
            nodes,
            Lrn,
            nothing,
            nothing
        )
    end

    MAP!(
        d,
        nodes2checks,
        ΔLf,
        Lr
    )

end