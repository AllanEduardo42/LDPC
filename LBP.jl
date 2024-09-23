################################################################################
# Allan Eduardo Feitosa
# 17 set 2024
# LBP Sum-Product Algorithm

include("llr_horizontal_update.jl")
include("MAP.jl")

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
            @inbounds Lq[check,node] = ΔLf[node]
            for c in nodes2checks[node]
                if c != check
                    @inbounds @fastmath Lq[check,node] += Lr[c,node]
                end
            end
        end
        # horizontal update
        _llr_horizontal_update_tnh!(
            view(Lr,check,:),
            view(Lq,check,:),
            Lrn,
            nodes
        )
    end

    MAP!(
        d,
        nodes2checks,
        ΔLf,
        Lr
    )

end