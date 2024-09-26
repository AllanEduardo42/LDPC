################################################################################
# Allan Eduardo Feitosa
# 26 set 2024
# instantaneos LBP Sum-Product Algorithm

include("llr_horizontal_update.jl")
include("llr_vertical_update_and_MAP.jl")

function
    iLBP!(
        d::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        ΔLf::Vector{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
        syndrome::Vector{Bool}
    )

    check = 0
    for nodes in checks2nodes
        check += 1
        # vertical update        
        for node in nodes
            d[node] = _llr_vertical_update_and_MAP!(
                view(Lq,:,node),
                view(Lr,:,node),
                ΔLf[node],
                nodes2checks[node]
            )
        end
        # calc syndrome
        @inbounds syndrome[check] = _calc_syndrome(d,nodes)
        if iszero(syndrome)
            break
        end            
        # horizontal update
        _llr_horizontal_update!(
            view(Lr,check,:),
            view(Lq,check,:),
            nodes,
            Lrn,
            nothing,
            nothing
        )
    end

    # MAP!(
    #     d,
    #     nodes2checks,
    #     ΔLf,
    #     Lr
    # )

end