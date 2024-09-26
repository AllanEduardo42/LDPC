################################################################################
# Allan Eduardo Feitosa
# 26 set 2024
# instantaneos LBP Sum-Product Algorithm

include("check2node_llr.jl")
include("node2check_llr_and_MAP.jl")

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
            d[node] = node2check_llr_and_MAP!(
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
        check2node_llr!(
            view(Lr,check,:),
            view(Lq,check,:),
            nodes,
            Lrn,
            nothing,
            nothing
        )
    end

end