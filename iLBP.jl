################################################################################
# Allan Eduardo Feitosa
# 26 set 2024
# instantaneos LBP Sum-Product Algorithm

include("update_check2nodes_messages.jl")
include("update_node2checks_messages.jl")

function
    iLBP!(
        d::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
        syndrome::Vector{Bool},
        Ldn::Vector{<:AbstractFloat},
        visited_nodes::Vector{Bool}
    )

    visited_nodes .*= false
    check = 0
    for nodes in checks2nodes
        check += 1
        # vertical update   
        for node in nodes # every node in Neighborhood(check)
            if visited_nodes[node]
                Lq[check,node] = Ldn[node] - Lr[check,node]
            else
                Ldn[node], d[node] = update_node2checks_messages!(
                    Lq,
                    Lr,
                    Lf[node],
                    node,
                    nodes2checks[node]
                )
                visited_nodes[node] = true
            end
        end
        # calc syndrome
        @inbounds syndrome[check] = _calc_syndrome(d,nodes)
        if iszero(syndrome)
            break
        end            
        # horizontal update
        pLr = 1.0
        for node in nodes
            @inbounds @fastmath Lrn[node] = tanh(0.5*Lq[check,node])
            @inbounds @fastmath pLr *= Lrn[node]
        end
        for node in nodes
            @inbounds @fastmath x = pLr/Lrn[node]
            if abs(x) < 1 # controls divergent values of Lr
                Ldn[node] -= Lr[check,node]
                @inbounds @fastmath Lr[check,node] = 2*atanh(x)
                Ldn[node] += Lr[check,node]
                d[node] = signbit(Ldn[node])
            end
        end
    end

end