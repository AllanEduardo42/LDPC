################################################################################
# Allan Eduardo Feitosa
# 17 set 2024
# LBP Sum-Product Algorithm

include("update_check2nodes_messages.jl")
include("update_node2checks_messages.jl")

function
    LBP!(
        d::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
        Ldn::Vector{<:AbstractFloat},
        visited_nodes::Vector{Bool}
    )

    visited_nodes .*= false
    check = 0
    for nodes in checks2nodes
        check += 1
        # vertical update
        ################### original implementation ############################
        # for node in nodes # every node in Neighborhood(check)
        #     Lq[check,node] = Lf[node]
        #     for c in nodes2checks[node]
        #         if c ≠ check
        #             Lq[check,node] += Lr[c,node]
        #         end
        #     end
        # end
        ##################### new implementation ###############################        
        for node in nodes # every node in Neighborhood(check)
            if visited_nodes[node]
                Lq[check,node] = Ldn[node] - Lr[check,node]
            else
                Ldn[node], d[node] = update_node2checks_messages!(
                    view(Lq,:,node),
                    view(Lr,:,node),
                    Lf[node],
                    nodes2checks[node]
                )
                visited_nodes[node] = true
            end
        end
        ################### original implementation ############################
        # update_check2nodes_messages!(
        #     view(Lr,check,:),
        #     view(Lq,check,:),
        #     nodes,
        #     Lrn,
        #     nothing,
        #     nothing
        # )
        ##################### new implementation ###############################
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

    ####################### original implementation ############################
    # MAP!(
    #     d,
    #     nodes2checks,
    #     Lf,
    #     Lr
    # )
end