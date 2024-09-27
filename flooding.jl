################################################################################
# Allan Eduardo Feitosa
# 27 set 2024
# Flooding Sum-Product Algorithm

include("update_check2nodes_messages.jl")
include("update_node2checks_messages.jl")

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
        update_check2nodes_messages!(
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
        _, d[node] = update_node2checks_messages!(
                        view(Lq,:,node),
                        view(Lr,:,node),
                        Lf[node],
                        checks
                    )
    end
end

### if mode == "MKAY"

function
    flooding!(
        d::Vector{Bool},
        q::Array{<:AbstractFloat,3},
        r::Array{<:AbstractFloat,3},
        f::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lrn::Nothing,
        sn::Nothing,
        phi::Nothing
    )

    # horizontal update
    check = 0
    δq = q[:,:,1]-q[:,:,2]  
    for nodes in checks2nodes
        check += 1
        update_check2nodes_messages!(
            view(r,check,:,:),
            view(δq,check,:),
            nodes
        )
    end

    # vertical update
    node = 0
    Ld = zeros(2)
    for checks in nodes2checks
        node += 1
        @inbounds Ld = f[node,:]
        d[node] = update_node2checks_messages!(
            view(q,:,node,:),
            view(r,:,node,:),
            Ld,
            checks
        )
    end
end