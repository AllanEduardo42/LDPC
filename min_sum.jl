################################################################################
# Allan Eduardo Feitosa
# 22 set 2024
# Horizontal update of the LLR based MIN SUM Algorithm

function abs_sign!(vLq::AbstractFloat,s::Bool)
    sn = signbit(vLq)
    return abs(vLq), sn, s ⊻ sn
end

function _min_sum!(                           
    vLq::AbstractVector{<:AbstractFloat},
    sn::Vector{Bool},
    nodes::Vector{<:Integer},        
    )

    s = false
    minL = Inf
    minL2 = Inf
    max_idx = 0
    for node in nodes
        @inbounds β, sn[node], s = abs_sign!(vLq[node],s)
        if β < minL
            max_idx = node
            minL, minL2 = β, minL
        elseif β < minL2
            minL2 = β
        end
    end

    return minL, minL2, s, max_idx

end

function 
    __min_sum!(
        node::Integer,
        sn::Bool,
        minL::AbstractFloat,
        minL2::AbstractFloat,
        s::Bool,
        max_idx::Integer
    )

    if node == max_idx #(pick the second least vLq)
        return (1 - 2*(sn ⊻ s))*minL2
    else
        return (1 - 2*(sn ⊻ s))*minL
    end

end

