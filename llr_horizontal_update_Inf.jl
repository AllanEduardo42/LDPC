################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Horizontal update of the LLR based Sum-Product Algorithm using the tanh
# function and no "Inf" restriction

function 
    llr_horizontal_update_Inf!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer}
    )

    check = 0
    for nodes in checks2nodes
        check += 1
        for node in nodes
            @inbounds Lr[check,node] = 1.0
            for n in nodes
                if n ≠ node
                    @inbounds @fastmath Lr[check,node] *= tanh(0.5*Lq[check,n])
                end
            end
            @inbounds @fastmath Lr[check,node] = 2*atanh(Lr[check,node])
        end
    end
end