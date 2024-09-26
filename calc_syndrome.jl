################################################################################
# Allan Eduardo Feitosa
# 19 sep 2024
# Function to calculate the syndrome in the LDPC-SPA decode algorithm

"""Return the syndrome H*d, where H is the parity check matrix and d an estimate
of the transmited codeword."""
function 
    calc_syndrome!(
        syndrome::Vector{Bool},
        d::Vector{Bool},
        checks2nodes::Vector{Vector{T}} where {T<:Integer}
    )

    syndrome .*= false
    check = 0
    for nodes in checks2nodes
        check += 1
        @inbounds syndrome[check] = _calc_syndrome(d,nodes)
    end
    
end

function 
    calc_syndrome(
        d::Vector{Bool},
        checks2nodes::Vector{Vector{T}} where {T<:Integer}
    )

    syndrome = zeros(Bool,length(checks2nodes))
    check = 0
    for nodes in checks2nodes
        check += 1
        @inbounds syndrome[check] = _calc_syndrome(d,nodes)
    end

    return syndrome

end

function 
    _calc_syndrome(
        d::Vector{Bool},
        nodes::Vector{<:Integer}
    )

    syndrome = false
    for node in nodes
        @inbounds syndrome ⊻= d[node]
    end

    return syndrome
end
    