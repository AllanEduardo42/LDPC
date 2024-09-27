################################################################################
# Allan Eduardo Feitosa
# 19 sep 2024
# Function to calculate the syndrome in the LDPC-SPA decode algorithm

"""Return the syndrome H*d, where H is the parity cn matrix and d an estimate
of the transmited codeword."""
function 
    calc_syndrome!(
        syndrome::Vector{Bool},
        d::Vector{Bool},
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )

    syndrome .*= false
    for cn in eachindex(cn2vn)
        @inbounds syndrome[cn] = _calc_syndrome(d,cn2vn[cn])
    end
    
end

function 
    _calc_syndrome(
        d::Vector{Bool},
        varnodes_cn::Vector{<:Integer}
    )

    syndrome = false
    for cn in varnodes_cn
        @inbounds syndrome ⊻= d[cn]
    end

    return syndrome
end
    