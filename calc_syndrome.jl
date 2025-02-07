################################################################################
# Allan Eduardo Feitosa
# 19 sep 2024
# Function to calculate the syndrome in the LDPC-SPA decode algorithm

"""Return the syndrome H*bitvector, where H is the parity n matrix and bitvector an estimate
of the transmited codeword."""
function 
    calc_syndrome!(
        syndrome::Vector{Bool},
        bitvector::Vector{Bool},
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )

    syndrome .*= false
    for m in eachindex(cn2vn)
        @inbounds syndrome[m] = _calc_syndrome(bitvector,cn2vn[m])
    end
    
end

function 
    _calc_syndrome(
        bitvector::Vector{Bool},
        varnodes_cn::Vector{<:Integer}
    )

    syndrome = false
    for n in varnodes_cn
        @inbounds syndrome ⊻= bitvector[n]
    end

    return syndrome
end
    