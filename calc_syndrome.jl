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
        Nc::Vector{Vector{Int}}
    )

    syndrome .= false
    @inbounds for ci in eachindex(Nc)
        syndrome[ci] = _calc_syndrome(bitvector,Nc[ci])
    end
    
end

function 
    _calc_syndrome(
        bitvector::Vector{Bool},
        Nci::Vector{Int}
    )

    syndrome = false
    @inbounds for vj in Nci
        syndrome ⊻= bitvector[vj]
    end

    return syndrome
end
    