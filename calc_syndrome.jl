################################################################################
# Allan Eduardo Feitosa
# 19 sep 2024
# Function to calculate the syndrome in the LDPC-SPA decode algorithm

"""Return the syndrome H*d, where H is the parity n matrix and d an estimate
of the transmited codeword."""
function 
    calc_syndrome!(
        syndrome::Vector{Bool},
        d::Vector{Bool},
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )

    syndrome .*= false
    for n in eachindex(cn2vn)
        @inbounds syndrome[n] = _calc_syndrome(d,cn2vn[n])
    end
    
end

function 
    _calc_syndrome(
        d::Vector{Bool},
        varnodes_cn::Vector{<:Integer}
    )

    syndrome = false
    for n in varnodes_cn
        @inbounds syndrome ⊻= d[n]
    end

    return syndrome
end
    