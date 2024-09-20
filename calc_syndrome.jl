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
    m = 0
    for indices in checks2nodes
        m += 1
        for n in indices
            @inbounds syndrome[m] ⊻= d[n]
        end
    end

end