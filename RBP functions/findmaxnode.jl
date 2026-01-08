################################################################################
# Allan Eduardo Feitosa
# 31 Mar 2025
# Function to find the edge of the maximum residue

function
    findmaxnode(
        alpha::Vector{Float64}        
    )

    cimax = 0
    max_alp = 0.0
    @inbounds @fastmath for e in eachindex(alpha)
        alp = alpha[e]
        if alp > max_alp
            max_alp = alp
            cimax = e
        end
    end

    return cimax

end