################################################################################
# Allan Eduardo Feitosa
# 31 Mar 2025
# Function to find the edge of the maximum residue

function
    findmaxnode(
        residues::Vector{Float64}        
    )

    nmax = 0
    maxresidue = 0.0
    @inbounds @fastmath for e in eachindex(residues)
        residue = residues[e]
        if residue > maxresidue
            maxresidue = residue
            nmax = e
        end
    end

    return nmax

end