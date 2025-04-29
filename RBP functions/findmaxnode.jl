################################################################################
# Allan Eduardo Feitosa
# 31 Mar 2025
# Function to find the edge of the maximum residue

function
    findmaxnode(
        residues::Vector{<:AbstractFloat}        
    )

    nmax = 0
    maxresidue = 0.0
    @fastmath @inbounds for e in eachindex(residues)
        residue = residues[e]
        if residue > maxresidue
            maxresidue = residue
            nmax = e
        end
    end

    return nmax

end