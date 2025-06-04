################################################################################
# Allan Eduardo Feitosa
# 22 Feb 2024
# Function to find the edge of the maximum residue

function
    findmaxedge(
        residues::Vector{<:AbstractFloat}  
    )

    max_edge = 0
    maxresidue = 0.0
    @fastmath @inbounds for edge in eachindex(residues)
        residue = residues[edge]
        if residue > maxresidue
            maxresidue = residue
            max_edge = edge
        end
    end

    return max_edge

end
