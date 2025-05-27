################################################################################
# Allan Eduardo Feitosa
# 22 Feb 2024
# Function to find the edge of the maximum residue

# RBP
function
    findmaxedge(
        residues::Vector{<:AbstractFloat},
        ::Nothing     
    )

    max_edge = 0
    maxresidue = 0.0
    @fastmath @inbounds for edge in eachindex(residues)
        residue = abs(residues[edge])
        if residue > maxresidue
            maxresidue = residue
            max_edge = edge
        end
    end

    return max_edge, maxresidue

end

# List-RBP
function
    findmaxedge(
        residues::Vector{<:AbstractFloat},
        ::Vector{<:AbstractFloat}   
    )

    @inbounds return 1, residues[1]

end
