################################################################################
# Allan Eduardo Feitosa
# 22 Feb 2024
# Function to find the edge of the maximum residue

function
    findmaxedge(
        residues::Vector{Float64}  
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

function
    findmaxedge_SVNF(
        residues::Matrix{Float64},
        vj::Integer,
        Nvj::Vector{Int},
        Nc::Vector{Vector{Int}}
    )

    cimax = 0
    vjmax = 0
    maxresidue = 0.0
    @fastmath @inbounds for ca in Nvj
        Nca = Nc[ca]
        for vk in Nca
            if vk ≠ vj
                residue = residues[ca,vk]
                if residue > maxresidue
                    maxresidue = residue
                    cimax = ca
                    vjmax = vk
                end
            end
        end
    end
    
    return cimax, vjmax

end

