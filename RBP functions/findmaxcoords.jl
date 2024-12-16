################################################################################
# Allan Eduardo Feitosa
# 05 dez 2024
# Function to find the coordenates of the maximum residue

# RBP
function
    findmaxcoords!(
        address::Matrix{<:Integer},
        residues::Vector{<:AbstractFloat},
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer}        
    )

    maxedge = 0
    @fastmath @inbounds for e in eachindex(residues)
        residue = residues[e]
        if residue > listres[1]
            listres[1] = residue
            maxedge = e
        end
    end
    if maxedge ≠ 0
        @inbounds listm[1] = address[1,maxedge]
        @inbounds listn[1] = address[2,maxedge]
    end

end

# Local-RBP and List-RBP
function 
    findmaxcoords!(
        ::Nothing,
        ::Nothing,
        ::Vector{<:AbstractFloat},
        ::Vector{<:Integer},
        ::Vector{<:Integer}        
    )

end