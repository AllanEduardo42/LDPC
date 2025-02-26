################################################################################
# Allan Eduardo Feitosa
# 13 jan 2025
# Calculate the residues for the RBP algorithm

include("../update_Lr.jl")
include("_calc_residue.jl")
include("calc_residues.jl")

function
    find_local_maxresidue!(
        max_residue::Vector{<:AbstractFloat},
        max_coords::Vector{<:Integer},
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        Ms::Matrix{<:AbstractFloat},
        Lr::Union{Matrix{<:AbstractFloat},Nothing},                      
        Lq::Matrix{<:AbstractFloat},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        vnmax::Integer,
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )
    
    # calculate the new check to node messages
    update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)

    # calculate the residues
    @fastmath @inbounds for n in cn2vn[m]
        if n ≠ vnmax
            index = LinearIndices(Ms)[m,n]
            residue = _calc_residue(Ms,Lr,index,Lrn,Lq)
            residue *= Factors[index]
            if residue > max_residue[1]
                max_residue[1], max_residue[2] = residue, max_residue[1]
                max_coords[3] = max_coords[1]
                max_coords[4] = max_coords[2]
                max_coords[1] = m
                max_coords[2] = n
            elseif residue > max_residue[2]
                max_residue[2] = residue
                max_coords[3] = m
                max_coords[4] = n
            end         
        end
    end
end