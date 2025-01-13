################################################################################
# Allan Eduardo Feitosa
# 13 jan 2025
# Calculate the residues for the RBP algorithm

include("findmaxresidue.jl")

function
    find_local_maxresidue!(
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        Ms::Matrix{<:AbstractFloat},
        Lr::Union{Matrix{<:AbstractFloat},Nothing},                      
        Lq::Matrix{<:AbstractFloat},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        vnmax::Integer,
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        maxcoords::Vector{<:Integer}
    )
    
    # calculate the new check to node messages
    update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)

    # calculate the residues
    maxresidue = 0
    @fastmath @inbounds for n in cn2vn[m]
        if n ≠ vnmax
            l = LinearIndices(Ms)[m,n]
            x = calc_residue(Ms,Factors,Lr,l)
            if x != 0.0
                if x > maxresidue
                    maxresidue = x
                    maxcoords[1] = m
                    maxcoords[2] = n
                end
            end            
        end
    end

    if maxresidue == 0 # no update: this will triger random selection of a check
        maxresidue = -1
    end

    return maxresidue

end