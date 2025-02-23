################################################################################
# Allan Eduardo Feitosa
# 13 jan 2025
# Calculate the residues for the RBP algorithm

include("../update_Lr.jl")
include("_calc_residue.jl")
include("calc_residues.jl")

function
    find_local_maxresidue!(
        largest_res::Vector{<:AbstractFloat},
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
        largestcoords::Vector{<:Integer}
    )
    
    # calculate the new check to node messages
    update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)

    # calculate the residues
    @fastmath @inbounds for n in cn2vn[m]
        if n ≠ vnmax
            l = LinearIndices(Ms)[m,n]
            x = _calc_residue(Ms,Lr,l,Lrn)
            x *= Factors[l]
            if x > largest_res[1]
                largest_res[1], largest_res[2] = x, largest_res[1]
                largestcoords[3] = largestcoords[1]
                largestcoords[4] = largestcoords[2]
                largestcoords[1] = m
                largestcoords[2] = n
            elseif x > largest_res[2]
                largest_res[2] = x
                largestcoords[3] = m
                largestcoords[4] = n
            end         
        end
    end
end