################################################################################
# Allan Eduardo Feitosa
# 13 jan 2025
# Calculate the residues for the RBP algorithm

function
    find_local_maxresidue!(
        maxresidues::Vector{<:AbstractFloat},
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
    @fastmath @inbounds for n in cn2vn[m]
        if n ≠ vnmax
            l = LinearIndices(Ms)[m,n]
            x = calc_residue(Ms,Lr,l)
            x *= Factors[l]
            if x > maxresidues[1]
                maxresidues[1], maxresidues[2] = x, maxresidues[1]
                maxcoords[3] = maxcoords[1]
                maxcoords[4] = maxcoords[2]
                maxcoords[1] = m
                maxcoords[2] = n
            elseif x > maxresidues[2]
                maxresidues[2] = x
                maxcoords[3] = m
                maxcoords[4] = n
            end         
        end
    end

end