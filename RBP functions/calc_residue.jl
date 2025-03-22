################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# Calculate the residue for the RBP algorithm

# FAST
function 
    calc_residue(
        Ms::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},
        ::Vector{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        m::Integer,
        n::Integer
    )

    @fastmath @inbounds begin
        li = LinearIndices(Ms)[m,n]        
        residue = Ms[li] - Lr[li]
        old = Lr[li] + Lq'[li]
        residue /= old
        residue *= Factors[li]
        if signbit(residue)
            return -residue, li
        else
            return residue, li
        end
    end
end

#TANH
function 
    calc_residue(
        Ms::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},
        ::Nothing,
        Lq::Matrix{<:AbstractFloat},
        m::Integer,
        n::Integer
    )

    @fastmath @inbounds begin
        li = LinearIndices(Ms)[m,n]
        residue = Ms[li] - Lr[li]
        # old = Lr[li] + Lq'[li]
        # residue /= old
        residue *= Factors[li]
        if isnan(residue)
            return 0.0
        else
            if signbit(residue)
                return -residue, li
            else
                return residue, li
            end
        end
    end
end