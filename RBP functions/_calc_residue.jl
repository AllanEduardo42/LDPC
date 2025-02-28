################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# Calculate the residue for the RBP algorithm

# FAST
function 
    _calc_residue(
        Ms::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},
        ::Vector{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        m::Integer,
        n::Integer
    )

    @fastmath @inbounds begin
        index = LinearIndices(Ms)[m,n]
        # Ld = Lr[index] + Lq'[index]
        residue = Ms[index] - Lr[index]
        # residue /= Ld
        residue *= Factors[index]
        if signbit(residue)
            return -residue, index
        else
            return residue, index
        end
    end
end

#TANH
function 
    _calc_residue(
        Ms::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},
        ::Nothing,
        Lq::Matrix{<:AbstractFloat},
        m::Integer,
        n::Integer
    )

    @fastmath @inbounds begin
        index = LinearIndices(Ms)[m,n]
        # Ld = Lr[index] + Lq'[index]
        residue = Ms[index] - Lr[index]
        # residue /= Ld
        residue *= Factors[index]
        if isnan(residue)
            return 0.0
        else
            if signbit(residue)
                return -residue, index
            else
                return residue, index
            end
        end
    end
end