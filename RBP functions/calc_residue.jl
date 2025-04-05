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
        li::Integer
    )

    @fastmath @inbounds begin      
        residue = Ms[li] - Lr[li]
        old = Lr[li] + Lq[li]
        residue /= old
        residue *= Factors[li]
        if signbit(residue)
            return -residue
        else
            return residue
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
        residue = Ms[m,n] - Lr[m,n]
        old = Lr[m,n] + Lq[n,m]
        residue /= old
        residue *= Factors[m,n]
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