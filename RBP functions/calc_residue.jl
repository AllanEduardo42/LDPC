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
        li::Integer,
        relative::Bool
    )

    @fastmath @inbounds begin      
        residue = Ms[li] - Lr[li]
        if relative
            old = Lr[li] + Lq[li]
            residue /= old
        end
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
        li::Integer,
        relative::Bool
    )

    @inbounds begin
        residue = Ms[li] - Lr[li]
        if isnan(residue)
            return 0.0
        else
            @fastmath begin 
                if relative
                old = Lr[li] + Lq[li]
                residue /= old
                end
                residue *= Factors[li]
                if signbit(residue)
                    return -residue
                else
                    return residue
                end
            end
        end
    end
end