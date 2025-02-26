################################################################################
# Allan Eduardo Feitosa
# 23 Feb 2025
# Calculate the residue for the RBP algorithm

# FAST
function 
    _calc_residue(
        Ms::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        l::Integer,
        ::Vector{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat}
    )
    # @fastmath @inbounds Ld = Lr[l] + Lq'[l]
    # @fastmath @inbounds x = (Ms[l] - Lr[l])/Ld
    @fastmath @inbounds x = Ms[l] - Lr[l]
    @fastmath if signbit(x)
        return -x
    else
        return x
    end
end

#TANH

function 
    _calc_residue(
        Ms::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        l::Integer,
        ::Nothing,
        Lq::Matrix{<:AbstractFloat}
    )

    # @fastmath @inbounds Ld = Lq'[l]
    # @inbounds x = (Ms[l] - Lr[l])/Ld
    @inbounds x = Ms[l] - Lr[l]
    if isnan(x)
        return 0.0
    else
        if signbit(x)
            return -x
        else
            return x
        end        
    end
end