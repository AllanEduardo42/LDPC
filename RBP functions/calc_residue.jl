################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# Calculate the residue for the RBP algorithm

# FAST
function 
    calc_residue(
        li::Integer,
        Lq::Matrix{<:AbstractFloat},
        newLr::AbstractFloat,
        oldLr::AbstractFloat,
        Factors::Matrix{<:AbstractFloat},
        relative::Bool,
        ::Vector{<:AbstractFloat}
    )


    @fastmath @inbounds begin      
        residue = newLr - oldLr
        return _calc_residue(residue,li,Lq,oldLr,Factors,relative)
    end
end

#TANH
function
    calc_residue(
        li::Integer,
        Lq::Matrix{<:AbstractFloat},
        newLr::AbstractFloat,
        oldLr::AbstractFloat,
        Factors::Matrix{<:AbstractFloat},
        relative::Bool,
        ::Nothing
    )

    @inbounds begin
        residue = newLr - oldLr
        if isnan(residue)
            return 0.0
        else
            return _calc_residue(residue,li,Lq,oldLr,Factors,relative)
        end
    end
end

# core
function _calc_residue(
    residue::AbstractFloat,
    li::Integer,
    Lq::Matrix{<:AbstractFloat},
    oldLr::AbstractFloat,
    Factors::Matrix{<:AbstractFloat},
    relative::Bool
)

    if relative
        old = oldLr + Lq[li]
        residue /= old
    end
    residue *= Factors[li]
    if signbit(residue)
        return -residue
    else
        return residue
    end

end

# NW-RBP and VN-RBP
function calc_residue(
    newLr::AbstractFloat,
    oldLr::AbstractFloat,
    factor::AbstractFloat
)

    @inbounds begin
        residue = newLr - oldLr
        if isnan(residue)
            return 0.0
        else
            @fastmath return abs(residue)*factor
        end
    end
end