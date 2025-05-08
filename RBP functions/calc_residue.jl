################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# Calculate the residue for the RBP algorithm

# FAST
function 
    calc_residue(
        newlr::AbstractFloat,
        oldlr::AbstractFloat,
        factor::AbstractFloat,
        relative::Bool,
        lq::AbstractFloat,
        ::Vector{<:AbstractFloat}
    )


    @fastmath @inbounds begin      
        residue = newlr - oldlr
        return _calc_residue(residue,oldlr,factor,relative,lq)
    end
end

#TANH
function
    calc_residue(
        newlr::AbstractFloat,
        oldlr::AbstractFloat,
        factor::AbstractFloat,
        relative::Bool,
        lq::AbstractFloat,
        ::Nothing
    )

    @inbounds begin
        residue = newlr - oldlr
        if isnan(residue)
            return 0.0
        else
            return _calc_residue(residue,oldlr,factor,relative,lq)
        end
    end
end

# core
function _calc_residue(
    residue::AbstractFloat,
    oldlr::AbstractFloat,
    factor::AbstractFloat,
    relative::Bool,
    lq::AbstractFloat,
)

    @fastmath @inbounds begin
        if relative
            rLd = oldlr + lq
            residue /= rLd
        end
        residue *= factor
        if signbit(residue)
            return -residue
        else
            return residue
        end
    end

end

# NW-RBP and VN-RBP
function calc_residue(
    newlr::AbstractFloat,
    oldlr::AbstractFloat,
    factor::AbstractFloat
)

    @inbounds begin
        residue = newlr - oldlr
        if isnan(residue)
            return 0.0
        else
            @fastmath return abs(residue)*factor
        end
    end
end