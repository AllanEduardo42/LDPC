################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# Calculate the residue for the RBP algorithm

# FAST, TABL and MSUM
function 
    calc_residue(
        newlr::AbstractFloat,
        oldlr::AbstractFloat,
        factor::AbstractFloat,
        relative::Bool,
        lq::AbstractFloat
    )


    @fastmath begin      
        residue = newlr - oldlr
        return _calc_residue(residue,oldlr,factor,relative,lq)
    end
end

#TANH
function
    calc_residue_raw(
        newlr::AbstractFloat,
        oldlr::AbstractFloat,
        factor::AbstractFloat,
        relative::Bool,
        lq::AbstractFloat
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

    @fastmath begin
        if relative
            rLd = oldlr + lq
            residue /= rLd
        end
        return abs(residue)*factor
    end
end

# NW-RBP and VN-RBP (FAST, TABL and MSUM)
function 
    calc_residue(
        newlr::AbstractFloat,
        oldlr::AbstractFloat,
        factor::AbstractFloat,
    )

    @fastmath abs(newlr - oldlr)*factor

end

# NW-RBP and VN-RBP (TANH)
function
    calc_residue_raw(
        newlr::AbstractFloat,
        oldlr::AbstractFloat,
        factor::AbstractFloat
    )
    residue = newlr - oldlr
    if isnan(residue)
        return 0.0
    else
        @fastmath return abs(residue)*factor
    end
end