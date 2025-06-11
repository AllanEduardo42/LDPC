################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# Calculate the residue for the RBP algorithm

# FAST, TABL and MSUM
function 
    calc_residue(
        newlr::Float64,
        oldlr::Float64,
        factor::Float64,
        relative::Bool,
        lq::Float64
    )


    @fastmath begin      
        residue = newlr - oldlr
        return _calc_residue(residue,oldlr,factor,relative,lq)
    end
end

#TANH
function
    calc_residue_raw(
        newlr::Float64,
        oldlr::Float64,
        factor::Float64,
        relative::Bool,
        lq::Float64
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
    residue::Float64,
    oldlr::Float64,
    factor::Float64,
    relative::Bool,
    lq::Float64,
)

    @fastmath begin
        if relative
            rLd = oldlr + lq
            if rLd ≠ 0.0
                residue /= rLd
            else
                residue = INFFLOAT
            end
        end
        return abs(residue)*factor
    end
end

# NW-RBP and VN-RBP (TANH)
function
    calc_residue_VN_NW_raw(
        newlr::Float64,
        oldlr::Float64,
        factor::Float64
    )
    residue = newlr - oldlr
    if isnan(residue)
        return 0.0
    else
        @fastmath return abs(residue)*factor
    end
end