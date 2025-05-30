################################################################################
# Allan Eduardo Feitosa
# 27 mai 202
# Function to obtain the proper CRC polynomial according to the 3GPP TS 38.212

include("GF2_poly.jl")

function get_CRC_poly(A::Int)::Tuple{Int,Vector{Bool}}
     if A > 3824
        # CRC24A:
        B = A + 24
        p_CRC = "x^24 + x^23 + x^18 + x^17 + x^14 + x^11 + x^10 + x^7 + x^6 + x^5 + x^4 + x^3 + x + 1"
    else
        # CRC16:
        B = A + 16
        p_CRC = "x^16 + x^12 + x^5 + 1"
    end

    return B, gf2_poly(p_CRC)
end