################################################################################
# Allan Eduardo Feitosa
# 21 Mai 2024
# Auxiliary functions of the NR-LDPC encoding

function 
    get_Zc(
        K_prime::Integer,
        Kb::Integer
    )

    Z = K_prime/Kb

    liftSizeMtx = [ 2   4   8  16  32  64 128 256;
                    3   6  12  24  48  96 192 384;
                    5  10  20  40  80 160 320   0;
                    7  14  28  56 112 224   0   0;
                    9  18  36  72 144 288   0   0;
                   11  22  44  88 176 352   0   0;
                   13  26  52 104 208   0   0   0;
                   15  30  60 120 240   0   0   0]



    _,coords = findmin(x -> (x ≥ 0) ? x : Inf, liftSizeMtx .- Z)

    Zc = liftSizeMtx[coords]
    iLS = coords[1]-1

    if Kb * Zc < K_prime
        throw(error(
                lazy"""Kb * Zc must be greater than K_prime"."""
            ))
    end

    return Zc, iLS    

end

function 
    get_Er(
        C::Integer,
        G::Integer,
        CBGTI::Vector{<:Any},
        N_L::Integer,
        Q_m::Integer
    )

    CBGTI_flags = get_CBGTI_flags(C, CBGTI)
    C_prime = sum(CBGTI_flags)
    j=0
    E_r = zeros(Int,C)
    @inbounds for r = 1:C
        if CBGTI_flags[r] == 0
            E_r[r] = 0
        else
            if j <= C_prime - rem(G/(N_L*Q_m),C_prime)-1
                E_r[r] = N_L*Q_m*fld(G,N_L*Q_m*C_prime)
            else
                E_r[r] = N_L*Q_m*cld(G,N_L*Q_m*C_prime)
            end
            j += 1
        end
    end
    return E_r
end

function 
    get_CBGTI_flags(
        C::Integer,
        CBGTI::Vector{<:Any}
    )
    
    CBGTI_flags = ones(Bool,C)
    CBGTI_flags[CBGTI[CBGTI .< C] .+ 1] .= false

    return CBGTI_flags
end

function
    get_k0(
        rv::Integer,
        Zc::Integer,
        N_cb::Integer,
        bg::String
    )

    if rv == 0
        k0 = 0
    else
        if bg == "1"
            den = 66       
            if rv == 1
                num = 17
            elseif rv == 2
                num = 33        
            elseif rv == 3
                num = 56        
            else
                throw(ArgumentError(
                    lazy"""rv must be a integer between "0" and "3"."""
                ))
            end
        elseif bg == "2"
            den = 50
            if rv == 1
                num = 13
            elseif rv == 2
                num = 25
            elseif rv == 3
                num = 43
            else
                throw(ArgumentError(
                    lazy"""rv must be a integer between "0" and "3"."""
                ))
            end
        else
            throw(ArgumentError(
                lazy"""bg must be "1" or "2"."""
            ))
        end
        k0 = fld(num*N_cb,den*Zc)*Zc
    end

    return k0
end