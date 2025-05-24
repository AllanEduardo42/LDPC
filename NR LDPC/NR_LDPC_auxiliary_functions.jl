################################################################################
# Allan Eduardo Feitosa
# 21 Mai 2024
# Auxiliary functions of the NR-LDPC encoding

function 
    get_TBS(
        G::Integer,
        R::AbstractFloat
    )

    TBS = [  24,   32,   40,   48,   56,   64,   72,   80,   88,   96,  104,  112,  
            120,  128,  136,  144,  152,  160,  168,  176,  184,  192,  208,  224,  
            240,  256,  272,  288,  304,  320,  336,  352,  368,  384,  408,  432,  
            456,  480,  504,  528,  552,  576,  608,  640,  672,  704,  736,  768,  
            808,  848,  888,  928,  984, 1032, 1064, 1128, 1160, 1192, 1224, 1256,  
           1288, 1320, 1352, 1416, 1480, 1544, 1608, 1672, 1736, 1800, 1864, 1928,  
           2024, 2088, 2152, 2216, 2280, 2408, 2472, 2536, 2600, 2664, 2728, 2792, 
           2856, 2976, 3104, 3240, 3368, 3496, 3624, 3752, 3824]

    Ninfo = G*R

    if Ninfo ≤ 3824
        n = max(3,floor(Int,log2(Ninfo)) - 6)
        Ninfo_prime = max(24,2^n*(Ninfo ÷ (2^n)))
        _ , idx =  findmin(x -> (x ≥ 0) ? x : Inf, TBS .- Ninfo_prime)
        A = TBS[idx]
    else
        n = floor(Int,log2(Ninfo-24)) - 5
        Ninfo_prime = max(3840,2^n*round(Int,(Ninfo-24)\(2^n)))
        if R ≤ 0.25
            C = cld(Ninfo_prime + 24,3816)
        elseif Ninfo_prime > 8424
            C = cld(Ninfo_prime + 24,3424)
        else
            C = 1
        end
        A = 8*C*cld(Ninfo_prime+24,8*C)-24
    end

    return A

end

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