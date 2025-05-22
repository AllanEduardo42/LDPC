################################################################################
# Allan Eduardo Feitosa
# 13 NOV 2024
# Function to encode a message using the NR-LDPC

using SparseArrays

include("../GF2_poly.jl")
include("../GF2_functions.jl")
include("NR_LDPC_functions.jl")

const R_LBRM = 2//3

struct nr_ldpc_data
    B::Integer
    K_prime::Integer
    K::Integer
    Zc::Integer
    bg::String
    N_cb::Integer
    E_r::Vector{<:Integer}
    k0::Integer
    g_CRC::Vector{Bool}
    P_Zc::Integer
end


function
    NR_LDPC_encode(
        a::Vector{Bool},
        R::Rational,
        rv::Integer,
        show_nr_par::Bool;
        E_H = nothing,
        I_LBRM = 0,
        TBS_LBRM = Inf,
        CBGTI = [],
        N_L = 1,
        Q_m = 1
    )

    A = length(a)

    G = round(Int,(A/R)/Q_m)*Q_m

    ### 1) Transport Block CRC attachment (TS38212 Clause 6.2.1, 7.2.1 and 5.1)

    if A > 3824
        # CRC24A:
        L_1 = 24
        p_CRC = "x^24 + x^23 + x^18 + x^17 + x^14 + x^11 + x^10 + x^7 + x^6 + x^5 + x^4 + x^3 + x + 1"
    else
        # CRC16:
        L_1 = 16
        p_CRC = "x^16 + x^12 + x^5 + 1"
    end

    g_CRC = gf2_poly(p_CRC)

    B = A + L_1

    ### 1.5) LDPC base graph selection (TS38212 Clause 7.2.2 and 6.2.2)

    if A ≤ 292 || (A ≤ 3824 && R ≤ 0.67) || R ≤ 0.25
        bg = "2"
    else
        bg = "1"
    end

    ### 2) Code block segmentation and code block CRC attachment
    ### (TS38212 Clauses 7.2.3, 6.2.3 and 5.2.2)

    if bg == "1"
        Kcb = 8448
    else
        Kcb = 3840
    end

    if B ≤ Kcb
        L_2 = 0
        C = 1
        B_prime = B
    else
        L_2 = 24
        C = cld(B,Kcb - L_2)
        B_prime = (B + C*L_2)
    end

    K_prime = B_prime ÷ C # if B ≤ Kcb, K_prime = B_prime = B

    if bg == "1"
        Kb = 22
    else
        if B > 640
            Kb = 10
        elseif B > 560
            Kb = 9
        elseif B > 192
            Kb = 8
        else
            Kb = 6
        end
    end

    Zc, iLS = get_Zc(K_prime, Kb)

    twoZc = 2*Zc

    if K_prime < twoZc
        throw(error(
            """K_prime must be greater than 2 Zc (K_prime = $K_prime, Zc = $Zc)"""
        ))
    end

    b = zeros(Bool,B)

    @inbounds b[1:A] = a
    @inbounds _,b[A+1:end] = divide_poly(b,g_CRC)

    if bg == "1"
        K = 22*Zc
        N = 66*Zc
    else
        K = 10*Zc
        N = 50*Zc
    end  
    
    # Parity bit puncturing size
    P = N - G÷C - K + K_prime

    if show_nr_par
        display("B = $B")
        display("bg = $bg")
        display("L_1 = $L_1")
        display("L_2 = $L_2")
        display("C = $C")
        display("K_prime = $K_prime")
        display("K = $K")
        display("Zc = $Zc")
        display("iLS = $iLS")
        display("N = $N")
        display("G = $G")
        display("P = $P")
    end

    recalc_P = false
    if P > 42*Zc && bg == "1"
        G = C*(N - 42*Zc - K + K_prime)
        recalc_P = true  
    elseif P > 38*Zc && bg == "2"
        G = C*(N - 38*Zc - K + K_prime)
        recalc_P = true
    elseif P < 0
        G = C*(N - K + K_prime)
        recalc_P = true
    end
    
    if recalc_P
        P = N - G÷C - K + K_prime
        R = A//G
        if show_nr_par
            display("new G = $G")
            display("new P = $P")
        end
        display("Atention: NR-LDPC new rate R = $(round(R,digits=3))")
    end

    # calc the length of d
    P_Zc = P÷Zc
    if C == 1
        if bg == "1"
            Cw = zeros(Bool,Zc,68 - P_Zc)
        else
            Cw = zeros(Bool,Zc,52 - P_Zc)
        end
        Cw[1:K_prime] = b
    else
        Ld = N - (P÷Zc)*Zc
        cw = zeros(Int8,Ld+twoZc,C)
        code_block_segmentation!(cw,b,C,K_prime,L_2)
    end

    ### 3) Channel coding (TS38212 Clauses 7.2.4, 6.2.4 and 5.3.2)    

    if E_H === nothing    
        H, E_H = NR_LDPC_make_parity_check_matrix(Zc,iLS,bg)
    else
        H = nothing
    end
    
    W = zeros(Bool,Zc,4)
    Z = zeros(Bool,Zc,4)
    Sc = zeros(Bool,Zc)
    
    if bg == "1"
        J = 22
        I = 46 - P_Zc
    else
        J = 10
        I = 42 - P_Zc
    end   

    if C == 1
        NR_LDPC_parity_bits!(Cw,W,Sc,Z,E_H,I,J)
    else
        for r = 1:C
            Cw .= false
            W .= false
            Sc .= false
            Z .= false
            Cw[1:K,r] = cw[1:K,r]
            NR_LDPC_parity_bits!(Cw,W,Sc,Z,E_H,I,J)
            for k = K_prime + 1 : K
                cw[k,r] = -1
            end
            cw[K+1:end,r] = Cw[K+1:end,r]
        end
    end

    # cw = cw[twoZc+1:end]    

    ### 4) Rate Matching (TS38212 Clauses 7.2.5, 6.2.5 and 5.4.2)

    N_ref = fld(TBS_LBRM,C*R_LBRM)

    if I_LBRM == 0
        N_cb = N
    else
        N_cb = min(N, N_ref)
    end

    E_r = get_Er(C,G,CBGTI,N_L,Q_m)

    k0 = get_k0(rv,Zc,N_cb,bg)

    if C == 1
        e = zeros(Bool,maximum(E_r))
        rate_matching!(e,Cw,twoZc,N_cb,E_r[1],k0,K,K_prime)
    else
        e = zeros(Bool,maximum(E_r),C)
        for r = 1:C
            rate_matching!(e[:,r],cw[:,r],twoZc,N_cb,E_r[r],k0)
        end
    end

    # bit interleaving
    if Q_m ≠ 1
        if C == 1
            f = zeros(Bool,maximum(E_r),1)
            bit_interleaving!(f,e,E_r[1],Q_m)
        else
            f = zeros(Bool,maximum(E_r),C)
            for r = 1:C
                bit_interleaving!(f[:,r],e[:,r],E_r[r],Q_m)
            end
        end
        
    else
        f = e
    end

    ### 5) Code block concatenation (TS38212 Clauses 7.2.6, 6.2.6 and 5.5)
    
    if C == 1
        g = f[:]
    else
        g = zeros(Bool,G)
        code_concatenation!(g,f,G,E_r)
    end

    # test inv functions

    # f_prime = inv_code_concatenation(g,C,E_r)

    # display("f: $(f_prime == f)")

    # e_prime = inv_bit_interleaving(f_prime,C,E_r,Q_m)

    # display("e: $(e_prime == e)")

    # range = (K_prime-twoZc+1):(K-twoZc)

    # d_prime, cw_prime = inv_rate_matching(e_prime,C,Zc,N,N_cb,E_r,k0,range)

    # display("d: $(sum(d_prime[k0+1:k0+E_r[1],1] .=== d[k0+1:k0+E_r[1],1]) === E_r[1])")
    # display("d: $(d_prime[k0+1:k0+E_r[1],1] == d[k0+1:k0+E_r[1],1])")

    if H !== nothing
        H = H[1:end-P,1:end-P]
        H = [H[:,1:K_prime] H[:,K+1:end]]

        if !iszero(H*[a[1:twoZc];g])
            throw(error(
                        lazy"""Wrong encoding"."""
                    ))
        end
    end

    return g, H, E_H, R, nr_ldpc_data(B,K_prime,K,Zc,bg,N_cb,E_r,k0,g_CRC,P_Zc)

end

# KK = 256
# MSG = rand(Bool,KK)
# RR = 1//2

# CWORD, HH, E_H, RR, NR_LDPC_DATA = NR_LDPC_encode(MSG,RR,0,true)
# twoZc = 2*NR_LDPC_DATA.Zc

# iszero(HH*[MSG[1:twoZc];CWORD])


