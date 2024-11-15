################################################################################
# Allan Eduardo Feitosa
# 13 NOV 2024
# Function to encode a message using the NR-LDPC

using SparseArrays

include("GF2_poly.jl")
include("GF2_functions.jl")
include("NR_LDPC_functions.jl")

const R_LBRM = 2//3

# mutable struct NR_LDPC
#     E_r::Vector{Int}
#     C::Int
#     N_cb::Int
#     N::Int
#     k0::Int
#     K::Int
#     K_prime::Int
#     Zc::Int
# end

function
    NR_LDPC_encode(
        a::Vector{Bool},
        R::Rational,
        rv::Integer;
        I_LBRM = 0,
        TBS_LBRM = Inf,
        CBGTI = [],
        N_L = 1,
        Q_m = 1
    )

    ### for testing
    # a = rand(Bool,1000)
    # a[1] = true
    # R = 1//5
    # rv = 0
    # I_LBRM = 0
    # TBS_LBRM = Inf
    # CBGTI = []
    # N_L = 1
    # Q_m = 1

    A = length(a)

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

    b = zeros(Bool,B)

    b[1:A] = a
    _,b[A+1:end] = divide_poly(b,g_CRC)

    ### 1.5) LDPC base graph selection (TS38212 Clause 7.2.2 and 6.2.2)

    if A ≤ 292 || (A ≤ 3824 && R ≤ 0.67) || R ≤ 0.25
        bg = "2"
    else
        bg = "1"
    end
    display("bg = $bg")

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
        C = cld(A,Kcb - L_2)
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

    if bg == "1"
        K = 22*Zc
    else
        K = 10*Zc
    end

    display("L_1 = $L_1")
    display("L_2 = $L_2")
    display("C = $C")
    display("K_prime = $K_prime")
    display("Kb = $Kb")
    display("Zc = $Zc")
    display("iLS = $iLS")    
    display("K = $K")    

    c = code_block_segmentation(b,C,K_prime,K,L_2)

    ### 3) Channel coding (TS38212 Clauses 7.2.4, 6.2.4 and 5.3.2)

    if bg == "1"
        N = 66*Zc
    else
        N = 50*Zc
    end
    display("N = $N")

    d, H = channel_coding(c,C,K_prime,K,Zc,iLS,N,bg)

    # d[1:(K_prime-2*Zc),:] == c[(2*Zc+1):K_prime,:] #(payload + CRC)
    # d[(K_prime-2*Zc+1):(K-2*Zc),:] == c[(K_prime+1):K,:] #(filler bits)

    ### 4) Rate Matching (TS38212 Clauses 7.2.5, 6.2.5 and 5.4.2)

    N_ref = fld(TBS_LBRM,C*R_LBRM)

    if I_LBRM == 0
        N_cb = N
    else
        N_cb = min(N, N_ref)
    end

    G = round(Int,(A/R)/Q_m)*Q_m

    display("G = $G")

    E_r = get_Er(C,G,CBGTI,N_L,Q_m)

    k0 = get_k0(rv,Zc,N_cb,bg)

    e =  rate_matching(d,C,N_cb,E_r,k0)

    f =  bit_interleaving(e,C,E_r,Q_m)

    ### 5) Code block concatenation (TS38212 Clauses 7.2.6, 6.2.6 and 5.5)

    g = code_concatenation(f,C,G,E_r)

    # nr_ldpc = NR_LDPC(E_r,C,N_cb,N,k0,K,K_prime,Zc)

    # test inv functions

    f_prime = inv_code_concatenation(g,C,E_r)

    display("f: $(f_prime == f)")

    e_prime = inv_bit_interleaving(f_prime,C,E_r,Q_m)

    display("e: $(e_prime == e)")

    range = (K_prime-2*Zc+1):(K-2*Zc)

    d_prime = inv_rate_matching(e_prime,C,N,N_cb,E_r,k0,range)

    display("d: $(sum(d_prime[1:E_r[1],1] .=== d[1:E_r[1],1]) === E_r[1])")


    P = N - G÷C - K + K_prime
    H = H[1:end-P,1:end-P]
    H = [H[:,1:K_prime] H[:,K+1:end]]

    return H, g, Zc

end



