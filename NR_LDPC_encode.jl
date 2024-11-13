################################################################################
# Allan Eduardo Feitosa
# 28 Out 2024
# Function to encode a message using the NR-LDPC

using SparseArrays

include("GF2_poly.jl")
include("make_parity_check_matrix.jl")
include("GF2_functions.jl")
include("rate_matching.jl")
include("encode.jl")

const R_LBRM = 2//3

function
    NR_LDPC_encode(
        a::Vector{Bool},
        R::Rational;
        I_LBRM = 0,
        TBS_LBRM = Inf,
        rv = 0,
        CBGTI = [],
        N_L = 1,
        Q_m = 1
    )   

    # a = rand(Bool,2000)
    # R = Float64(1//2)
    # I_LBRM = 0
    # TBS_LBRM = Inf
    # rv = 0
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

    if b[1:A] ≠ a
        throw(error(
                lazy"""Encode error at 1."""
            ))
    end

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

    # c = Matrix{Union{Bool,Nothing}}(undef,K,C)
    c = zeros(Bool,K,C)

    if C > 1

        # CRC24B:
        p_CRC =  "x^24 + x^23 + x^6 + x^5 + x + 1"

        g_CRC = gf2_poly(p_CRC)

    end

    s = 1
    for r = 1:C
        for k = 1 : (K_prime - L_2)
            c[k,r] = b[s]
            s += 1
        end
        if C > 1
            _, p = divide_poly(c[1:K_prime,r],g_CRC)
            for k = (K_prime - L_2 + 1) : K_prime
                c[k,r] = p[k + L_2 - K_prime]
            end
        end
        # for k = (K_prime + 1) : K #(filler bits)
        #     c[k,r] = nothing
        # end        
    end

    b_prime = c[1:K_prime - L_2,1]

    for r = 2:C
        b_prime = [b_prime;c[1:K_prime-L_2,r]]
    end
    if b_prime ≠ b
        throw(error(
                lazy"""Encode error at 2 (b)."""
            ))
    end

    ### 3) Channel coding (TS38212 Clauses 7.2.4, 6.2.4 and 5.3.2)

    if bg == "1"
        N = 66*Zc
    else
        N = 50*Zc
    end
    display("N = $N")

    d = Matrix{Union{Bool,Nothing}}(undef,N,C)

    for r = 1:C
        for k = 2*Zc +1 : K_prime
            d[k-2*Zc,r] = c[k,r]
        end
        for k = K_prime + 1 : K #(filler bits)
            d[k-2*Zc,r] = nothing
        end
    end

    # d[1:(K_prime-2*Zc),:] == c[(2*Zc+1):K_prime,:] #(payload + CRC)
    # d[(K_prime-2*Zc+1):(K-2*Zc),:] == c[(K_prime+1):K,:] #(filler bits)

    H,E_H = make_parity_check_matrix(bg,Zc,iLS)

    for r = 1:C
        global w = parity_bits(E_H,c[:,r],Zc,bg)
        for k = (K + 1) : N + 2*Zc
            d[k - 2*Zc,r] = w[k - K]
        end
        if !iszero(gf2_mat_mult(H,[c[:,r];w]))
            throw(error(
                    lazy"""Wrong encoding"."""
                ))
        end
    end

    c_prime = zeros(Bool,K,C)
    for r = 1:C
        c_prime[1:(2*Zc),r] = c[1:(2*Zc),r]
        for k=(2*Zc+1):K
            if d[k - 2*Zc,r] === nothing
                c_prime[k,r] = false
            else
                c_prime[k,r] = d[k - 2*Zc,r]
            end
        end
    end
    if c_prime ≠ c
        throw(error(
                lazy"""Encode error at 3."""
            ))
    end

    ### 4) Rate Matching (TS38212 Clauses 7.2.5, 6.2.5 and 5.4.2)

    N_ref = fld(TBS_LBRM,C*R_LBRM)

    if I_LBRM == 0
        N_cb = N
    else
        N_cb = min(N, N_ref)
    end

    G = round(Int,(A/R)/Q_m)*Q_m

    E_r = get_Er(C,CBGTI,G,N_L,Q_m)

    k0 = get_k0(rv,bg,N_cb,Zc)

    e =  rate_matching(E_r,d,N_cb,k0,C)

    f =  bit_interleaving(E_r,Q_m,e,C)

    ### 5) Code block concatenation (TS38212 Clauses 7.2.6, 6.2.6 and 5.5)

    g = code_concatenation(G,E_r,f,C)

    # x = c[1:2*Zc,1]
    # x = [x;g[1:K_prime-2*Zc]]
    # x = [x; zeros(Bool,K-K_prime)]
    # x = [x;g[K_prime-2*Zc+1:end]]

    # Nx = length(x)
    # # Mx = Int(Nx - length(g)*R)
    # Mx = Nx - K
    # iszero(gf2_mat_mult(H[1:Mx,1:Nx],x))

    # x = c[1:2*Zc,1]
    # x = [x;g]
    # Nx = length(x)
    # Mx = Nx - K
    # Hx = H[1:Mx,1:K_prime]
    # Hx = [Hx H[1:Mx,K+1:length(g)+(2*Zc+(K-K_prime))]]
    # iszero(gf2_mat_mult(Hx,x))

    return H, E_r, C, g, N_cb, N, k0, K, K_prime, Zc, d

end



