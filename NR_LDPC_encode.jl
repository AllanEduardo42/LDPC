################################################################################
# Allan Eduardo Feitosa
# 13 NOV 2024
# Function to encode a message using the NR-LDPC

using SparseArrays

include("GF2_poly.jl")
include("GF2_functions.jl")
include("NR_LDPC_functions.jl")

const R_LBRM = 2//3

struct nr_ldpc_data
    A::Integer
    B::Integer
    C::Integer
    G::Integer
    K_prime::Integer
    K::Integer
    L_2::Integer
    Zc::Integer
    iLS::Integer
    N::Integer
    bg::String
    N_cb::Integer
    E_r::Vector{<:Integer}
    k0::Integer
    g_CRC::Vector{Bool}
end


function
    NR_LDPC_encode(
        a::Vector{Bool},
        R::Float64,
        rv::Integer;
        E_H = nothing,
        I_LBRM = 0,
        TBS_LBRM = Inf,
        CBGTI = [],
        N_L = 1,
        Q_m = 1
    )

    ### for testing
    # a = rand(Bool,2064)
    # R = 1/2
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

    # display("B = $B")

    b = zeros(Bool,B)

    @inbounds b[1:A] = a
    @inbounds _,b[A+1:end] = divide_poly(b,g_CRC)

    ### 1.5) LDPC base graph selection (TS38212 Clause 7.2.2 and 6.2.2)

    if A ≤ 292 || (A ≤ 3824 && R ≤ 0.67) || R ≤ 0.25
        bg = "2"
    else
        bg = "1"
    end
    # display("bg = $bg")

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

    # display("L_1 = $L_1")
    # display("L_2 = $L_2")
    # display("C = $C")
    # display("K_prime = $K_prime")
    # display("K = $K")
    # display("Zc = $Zc")
    # display("iLS = $iLS")

    c = code_block_segmentation(b,C,K_prime,K,L_2)

    ### 3) Channel coding (TS38212 Clauses 7.2.4, 6.2.4 and 5.3.2)

    if bg == "1"
        N = 66*Zc
    else
        N = 50*Zc
    end
    # display("N = $N")

    if R == 0
        P = -1
    elseif R == 1
        P = 100*Zc
    else
        G = round(Int,(A/R)/Q_m)*Q_m
        # display("G = $G")
        P = N - G÷C - K + K_prime
        # display("P = $P")
    end

    if P > 42*Zc && bg == "1"
        G = C*(N - 42*Zc - K + K_prime)
        # display("new G = $G")
        P = N - G÷C - K + K_prime
        # display("new P = $P")
        R = A/G
        # display("new R = $(round(R*1000)/1000)")
    elseif P > 38*Zc && bg == "2"
        G = C*(N - 38*Zc - K + K_prime)
        # display("new G = $G")
        P = N - G÷C - K + K_prime
        # display("new P = $P")
        R = A/G
        # display("new R = $(round(R*1000)/1000)")
    end

    if P < 0
        G = C*(N -  K + K_prime)
        # display("new G = $G")
        P = N - G÷C - K + K_prime
        # display("new P = $P")
        R = A/G
        # display("new R = $(round(R*1000)/1000)")
    end

    d, H, E_H = channel_coding(c,C,K_prime,K,Zc,iLS,N,bg,E_H)

    # d[1:(K_prime-2*Zc),:] == c[(2*Zc+1):K_prime,:] #(payload + CRC)
    # d[(K_prime-2*Zc+1):(K-2*Zc),:] == c[(K_prime+1):K,:] #(filler bits)

    ### 4) Rate Matching (TS38212 Clauses 7.2.5, 6.2.5 and 5.4.2)

    N_ref = fld(TBS_LBRM,C*R_LBRM)

    if I_LBRM == 0
        N_cb = N
    else
        N_cb = min(N, N_ref)
    end

    E_r = get_Er(C,G,CBGTI,N_L,Q_m)

    k0 = get_k0(rv,Zc,N_cb,bg)

    e =  rate_matching(d,C,N_cb,E_r,k0)

    f =  bit_interleaving(e,C,E_r,Q_m)

    ### 5) Code block concatenation (TS38212 Clauses 7.2.6, 6.2.6 and 5.5)

    g = code_concatenation(f,C,G,E_r)

    # test inv functions

    # f_prime = inv_code_concatenation(g,C,E_r)

    # display("f: $(f_prime == f)")

    # e_prime = inv_bit_interleaving(f_prime,C,E_r,Q_m)

    # display("e: $(e_prime == e)")

    # range = (K_prime-2*Zc+1):(K-2*Zc)

    # d_prime, cw_prime = inv_rate_matching(e_prime,C,Zc,N,N_cb,E_r,k0,range)

    # display("d: $(sum(d_prime[k0+1:k0+E_r[1],1] .=== d[k0+1:k0+E_r[1],1]) === E_r[1])")
    # display("d: $(d_prime[k0+1:k0+E_r[1],1] == d[k0+1:k0+E_r[1],1])")

    if H !== nothing
        H = H[1:end-P,1:end-P]
        H = [H[:,1:K_prime] H[:,K+1:end]]

        if !iszero(H*[b[1:2*Zc];g])
            throw(error(
                        lazy"""Wrong encoding"."""
                    ))
        end
    end

    return g, H, E_H, nr_ldpc_data(A, B, C, G, K_prime, K, L_2, Zc, iLS, N, bg, N_cb, E_r, k0, g_CRC)

end

function 
    NR_LDPC_encode(
        E_H::Matrix{<:Integer},
        a::Vector{Bool},
        nr_ldpc_data::Any;
        Q_m = 1
    )

    A = nr_ldpc_data.A
    B = nr_ldpc_data.B
    C = nr_ldpc_data.C
    G = nr_ldpc_data.G
    K_prime = nr_ldpc_data.K_prime
    K = nr_ldpc_data.K
    L_2 = nr_ldpc_data.L_2
    Zc = nr_ldpc_data.Zc
    iLS = nr_ldpc_data.iLS
    N = nr_ldpc_data.N
    bg = nr_ldpc_data.bg
    N_cb = nr_ldpc_data.N_cb
    E_r = nr_ldpc_data.E_r
    k0 = nr_ldpc_data.k0
    g_CRC = nr_ldpc_data.g_CRC
    

    b = zeros(Bool,B)
    @inbounds b[1:A] = a
    @inbounds _,b[A+1:end] = divide_poly(b,g_CRC)

    c = code_block_segmentation(b,C,K_prime,K,L_2)

    d, ___ = channel_coding(c,C,K_prime,K,Zc,iLS,N,bg,E_H)

    e =  rate_matching(d,C,N_cb,E_r,k0)

    f =  bit_interleaving(e,C,E_r,Q_m)

    g = code_concatenation(f,C,G,E_r)

    return g

end

# # Message (Payload) size
# A::Int = 2000
# # Rate
# R::Float64 = 1/2
# a = rand(Bool,A)
# rv = 0
# g, H, E_H, A, B, C, G, K_prime, K, L, Zc, iLS, N, bg, N_cb, E_r, k0, g_CRC = NR_LDPC_encode(a,R,rv)
# if !iszero(H*[a[1:2*Zc];g])
#     throw(error(
#                 lazy"""Wrong encoding"."""
#             ))
# end
# for i = 1:10
#     global a = rand(Bool,A)
#     global g = NR_LDPC_encode(E_H,a,A,B,C,G,K_prime,K,L,Zc,iLS,N,bg,N_cb,E_r,k0,g_CRC)
#     if !iszero(H*[a[1:2*Zc];g])
#         throw(error(
#                     lazy"""Wrong encoding"."""
#                 ))
#     end
# end


