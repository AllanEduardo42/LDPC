################################################################################
# Allan Eduardo Feitosa
# 13 NOV 2024
# Function to encode a message using the NR-LDPC

using SparseArrays

include("../get_CRC_poly.jl")
include("../GF2_functions.jl")
include("NR_LDPC_auxiliary_functions.jl")

const R_LBRM = 2//3

struct nr_ldpc_data
    A::Int
    B::Int
    C::Int
    bg::String
    K_prime::Int
    K::Int
    N::Int
    Zc::Int
    iLS::Int
    P::Int
    P_Zc::Int
    N_cb::Int
    E_r::Vector{Int}
    k0::Int    
end

function 
    NR_LDPC_parameters(
        G::Int,
        R::Float64,
        rv::Int,
        show_nr_par::Bool;
        Q_m = 1,
        I_LBRM = 0,
        TBS_LBRM = Inf,
        CBGTI = [],
        N_L = 1,
    )::Tuple{Int,Int,Float64,Vector{Bool},Int,nr_ldpc_data}

    A = get_TBS(G,R)

    B, g_CRC = get_CRC_poly(A)   

    if A ≤ 292 || (A ≤ 3824 && R ≤ 0.67) || R ≤ 0.25
        bg = "2"
        Kcb = 3840
        if B > 640
            Kb = 10
        elseif B > 560
            Kb = 9
        elseif B > 192
            Kb = 8
        else
            Kb = 6
        end
    else
        bg = "1"
        Kcb = 8448
        Kb = 22
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

    Zc, iLS = get_Zc(K_prime, Kb)

    twoZc = 2*Zc

    if K_prime < twoZc
        throw(error(
            """K_prime must be greater than 2 Zc (K_prime = $K_prime, Zc = $Zc)"""
        ))
    end

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
        display("A = $A")
        display("B = $B")
        display("bg = $bg")
        display("C = $C")
        display("K_prime = $K_prime")
        display("K = $K")
        display("Zc = $Zc")
        display("iLS = $iLS")
        display("N = $N")
        display("G = $G")
        display("P = $P")
    end

    # recalc_P = false
    # if P > 42*Zc && bg == "1"
    #     G = C*(N - 42*Zc - K + K_prime)
    #     recalc_P = true  
    # elseif P > 38*Zc && bg == "2"
    #     G = C*(N - 38*Zc - K + K_prime)
    #     recalc_P = true
    # elseif P < 0
    #     G = C*(N - K + K_prime)
    #     recalc_P = true
    # end
    
    # if recalc_P
    #     P = N - G÷C - K + K_prime
    #     R = A//G
    #     if show_nr_par
    #         display("new G = $G")
    #         display("new P = $P")
    #     end
    #     display("Atention: NR-LDPC new rate R = $(round(R,digits=3))")
    # end

    if P > 42*Zc && bg == "1" || P > 38*Zc && bg == "2"
        throw(error(
            """rate is too high or G is too low"""
        ))
    end
    if P < 0
        throw(error(
            """rate is too low"""
        ))
    end

    P_Zc = P÷Zc

    N_ref = fld(TBS_LBRM,C*R_LBRM)

    if I_LBRM == 0
        N_cb = N
    else
        N_cb = min(N, N_ref)
    end

    E_r = get_Er(C,G,CBGTI,N_L,Q_m)

    k0 = get_k0(rv,Zc,N_cb,bg)

    return A, K_prime, R, g_CRC, Zc, nr_ldpc_data(A,B,C,bg,K_prime,K,N,Zc,iLS,P,P_Zc,N_cb,E_r,k0) 

end


