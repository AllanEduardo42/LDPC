################################################################################
# Allan Eduardo Feitosa
# 28 Out 2024
# Function to encode a message using the NR-LDPC

using SparseArrays

include("GF2_poly.jl")
include("make_parity_check_matrix.jl")
include("GF2_functions.jl")

function NR_LDPC_encode(B::Integer,bg::String)

    if bg == "1"
        Kcb = 8448
        Kb = 22
    elseif bg == "2"
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
        throw(ArgumentError(
                lazy"""bg must be "1" or "2"."""
            ))
    end

    b = rand(Bool,B)

    if B <= Kcb
        L = 0
        C = 1
        K0 = B
    else
        L = 24
        C = cld(B,Kcb - L)
        K0 = (B + C*L) ÷ C
    end
    display("L = $L")
    display("C = $C")
    display("K0 = $K0")
    display("Kb = $Kb")

    Z = K0/Kb

    display("Z = $Z")

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

    display("Zc = $Zc")
    display("iLS = $iLS")

    if Kb * Zc < K0
        throw(error(
                lazy"""Kb * Zc must be greater than K0"."""
            ))
    end

    if bg == "1"
        K = 22*Zc
        N = 66*Zc
    else
        K = 10*Zc
        N = 50*Zc
    end
    display("K = $K")
    display("N = $N")

    c = zeros(Bool,K,C)

    # CRC24B:
    # g(x) = x^24 + x^23 + x^6 + x^5 + x + 1

    if C > 1
        g = zeros(Bool,L+1)
        g[L - 24 + 1] = 1
        g[L - 23 + 1] = 1
        g[L -  6 + 1] = 1
        g[L -  5 + 1] = 1
        g[L -  1 + 1] = 1
        g[L -  0 + 1] = 1
    end

    # g[1]*x^L + g[2]*x^(L-1) + ... + g[L]*x + g[L+1] 

    L2 = K0-L
    for r = 1:C
        c[1:L2,r] = b[1 + (r-1)*L2:r*L2]
        if C > 1
            _,c[L2+1:K0,r] = divide_poly([c[1:L2,r];zeros(Bool,L)],g)
        end
    end

    d = zeros(Bool,N,C)

    # for r = 1:C
    #     d[1:K-2*Zc,r] = c[2*Zc+1:K,r]
    # end

    H,_ = make_parity_check_matrix(bg,Zc,iLS)

    G = gf2_nullspace(H)
    gf2_reduce!(G)

    for r=1:C
        x = gf2_mat_mult(G,c[:,r])
        d[:,r] = x[2*Zc+1:end]
    end

    if C == 1
        return b, d[:], H
    else
        return b,d,H
    end
end
