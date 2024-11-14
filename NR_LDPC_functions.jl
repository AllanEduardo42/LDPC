################################################################################
# Allan Eduardo Feitosa
# 14 Nov 2024
# Auxiliary functions of the NR-LDPC encoding

using DelimitedFiles

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
    code_block_segmentation(
        b::Vector{Bool},
        C::Integer,
        K_prime::Integer,
        K::Integer,
        L::Integer
    )

    if C > 1
        # CRC24B:
        p_CRC =  "x^24 + x^23 + x^6 + x^5 + x + 1"

        g_CRC = gf2_poly(p_CRC)

    end

    c = zeros(Union{Bool,Missing},K,C)
    s = 1
    for r = 1:C
        for k = 1 : (K_prime - L)
            c[k,r] = b[s]
            s += 1
        end
        if C > 1
            _, p = divide_poly(c[1:K_prime,r],g_CRC)
            for k = (K_prime - L + 1) : K_prime
                c[k,r] = p[k + L - K_prime]
            end
        end
        for k = (K_prime + 1) : K #(filler bits)
            c[k,r] = missing
        end        
    end

    return c
    
end

function 
    channel_coding(
        c::Matrix{Union{Bool,Missing}},
        C::Integer,
        K_prime::Integer,
        K::Integer,
        Zc::Integer,
        iLS::Integer,
        N::Integer,
        bg::String
    )

    d = zeros(Union{Bool,Missing},N,C)

    for r = 1:C
        for k = 2*Zc +1 : K_prime
            d[k-2*Zc,r] = c[k,r]
        end
        for k = K_prime + 1 : K #(filler bits)
            c[k,r] = false
            d[k-2*Zc,r] = missing
        end
    end
    
    H, E_H = make_parity_check_matrix(Zc,iLS,bg)

    for r = 1:C
        w = parity_bits(c[:,r],Zc,E_H,bg)
        for k = (K + 1) : N + 2*Zc
            d[k - 2*Zc,r] = w[k - K]
        end
        if !iszero(gf2_mat_mult(H,[Bool.(c[1:K,r]); w]))
            throw(error(
                    lazy"""Wrong encoding"."""
                ))
        end
    end

    return d, H

end


function 
    rate_matching(
        d::Matrix{Union{Bool,Missing}},
        C::Integer,
        N_cb::Integer,
        E_r::Vector{<:Integer},
        k0::Integer
    )
    e = zeros(Bool,maximum(E_r),C)
    
    for r = 1:C
        j = 1
        k = 1
        while k ≤ E_r[r]
            x = rem(k0+j,N_cb)
            if d[x,r] !== missing
                e[k,r] = d[x,r]
                k += 1
            end
            j += 1
        end
    end    

    return e
end

function 
    inv_rate_matching(
        e::Matrix{Bool},
        C::Integer,
        N::Integer,
        N_cb::Integer,
        E_r::Vector{<:Integer},        
        k0::Integer,
        range::UnitRange{Int}
        # I_HARQ,
        # d_buffer
    )

    d = zeros(Union{Bool,Missing},N,C)
    d[range,:] .= missing

    for r = 1:C        
        j = 1
        k = 1
        while k ≤ E_r[r]
            x = rem(k0+j,N_cb)
            if d[x,r] !== missing
                d[x,r] ⊻= e[k,r]
                k += 1
            end
            j += 1
        end

        # if I_HARQ != 0
        #     d[1:N_cb,r] .⊻= d_buffer[:,r]
        #     d_buffer[:,r] = d[1:N_cb,r]
        # end
    end

    return d
end

function 
    bit_interleaving(
        e::Matrix{Bool},
        C::Integer,
        E_r::Vector{<:Integer},
        Q_m::Integer
    )

    f = zeros(Bool,maximum(E_r),C)
    for r = 1:C
        for j = 0:(E_r[r]÷Q_m - 1)
            for i = 0:(Q_m - 1)
                f[i + j*Q_m + 1,r] = e[i*E_r[r]÷Q_m + j + 1,r]
            end
        end
    end    

    return f
end

function 
    inv_bit_interleaving(
        f::Matrix{Bool},
        C::Integer,
        E_r::Vector{<:Integer},
        Q_m::Integer
    )

    e = zeros(Bool,maximum(E_r),C)

    for r = 1:C        
        for j = 0:(E_r[r]÷Q_m - 1)
            for i = 0:(Q_m - 1)
                e[i*E_r[r]÷Q_m + j + 1,r] = f[i + j*Q_m + 1,r]
            end
        end
    end

    return e
end

function
    code_concatenation(
        f::Matrix{Bool},
        C::Integer,
        G::Integer,
        E_r::Vector{<:Integer}
    )

    g = zeros(Bool,Int(G))
    k = 1
    r = 1
    while r ≤ C
        j = 1
        while j ≤ E_r[r]
            g[k] = f[j,r]
            k += 1
            j += 1
        end
        r += 1
    end

    return g
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
    for r = 1:C
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
    inv_code_concatenation(
        g::Vector{Bool},
        C::Integer,
        E_r::Vector{<:Integer}
    )

    f = zeros(Bool,maximum(E_r),C)

    k = 1
    r = 1    
    while r ≤ C
        j = 1
        while j ≤ E_r[r]
            f[j,r] = g[k]
            k += 1
            j += 1
        end
        r += 1
    end

    return f
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

function 
    parity_bits(
        c::Vector{Union{Bool,Missing}},
        Zc::Integer,
        E_H::Matrix{<:Integer},
        bg::String
    )

    if bg == "1"
        I = 22
        J = 68
        M = 46
    else
        I = 10
        J = 52
        M = 42
    end

    E_A = E_H[1:4,1:I]
    E_C = E_H[5:end,1:I+4]

    a1 = circ_mult(E_A[1,:],I,Zc,c)
    a2 = circ_mult(E_A[2,:],I,Zc,c)
    a3 = circ_mult(E_A[3,:],I,Zc,c)
    a4 = circ_mult(E_A[4,:],I,Zc,c)

    p1 = a1 .⊻ a2 .⊻ a3 .⊻ a4

    if bg == "1"
        zp1 = circshift(p1,-1)
        p2 = circ_mult(E_H[1,:],I,Zc,c) .⊻ zp1
        p4 = circ_mult(E_H[4,:],I,Zc,c) .⊻ zp1
        p3 = circ_mult(E_H[3,:],I,Zc,c) .⊻ p4
        
    else
        p1 = circshift(p1,1)
        p2 = circ_mult(E_H[1,:],I,Zc,c) .⊻ p1
        p3 = circ_mult(E_H[2,:],I,Zc,c) .⊻ p2
        p4 = circ_mult(E_H[4,:],I,Zc,c) .⊻ p1
    end

    p = [p1;p2;p3;p4]

    c = [c;p] # K + 4*Zc

    q = zeros(Bool,(J-I-4)*Zc) # N - K - 4*Zc

    for m = 1:M-4

        q[1+Zc*(m-1):Zc*m] .= circ_mult(E_C[m,:],I+4,Zc,c)

    end

    return [p;q]

end

function circ_mult(
    e_A::AbstractArray{<:Integer},
    I::Integer,
    Zc::Integer,
    c::Vector{Union{Bool,Missing}})

    q = zeros(Bool,Zc)

    for i = 1:I
        if e_A[i] != -1
            range = (1 + (i-1)*Zc) : (i*Zc)
            q .⊻= circshift(c[range],-e_A[i])
        end
    end

    return q
end

function
    make_parity_check_matrix(
        Zc::Integer,
        iLS::Integer,
        bg::String
    )

    E_H = readdlm("./exponent_matrices/EM_$(bg)_$(iLS)_$(Zc).txt",'\t', Int,'\n')

    m, n = size(E_H)

    H = zeros(Bool,m*Zc, n*Zc)

    I_matrix = Matrix(I(Zc))

    for i = 1:m
        for j = 1:n 
            row_range = Zc*(i-1)+1 : Zc*i
            col_range = Zc*(j-1)+1 : Zc*j
            if E_H[i,j] != -1
                H[row_range,col_range] = circshift(I_matrix,-E_H[i,j])
            end  
        end   
    end

    return BitMatrix(H), E_H

end