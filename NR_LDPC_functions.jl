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
    @inbounds for r = 1:C
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
        bg::String,
        E_H::Union{Matrix{<:Integer},Nothing}
    )

    d = zeros(Union{Bool,Missing},N,C)
    cw = zeros(Bool,N+2*Zc,C)
    cw[1:K_prime,:] = c[1:K_prime,:]

    @inbounds for r = 1:C
        for k = 2*Zc +1 : K_prime
            d[k-2*Zc,r] = c[k,r]
        end
        for k = K_prime + 1 : K #(filler bits)
            d[k-2*Zc,r] = missing
        end
    end

    if E_H === nothing    
        H, E_H = make_parity_check_matrix(Zc,iLS,bg)
    else
        H = nothing
    end

    @inbounds for r = 1:C
        w = parity_bits(cw[1:K,r],bg,Zc,K,E_H)
        cw[K+1:end,r] = w
        if H !== nothing && !iszero(H*cw[:,r])
            throw(error(
                    lazy"""Wrong encoding"."""
                ))
        end
        for k = (K + 1) : N + 2*Zc
            d[k - 2*Zc,r] = w[k - K]
        end
    end
    return d, H, E_H

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
    
    @inbounds for r = 1:C
        j = 0
        k = 1
        while k ≤ E_r[r]
            x = rem(k0+j,N_cb)+1
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
        Zc::Integer,
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

    cw = zeros(Bool,N+2*Zc,C)

    @inbounds for r = 1:C        
        j = 0
        k = 1
        while k ≤ E_r[r]
            x = rem(k0+j,N_cb)+1
            if d[x,r] !== missing
                d[x,r] = e[k,r]
                cw[2*Zc + x,r] = e[k,r]
                k += 1
            end
            j += 1
        end

        # if I_HARQ != 0
        #     d[1:N_cb,r] .⊻= d_buffer[:,r]
        #     d_buffer[:,r] = d[1:N_cb,r]
        # end
    end

    return d,cw
end

function 
    bit_interleaving(
        e::Matrix{Bool},
        C::Integer,
        E_r::Vector{<:Integer},
        Q_m::Integer
    )

    f = zeros(Bool,maximum(E_r),C)
    @inbounds for r = 1:C
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

    @inbounds for r = 1:C        
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
    @inbounds while r ≤ C
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
    inv_code_concatenation(
        g::Vector{Bool},
        C::Integer,
        E_r::Vector{<:Integer}
    )

    f = zeros(Bool,maximum(E_r),C)

    k = 1
    r = 1    
    @inbounds while r ≤ C
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
    parity_bits_G(
        c::Vector{Bool},
        Zc::Integer,
        K::Integer,
        H::BitMatrix,
    )

    G = gf2_nullspace(H[1:(4*Zc),1:(K+4*Zc)])
    gf2_reduce!(G)

    p = G[K+1:end,:]*c

    q = H[4*Zc+1:end,1:K+4*Zc]*[c;p]

    return [p;q]

end

function 
    parity_bits(
        c::Vector{Bool},
        bg::String,
        Zc::Integer,
        K::Integer,
        E_H::Matrix{<:Integer}
    )

    if bg == "1"
        J = 22
        N = 68
        M = 46
    else
        J = 10
        N = 52
        M = 42
    end

    cw = zeros(Bool,Zc,N)
    cw[1:K] = c

    a = zeros(Bool,Zc,4)
    Sc = zeros(Bool,Zc)

    @inbounds for i = 1:4
        for j = 1:J            
            if E_H[i,j] ≠ -1
                a[:,i] .⊻= circshift(cw[:,j],-E_H[i,j])
            end
        end
        Sc .⊻= a[:,i]
    end

    @inbounds for i = 1:4
        if E_H[i,J+1] ≠ -1
            cw[:,J+1] .⊻= circshift(Sc,E_H[i,J+1])
        end
    end
    
    z = zeros(Bool,Zc,4)
    @inbounds for i=1:4
        if E_H[i,J+1] ≠ -1
            z[:,i] = circshift(cw[:,J+1],-E_H[i,J+1])
        end
    end

    @inbounds for i = 4:-1:2
        cw[:,J+i] = a[:,i] .⊻ z[:,i] .⊻ cw[:,J+i+1]
    end

    @inbounds for i=5:M
        for j=1:J+4
            if E_H[i,j] ≠ -1
                cw[:,J+i] .⊻= circshift(cw[:,j],-E_H[i,j])
            end
        end
    end    

    return cw[K+1:end]

end

# function circ_mult(
#         e_H::AbstractArray{<:Integer},
#         Zc::Integer,
#         u::Vector{Bool}
#     )

#     L = size(e_H,1)

#     q = zeros(Bool,Zc,L)

#     @inbounds for i in axes(e_H,1)
#         for j in axes(e_H,2)
#             if e_H[i,j] != -1
#                 range = (1 + (j-1)*Zc) : (j*Zc)
#                 q[:,i] .⊻= circshift(u[range],-e_H[i,j])
#             end
#         end
#     end

#     return q
# end

function
    make_parity_check_matrix(
        Zc::Integer,
        iLS::Integer,
        bg::String
    )

    E_H = readdlm("./5G_exponent_matrices/EM_$(bg)_$(iLS)_$(Zc).txt",'\t', Int,'\n')

    m, n = size(E_H)

    H = zeros(Bool,m*Zc, n*Zc)

    I_matrix = Matrix(I(Zc))

    @inbounds for i = 1:m
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