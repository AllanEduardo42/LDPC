function 
    rate_matching(
        E_r::Vector{<:Integer},
        d::Matrix{Union{Bool,Nothing}},
        N_cb::Integer,
        k0::Integer,
        C::Integer
    )

    e = zeros(Bool,maximum(E_r),C)
    
    for r = 1:C
        j = 1
        k = 1
        while k ≤ E_r[r]
            x = rem(k0+j,N_cb)
            if d[x,r] !== nothing
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
        E_r,
        N::Integer,
        C,
        e,
        N_cb::Integer,
        k0,
        K,
        K_prime,
        Zc
        # I_HARQ,
        # d_buffer
    )

    d = Matrix{Union{Bool,Nothing}}(undef,N,C)
    d .= false
    d[(K_prime-2*Zc+1):(K-2*Zc),:] .= nothing

    for r = 1:C        
        j = 1
        k = 1
        while k ≤ E_r[r]
            x = rem(k0+j,N_cb)
            if d[x,r] !== nothing
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

function get_Er(C, CBGTI, G, N_L, Q_m)

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
            j = j + 1
        end
    end
    return E_r
end

function 
    get_CBGTI_flags(C, CBGTI)
    
    CBGTI_flags = ones(Bool,C)
    CBGTI_flags[CBGTI[CBGTI .< C] .+ 1] .= false

    return CBGTI_flags
end

function 
    bit_interleaving(E_r,Q_m,e,C)

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
    inv_bit_interleaving(C,E_r,Q_m,f)

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
    code_concatenation(G,E_r,f,C)

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
    inv_code_concatenation(C,E_r,g)

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

function get_Zc(K_prime, Kb)

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
    get_k0(
        rv::Integer,
        bg::String,
        N_cb::Integer,
        Zc::Integer
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