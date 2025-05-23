################################################################################
# Allan Eduardo Feitosa
# 14 Nov 2024
# Encoding functions for the NR-LDPC encoding

include("NR_LDPC_make_parity_check_matrix.jl")
include("NR_LDPC_auxiliary_functions.jl")

using DelimitedFiles

function 
    code_block_segmentation!(
        cw::Matrix{<:Integer},
        b::Vector{Bool},
        C::Integer,
        K_prime::Integer,
        L::Integer
    )

    if C > 1
        # CRC24B:
        p_CRC =  "x^24 + x^23 + x^6 + x^5 + x + 1"

        g_CRC = gf2_poly(p_CRC)

    end

    s = 1
    @inbounds for r = 1:C
        for k = 1 : (K_prime - L)
            cw[k,r] = b[s]
            s += 1
        end
        if C > 1
            _, p = divide_poly(cw[1:K_prime,r],g_CRC)
            for k = (K_prime - L + 1) : K_prime
                cw[k,r] = p[k + L - K_prime]
            end
        end      
    end
    
end

function 
    NR_LDPC_parity_bits!(
        Cw::Matrix{Bool},
        W::Matrix{Bool},
        Sw::Vector{Bool},
        Z::Matrix{Bool},
        E_H::Matrix{<:Integer},
        I::Integer,
        J::Integer,
    )  
    
    Sw .= false

    # I, _ = size(E_H)

    @inbounds begin

        for i = 1:4
            W[:,i] .= false
            for j = 1:J   
                if E_H[i,j] ≠ -1
                    W[:,i] .⊻= circshift(Cw[:,j],-E_H[i,j])
                end
            end
            Sw .⊻= W[:,i]
        end

        # b1
        for i = 1:4
            if E_H[i,J+1] ≠ -1
                Cw[:,J+1] .⊻= circshift(Sw,E_H[i,J+1])
            end
        end
        
        # z
        for i=1:4
            if E_H[i,J+1] ≠ -1
                Z[:,i] = circshift(Cw[:,J+1],-E_H[i,J+1])
            end
        end

        # B = | I(z1) I 0 0 |
        #     | I(z2) I I 0 |
        #     | I(z3) 0 I I |
        #     | I(z4) 0 0 I |

        Cw[:,J+4] = W[:,4] .⊻ Z[:,4]
        Cw[:,J+3] = W[:,3] .⊻ Z[:,3] .⊻ Cw[:,J+4]
        Cw[:,J+2] = W[:,1] .⊻ Z[:,1]

        for i = 5:I
            for j=1:J+4
                if E_H[i,j] ≠ -1
                    Cw[:,J+i] .⊻= circshift(Cw[:,j],-E_H[i,j])
                end
            end
        end
    end

end

function 
    rate_matching!(
        e::Vector{Bool},
        Cw::Matrix{Bool},
        twoZc::Integer,
        N_cb::Integer,
        E_r::Integer,
        k0::Integer,
        K::Integer,
        K_prime::Integer
    )
    
    @inbounds begin
        j = 0
        k = 1
        while k ≤ twoZc + E_r
            x = rem(k0+j,N_cb + twoZc) + 1
            if x ≤ K_prime || x > K
                e[k] = Cw[x]
                k += 1
            end
            j += 1
        end
    end
end

function 
    bit_interleaving!(
        f::Vector{Bool},
        e::Vector{Bool},
        E_r::Integer,
        Q_m::Integer
    )

    @inbounds begin
        for j = 0:(E_r÷Q_m - 1)
            for i = 0:(Q_m - 1)
                f[i + j*Q_m + 1] = e[i*E_r÷Q_m + j + 1]
            end
        end
    end    

end

function
    code_concatenation!(
        g::Matrix{Bool},
        f::Matrix{Bool},
        C::Integer,
        E_r::Vector{<:Integer},
        twoZc::Integer
    )

    
    k = 1
    r = 1
    @inbounds while r ≤ C
        j = 1
        while j ≤ E_r[r]
            g[k] = f[j+twoZc,r]
            k += 1
            j += 1
        end
        r += 1
    end
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