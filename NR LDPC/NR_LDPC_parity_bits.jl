################################################################################
# Allan Eduardo Feitosa
# 19 May 2029
# Calculate the parity bits of the NR-LDPC encoding

function 
    NR_LDPC_parity_bits(
        c::Vector{Bool},
        bg::String,
        Zc::Integer,
        K::Integer,
        E_H::Matrix{<:Integer},
        P::Integer
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

    @inbounds for i = 5:(M - P÷Zc)
        for j=1:J+4
            if E_H[i,j] ≠ -1
                cw[:,J+i] .⊻= circshift(cw[:,j],-E_H[i,j])
            end
        end
    end    

    return cw[K+1:end]

end