################################################################################
# Allan Eduardo Feitosa
# 25 Mai 2025
# Functions to encode the LDPC code

# Encoding LDPC using base graphs
function 
    encode_LDPC_BG!(
        cword::Vector{Bool},
        Cw::Matrix{Bool},
        Z::Matrix{Bool},
        W::Matrix{Bool},
        sw::Vector{Bool},
        K_prime::Integer,
        K::Integer,
        E_H::Matrix{<:Integer},
        E_M::Integer,
        E_K::Integer,
        E_B::Integer,
        P::Integer,
    )

    calc_LDPC_BG_parity_bits!(Cw,W,sw,Z,E_H,E_M,E_K,E_B)

    @inbounds begin
        cword[1:K_prime] = Cw[1:K_prime]        # systematic bits
        cword[K_prime+1:end] = Cw[K+1:end-P]    # parity bits 
    end
end

# encoding LDPC using LU decomposition
function 
    encode_LDPC_LU!(
        Cw::Vector{Bool},
        H::Matrix{Bool},
        z::Vector{Bool},
        w::Vector{Bool},
        v::Vector{Bool},
        L::Matrix{Bool},
        U::Matrix{Bool},
        K::Integer
    )

    @inbounds begin
        z .= H[:,1:K]*Cw[1:K]
        gf2_solve_LU!(v,L,z)
        gf2_solve_LU!(w,U,v)
        Cw[K+1:end] = w
    end    
end

# parity bits calculation for base graph LDPC
function 
    calc_LDPC_BG_parity_bits!(
        Cw::Matrix{Bool},
        W::Matrix{Bool},
        sw::Vector{Bool},
        Z::Matrix{Bool},
        E_H::Matrix{<:Integer},
        E_M::Integer,
        E_K::Integer,
        E_B::Integer
    )  
    
    sw .= false

    @inbounds begin

        for i in 1:E_B
            W[:,i] .= false
            for j in 1:E_K   
                if E_H[i,j] ≠ -1
                    W[:,i] .⊻= circshift(Cw[:,j],-E_H[i,j])
                end
            end
            sw .⊻= W[:,i]
        end

        for i in 1:E_B
            if E_H[i,E_K+1] ≠ -1
                Cw[:,E_K+1] .⊻= circshift(sw,E_H[i,E_K+1])
            end
        end
        
        for i in 2:E_B  # Z[:,1] is not used
            if E_H[i,E_K+1] ≠ -1
                Z[:,i] = circshift(Cw[:,E_K+1],-E_H[i,E_K+1])
            else
                Z[:,i] .= false
            end
        end

        Cw[:,E_K + E_B] = W[:,E_B] .⊻ Z[:,E_B]
        for i in E_B-1:-1:2
            Cw[:,E_K+i] = W[:,i] .⊻ Z[:,i] .⊻ Cw[:,E_K+i+1]
        end

        if E_B != E_M  # Here WiMAX ≠ NR5G
            for i in 5:E_M
                for j in 1:E_K+E_B
                    if E_H[i,j] ≠ -1
                        Cw[:,E_K+i] .⊻= circshift(Cw[:,j],-E_H[i,j])
                    end
                end
            end
        end
    end

end