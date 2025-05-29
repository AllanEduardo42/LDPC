################################################################################
# Allan Eduardo Feitosa
# 25 Mai 2025
# Functions to encode the LDPC code

# Encoding LDPC using base graphs
function 
    encode_LDPC_BG!(
        cword::Vector{Bool},
        Cw::Matrix{Bool},
        W::Matrix{Bool},
        sw::Vector{Bool},
        circ_aux::Vector{Bool},
        K::Integer,
        N::Integer,
        E_H::Matrix{<:Integer},
        E_M::Integer,
        E_K::Integer,
        E_B::Integer,
        S::Integer,
    )

    calc_LDPC_BG_parity_bits!(Cw,W,sw,circ_aux,E_H,E_M,E_K,E_B)

    @inbounds begin
        for i in 1:K
            cword[i] = Cw[i]        # systematic bits
        end
        for i = K+1:N
            cword[i] = Cw[S+i]    # parity bits 
        end
    end
end

# encoding LDPC using LU decomposition
function 
    encode_LDPC_LU!(
        cw::Vector{Bool},
        H1::Matrix{Bool},
        w::Vector{Bool},
        L::Matrix{Bool},
        U::Matrix{Bool},
        M::Integer,
        K::Integer
    )

    _gf2_mat_mult!(w,H1,cw,M,K)
    _gf2_solve_LU!(w,L,U,w,M) 
end

# parity bits calculation for base graph LDPC
function 
    calc_LDPC_BG_parity_bits!(
        Cw::Matrix{Bool},
        W::Matrix{Bool},
        sw::Vector{Bool},
        circ_aux::Vector{Bool},
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
                    circ_aux .= Cw[:,j]
                    circshift!(circ_aux,-E_H[i,j])
                    W[:,i] .⊻= circ_aux
                end
            end
            sw .⊻= W[:,i]
        end

        Cw[:,E_K+1] .= false
        for i in 1:E_B
            if E_H[i,E_K+1] ≠ -1
                circ_aux .= sw
                circshift!(circ_aux,E_H[i,E_K+1])
                Cw[:,E_K+1] .⊻= circ_aux
            end
        end
        
        for i in 2:E_B
            if E_H[i,E_K+1] ≠ -1
                circ_aux .= Cw[:,E_K+1]
                circshift!(circ_aux,-E_H[i,E_K+1])
                Cw[:,E_K+i] = circ_aux
            else
                Cw[:,E_K+i] .= false
            end
        end

        Cw[:,E_K + E_B] .⊻= W[:,E_B]
        for i in E_B-1:-1:2
            Cw[:,E_K+i] .⊻= W[:,i]
            Cw[:,E_K+i] .⊻= Cw[:,E_K+i+1]
        end

        if E_B != E_M  # Here WiMAX ≠ NR5G
            for i in 5:E_M
                Cw[:,E_K+i] .= false
                for j in 1:E_K+E_B
                    if E_H[i,j] ≠ -1
                        circ_aux .= Cw[:,j]
                        circshift!(circ_aux,-E_H[i,j])
                        Cw[:,E_K+i] .⊻= circ_aux
                    end
                end
            end
        end
    end

end