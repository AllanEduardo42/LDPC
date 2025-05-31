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
        aux1::Vector{Bool},
        circ_aux::Vector{Bool},
        K::Int,
        N::Int,
        E_H::Matrix{Int},
        eM::Int,
        eK::Int,
        eB::Int,
        S::Int,
        LS::Int
    )

    calc_LDPC_BG_parity_bits!(Cw,W,sw,aux1,circ_aux,E_H,eM,eK,eB,LS)

    @inbounds begin
        for i in 1:K
            cword[i] = Cw[i]      # systematic bits
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
        M::Int,
        K::Int
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
        aux1::Vector{Bool},
        circ_aux::Vector{Bool},
        E_H::Matrix{Int},
        eM::Int,
        eK::Int,
        eB::Int,
        LS::Int
    )  

    # for reusing the vectors sw and circ_aux
    aux2 = sw
    aux3 = circ_aux

    @inbounds begin
    # W[:,i] = ⊻(circshift(Cw[:,j],-E_H[i,j])), 1 ≤ i ≤ eB, 1 ≤ j ≤ eK, if E_H[i,j]≠-1
    # sw = ⊻(W[:,i]), 1 ≤ i ≤ eB
        sw .= false
        for i in 1:eB
            aux1 .= false
            for j in 1:eK   
                if E_H[i,j] ≠ -1
                    my_copyto!(circ_aux,Cw,LS,j)
                    circshift!(circ_aux,-E_H[i,j])
                    my_inplace_xor!(aux1,circ_aux,LS)
                end
            end
            my_copyto!(W,LS,i,aux1)
            my_inplace_xor!(sw,aux1,LS)
        end

    # Cw[:,eK+1] = ⊻(circshift(sw,E[i,eK+1])), 1 ≤ i ≤ eB, if E_H[i,J] ≠ -1
        J = eK + 1
        aux1 .= false
        for i in 1:eB
            if E_H[i,J] ≠ -1
                copy!(circ_aux,sw)
                circshift!(circ_aux,E_H[i,J])
                my_inplace_xor!(aux1,circ_aux,LS)
            end
        end
        my_copyto!(Cw,LS,J,aux1)   


    # Cw[:,eK+i] = circshift(Cw[:,eK+1],-E_H[i,eK+1]), 2 ≤ i ≤ eB, if E_H[i,J] ≠ -1
        aux2 .= false   # now the vector sw can be reused as aux2
        for i in 2:eB
            Ji = eK + i 
            if E_H[i,J] ≠ -1
                copy!(circ_aux,aux1)            # aux1 = Cw[:,eK+1]
                circshift!(circ_aux,-E_H[i,J])
                my_copyto!(Cw,LS,Ji,circ_aux)
            else
                my_copyto!(Cw,LS,Ji,aux2)                
            end
        end

        # Solution for the system of equations:
        # Cw[:,eK +i] = Cw[:,eK +i] ⊻ W[:,i] ⊻ Cw[:,eK+i+1], 1 ≤ i ≤ E_B-1
        # Cw[:,eK+eB] = Cw[:,eK+eB] ⊻ W[:,eB]
        
        JB = eK + eB
        my_copyto!(aux3,Cw,LS,JB)       # we use the vector circ_shift here as aux3
        my_copyto!(aux2,W,LS,eB)        # we use the vector sw here as aux2
        my_inplace_xor!(aux3,aux2,LS)
        my_copyto!(Cw,LS,JB,aux3)

        for i in eB-1:-1:2
            Ji = eK + i
            my_copyto!(aux1,Cw,LS,Ji)
            my_copyto!(aux2,W,LS,i)
            for k in 1:LS
                aux1[k] ⊻= aux2[k]
                aux1[k] ⊻= aux3[k]      # aux3 = Cw[:,eK+i+1]
            end
            my_copyto!(Cw,LS,Ji,aux1)
            copy!(aux3,aux1)      
        end

        if eB != eM  # Here WiMAX ≠ NR5G
            for i in 5:eM
                aux1 .= false
                Ji = eK+i
                for j in 1:eK+eB
                    if E_H[i,j] ≠ -1
                        my_copyto!(circ_aux,Cw,LS,j)
                        circshift!(circ_aux,-E_H[i,j])
                        my_inplace_xor!(aux1,circ_aux,LS)
                    end
                end
                my_copyto!(Cw,LS,Ji,aux1)
            end
        end
    end

end

function 
    my_copyto!(
        A::Matrix{Bool},
        r::Int,
        i::Int,
        x::Vector{Bool}
    )

    copyto!(A,1:r,i:i,'N',x,1:r,1:1)

end

function 
    my_copyto!(
        x::Vector{Bool},
        A::Matrix{Bool},
        r::Int,
        i::Int)

    copyto!(x,1:r,1:1,'N',A,1:r,i:i)

end

function 
    my_inplace_xor!(
        x::Vector{Bool},
        y::Vector{Bool},
        L::Int
    )

    @inbounds for k in 1:L
        x[k] ⊻= y[k]
    end
end

