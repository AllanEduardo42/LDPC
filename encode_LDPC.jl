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
        aux1::Vector{Bool},
        aux2::Vector{Bool},
        aux3::Vector{Bool},
        K::Int,
        N::Int,
        E_H::Matrix{Int},
        eM::Int,
        eK::Int,
        eB::Int,
        S::Int,
        Ls::Int
    )

    calc_LDPC_BG_parity_bits!(Cw,W,aux1,aux2,aux3,E_H,eM,eK,eB,Ls)

    @inbounds begin
        for i in 1:K
            cword[i] = Cw[i]      # systematic bits
        end
        for i = K+1:N
            cword[i] = Cw[S+i]    # parity bits (remove the filler bits)
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
        aux1::Vector{Bool},
        aux2::Vector{Bool},
        aux3::Vector{Bool},
        E_H::Matrix{Int},
        eM::Int,
        eK::Int,
        eB::Int,
        Ls::Int
    )

    @inbounds begin
    # W[:,i] = ⊻(circshift(Cw[:,j],-E_H[i,j])), 1 ≤ i ≤ eB, 1 ≤ j ≤ eK, if E_H[i,j]≠-1
    # aux3 = circshift(Cw[:,j],-E_H[i,j])
    # aux2 = W[:,i] = ⊻(aux3)
    # aux1 = ⊻(W[:,i]) = ⊻(aux2)
        aux1 .= false
        for i in 1:eB
            aux2 .= false
            for j in 1:eK   
                if E_H[i,j] ≠ -1
                    my_copyto!(aux3,Cw,Ls,j)
                    circshift!(aux3,-E_H[i,j])
                    my_inplace_xor!(aux2,aux3,Ls)
                end
            end
            my_copyto!(W,Ls,i,aux2)
            my_inplace_xor!(aux1,aux2,Ls)
        end
    
    # Cw[:,eK+1] = ⊻(circshift(aux1,E[i,eK+1])), 1 ≤ i ≤ eB, if E_H[i,J] ≠ -1
    # aux1 = ⊻(W[:,i]), 1 ≤ i ≤ eB
    # aux3 = circshift(aux1,E[i,eK+1])
    # aux2 = ⊻(aux3) = Cw[:,eK+1] 
        J = eK + 1
        aux2 .= false
        for i in 1:eB
            if E_H[i,J] ≠ -1
                copy!(aux3,aux1)
                circshift!(aux3,E_H[i,J])
                my_inplace_xor!(aux2,aux3,Ls)
            end
        end
        my_copyto!(Cw,Ls,J,aux2)

    # Cw[:,eK+i] = circshift(Cw[:,eK+1],-E_H[i,eK+1]), 2 ≤ i ≤ eB, if E_H[i,J] ≠ -1
    # aux2 = Cw[:,eK+i]
    # aux3 = circshift(aux1,-E_H[i,eK+1])
    # aux2 = zeros(Bool,Ls)
        aux1 .= false
        for i in 2:eB
            Ji = eK + i 
            if E_H[i,J] ≠ -1
                copy!(aux3,aux2)            # aux2 = Cw[:,eK+1]
                circshift!(aux3,-E_H[i,J])
                my_copyto!(Cw,Ls,Ji,aux3)
            else
                my_copyto!(Cw,Ls,Ji,aux1)                
            end
        end

        # Solution for the system of equations:
        # Cw[:,eK +i] = Cw[:,eK +i] ⊻ W[:,i] ⊻ Cw[:,eK+i+1], 1 ≤ i ≤ E_B-1
        # Cw[:,eK+eB] = Cw[:,eK+eB] ⊻ W[:,eB]
        
        # aux3 = Cw[:,eK+eB]
        # aux2 = W[:,eB]
        JB = eK + eB    
        my_copyto!(aux2,W,Ls,eB) 
        my_copyto!(aux3,Cw,Ls,JB)          
        my_inplace_xor!(aux3,aux2,Ls)
        my_copyto!(Cw,Ls,JB,aux3)

        # for E_B-1 ≥ i ≥ 1  :
        # aux3 = Cw[:,eK+i+1]
        # aux2 = W[:,i]
        # aux1 = Cw[:,eK+i]
        for i in eB-1:-1:2
            Ji = eK + i
            my_copyto!(aux1,Cw,Ls,Ji)
            my_copyto!(aux2,W,Ls,i)
            my_inplace_xor!(aux1,aux2,Ls)
            my_inplace_xor!(aux1,aux3,Ls)
            my_copyto!(Cw,Ls,Ji,aux1)
            copy!(aux3,aux1)      
        end

        if eB != eM  # Here WiMAX ≠ NR5G
            for i in 5:eM
                aux1 .= false
                Ji = eK+i
                for j in 1:eK+eB
                    if E_H[i,j] ≠ -1
                        my_copyto!(aux3,Cw,Ls,j)
                        circshift!(aux3,-E_H[i,j])
                        my_inplace_xor!(aux1,aux3,Ls)
                    end
                end
                my_copyto!(Cw,Ls,Ji,aux1)
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

