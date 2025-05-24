################################################################################
# Allan Eduardo Feitosa
# 24 Mai 2025
# Calculate the parity bits from base parity check matrix (E_H)

function 
    calc_parity_bits!(
        Cw::Matrix{Bool},
        W::Matrix{Bool},
        Sw::Vector{Bool},
        Z::Matrix{Bool},
        E_H::Matrix{<:Integer},
        E_M::Integer,
        E_K::Integer,
        E_B::Integer
    )  
    
    Sw .= false

    @inbounds begin

        for i in 1:E_B
            W[:,i] .= false
            for j in 1:E_K   
                if E_H[i,j] ≠ -1
                    W[:,i] .⊻= circshift(Cw[:,j],-E_H[i,j])
                end
            end
            Sw .⊻= W[:,i]
        end

        for i in 1:E_B
            if E_H[i,E_K+1] ≠ -1
                Cw[:,E_K+1] .⊻= circshift(Sw,E_H[i,E_K+1])
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