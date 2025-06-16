################################################################################
# Allan Eduardo Feitosa
# 19 May 2029
# Make the Parity Check Matrix the NR-LDPC encoding

function
    NR_LDPC_make_parity_check_matrix(
        Zc::Int,
        iLS::Int,
        bg::String,
        P::Int,
        K_prime::Int,
        K::Int,
        P_Zc::Int
    )::Tuple{Matrix{Bool},Matrix{Int}}

    E_H = readdlm("./5G_exponent_matrices/EM_$(bg)_$(iLS)_$(Zc).txt",'\t', Int,'\n')

    m, n = size(E_H)

    H = zeros(Bool,m*Zc, n*Zc)

    I_matrix = Matrix(I(Zc))

    @inbounds begin
        for i = 1:m
            for j = 1:n 
                row_range = Zc*(i-1)+1 : Zc*i
                col_range = Zc*(j-1)+1 : Zc*j
                if E_H[i,j] != -1
                    H[row_range,col_range] = circshift(I_matrix,-E_H[i,j])
                end  
            end   
        end

        # Puncturing 
        H = H[1:end-P,1:end-P]
        H = [H[:,1:K_prime] H[:,K+1:end]]

        J1 = cld(K_prime,Zc)
        J2 = K÷Zc

        E_H = E_H[1:end-P_Zc,1:end-P_Zc]
        E_H = [E_H[:,1:J1] E_H[:,J2+1:end]]
    end

    return H, E_H

end