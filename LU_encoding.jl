################################################################################
# Allan Eduardo Feitosa
# 20 Mai 2025
# General encoding procedure for LDPC using LU decomposition
# Original Author: R. C. de Lamare (2007)

function 
    remake_H(
        H::Matrix{Bool},
        strategy::Integer
    )
    
    #  strategy: Strategy for finding the next non-zero diagonal elements
    #            {0} First  : First non-zero found by column search
    #            {1} Mincol : Minimum number of non-zeros in later columns
    #            {2} Minprod: Minimum product of:
    #                         - Number of non-zeros its column minus 1
    #                         - Number of non-zeros its row minus 1

    # Get the matric dimension
    M,N = size(H)
    newH = copy(H)
    # Set a new matrix F for LU decomposition
    F = copy(newH)
    # LU matrices
    L = zeros(Bool,M,M)
    U = zeros(Bool,M,M)

    # Re-order the M x (N - M) submatrix
    @inbounds for i = 1:M

        # strategy {0 = First 1 = Mincol 2 = Minprod}
        
        # Create diagonally structured matrix using 'First' strategy
        if strategy == 0
            
            # Find non-zero elements (1s) for the diagonal
            y = findall(F[:,i:end])
            r = getindex.(y,1)
            c = getindex.(y,2)
            
            # Find non-zero diagonal element candidates
            rowIndex = findall(r .== i)
            
            # Find the first non-zero column
            chosenCol = c[rowIndex[1]] + (i - 1)
            
        # Create diagonally structured matrix using 'Mincol' strategy
        elseif strategy == 1
            
            # Find non-zero elements (1s) for the diagonal
            y = findall(F[:, i:end])
            r = getindex.(y,1)
            c = getindex.(y,2)
            colWeight = sum(F[:, i:end],dims=1)
            
            # Find non-zero diagonal element candidates
            rowIndex = findall(r .== i)
            
            # Find the minimum column weight
            x, ix = findmin(colWeight[c[rowIndex]])
            # Add offset to the chosen row index to match the dimension of the... 
            # original matrix F
            chosenCol = c[rowIndex(ix)] + (i - 1)
                
        # Create diagonally structured matrix using 'Minprod' strategy   
        elseif strategy == 2
            
            # Find non-zero elements (1s) for the diagonal
            y = findall(F[:,i:end])
            r = getindex.(y,1)
            c = getindex.(y,2)
            colWeight = sum(F[:, i:end],dims=1) - 1
            rowWeight = sum(F[i, :]) - 1
            
            # Find non-zero diagonal element candidates
            rowIndex = findall(r .== i)
            
            # Find the minimum product
            x, ix = minimum(colWeight[c[rowIndex]]*rowWeight)
            # Add offset to the chosen row index to match the dimension of the... 
            # original matrix F
            chosenCol = c[rowIndex[x]] + (i - 1)
            
        else
            println("Please select a correct columns re-ordering strategy!")
        end 

        # Re-ordering columns of both newH and F
        tmp1 = F[:,i]
        tmp2 = newH[:,i]
        F[:,i] = F[:,chosenCol]
        newH[:,i] = newH[:,chosenCol]
        F[:,chosenCol] = tmp1
        newH[:,chosenCol] = tmp2
                            
        # Fill the LU matrices column by column
        L[i:end,i] = F[i:end,i]
        U[1:i,i] = F[1:i,i]
                
        # There will be no rows operation at the last row
        if i < M
            # Find the later rows with non-zero elements in column i
            r2 = findall(F[(i + 1):end, i])       
            # Add current row to the later rows which have a 1 in column i
            for j in r2
                @. F[i+j,:] ⊻= F[i,:]
            end                            
        end
            
    end

    K = N - M

    display("K = $K, M = $M, N = $N")

    @inbounds newH = [newH[:,M+1:end] newH[:,1:M]]

    return newH, L, U

end

function 
    LU_parity_bits!(
        cword::Vector{Bool},
        H::Matrix{Bool}, 
        L::Matrix{Bool}, 
        U::Matrix{Bool}, 
        K::Integer
    )

    @inbounds begin
        z = H[:,1:K]*cword[1:K]
        # Parity check vector found by solving sparse LU
        cword[K+1:end] = gf2_solve_LU(U,gf2_solve_LU(L,z))
    end
end
