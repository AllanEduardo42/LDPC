################################################################################
# Allan Eduardo Feitosa
# 20 Mai 2025
# General encoding procedure for LDPC using LU decomposition
# Original Author: R. C. de Lamare (2007)

function 
    remake_H(
        H::Matrix{Bool},
        strategy::Int
    )::Tuple{Matrix{Bool},Matrix{Bool},Matrix{Bool}}
    
    #  strategy: Strategy for finding the next non-zero diagonal elements
    #            {0} First  : First non-zero found by column search
    #            {1} Mincol : Minimum number of non-zeros in later columns
    #            {2} Minprod: Minimum product of:
    #                         - Number of non-zeros its column minus 1
    #                         - Number of non-zeros its row minus 1

    # Get the matric dimension
    M,N = size(H)
    K = N - M
    @inbounds newH = [H[:,K+1:end] H[:,1:K]]
    # Set a new matrix F for LU decomposition
    F = copy(newH)
    # L matrix
    L = zeros(Bool,M,M)

    # Re-order the M x (N - M) submatrix
    @inbounds for i = 1:M

        # strategy {0 = First 1 = Mincol 2 = Minprod}
        
        # Create diagonally structured matrix using 'First' strategy
        if strategy == 0
            
            # # Find non-zero elements (1s) for the diagonal
            chosenCol = i
            while !F[i,chosenCol]
                chosenCol += 1
            end
            
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
        if chosenCol ≠ i
            tmp = F[:,i]
            F[:,i] = F[:,chosenCol]
            F[:,chosenCol] = tmp

            tmp = newH[:,i]        
            newH[:,i] = newH[:,chosenCol]        
            newH[:,chosenCol] = tmp
        end
                            
        # Fill the L matrix column by column
        L[i:end,i] = F[i:end,i]
                
        # There will be no rows operation at the last row
        if i < M      
            for j in (i + 1):M
                if F[j,i]
                    @. F[j,i:end] ⊻= F[i,i:end]
                end
            end                            
        end            
    end

    # U = F[:,1:M]
    # L*U = newH[1:M]

    @inbounds return [newH[:,M+1:end] newH[:,1:M]], L, F[:,1:M]

end
