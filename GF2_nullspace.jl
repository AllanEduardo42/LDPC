################################################################################
# Allan Eduardo Feitosa
# 4 set 2024
# find GF2 matrix nullspace

using LinearAlgebra

function GF2_nullspace(H::Matrix{Int64})

    M,N = size(H)

    A = [H;I(N)]

    # Put the upper part (H) in column echelon form by column-wise Gaussian 
    # elimination
    for j=1:M
        if A[j,j] != 1
            p = j+1
            while p <= N && A[j,p] != 1
                p +=1
            end
            if p <= N
                column_j = A[:,j]
                column_p = A[:,p]
                A[:,j] = column_p
                A[:,p] = column_j
            end
        end
        for k=j+1:N
            if A[j,k] == 1
                A[:,k] = (A[:,k] + A[:,j]) .% 2
            end
        end
    end

    A_sup = A[1:M,:]
    A_inf = A[M+1:end,:]

    # find the zero columns of A_sup
    zero_columns = []

    for j = 1:N
        if sum(A_sup[:,j]) == 0
            append!(zero_columns, j)
        end
    end

    # The nullspace of H is the columns of A_inf corresponding to the zero 
    # columns of A_sup

    null_space_H = zeros(Int, N,length(zero_columns))
    for j in eachindex(zero_columns)
        null_space_H[:,j] = A_inf[:,zero_columns[j]]
    end

    return A_sup, A_inf, null_space_H

end

function GF2_nullspace(H::BitMatrix)

    M,N = size(H)

    A = [H;I(N)]

    # Put the upper part (H) in column echelon form by column-wise Gaussian 
    # elimination
    for j=1:M
        if A[j,j] != 1
            p = j+1
            while p <= N && A[j,p] != 1
                p +=1
            end
            if p <= N
                column_j = A[:,j]
                column_p = A[:,p]
                A[:,j] = column_p
                A[:,p] = column_j
            end
        end
        for k=j+1:N
            if A[j,k] == 1
                @. A[:,k] = A[:,k] ⊻ A[:,j]
            end
        end
    end

    A_sup = A[1:M,:]
    A_inf = A[M+1:end,:]

    # find the zero columns of A_sup
    zero_columns = []

    for j = 1:N
        if sum(A_sup[:,j]) == 0
            append!(zero_columns, j)
        end
    end

    # The nullspace of H is the columns of A_inf corresponding to the zero 
    # columns of A_sup

    null_space_H = BitMatrix(undef,N,length(zero_columns))
    null_space_H .= false
    for j in eachindex(zero_columns)
        null_space_H[:,j] = A_inf[:,zero_columns[j]]
    end

    return A_sup, A_inf, null_space_H

end