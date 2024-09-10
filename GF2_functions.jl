################################################################################
# Allan Eduardo Feitosa
# 9 set 2024
# GF(2) Matrix Functions

using Random
using LinearAlgebra

function
    GF2_mat_mult(
        A::BitMatrix,
        B::BitMatrix
    )::BitMatrix

    A_bool = Matrix(A)
    B_bool = Matrix(B)

    mA, nA = size(A)
    mB, nB = size(B)

    if nA == mB
        C = zeros(Bool,mA,nB)
        for i in 1:mA
            for j in 1:nB
                for k in 1:nA
                    @inbounds C[i,j] ⊻= A_bool[i,k] && B_bool[k,j]
                end
            end
        end
    else
        throw(
            DimensionMismatch(
                lazy"A has dimensions ($mA,$nA) but B has dimensions ($mB,$nB)"
            )
        )
    end

    return BitMatrix(C)

end

function 
    GF2_mat_mult(
        A::BitMatrix,
        x::Vector{Bool}
    )::Vector{Bool}

    A_bool = Matrix(A)

    mA, nA = size(A)
    L = length(x)
    if nA == L
        y = zeros(Bool,mA)
        for i in 1:mA
            for k in 1:nA
                @inbounds y[i] ⊻= A_bool[i,k] && x[k]
            end
        end
    else
        throw(
            DimensionMismatch(
                lazy"second dimension of A, $nA, does not match length of x, $L"
            )
        )
    end

    return y

end

function GF2_nullspace(A::BitMatrix)

    M,N = size(A)

    AA = [A;I]

    invertible = GF2_column_echelon_form!(AA,N)

    if invertible
        return falses(0,0)
    else
        AA_sup = view(AA,1:M,:)
        AA_inf = view(AA,M+1:M+N,:)

        # find the zero columns of A_sup
        zero_columns = []

        for j = 1:N
            if iszero(view(AA_sup,:,j))
                append!(zero_columns, j)
            end
        end

        # The nullspace of A is the columns of A_inf corresponding to the zero 
        # columns of A_sup

        null_space_A = falses(N,length(zero_columns))
        j = 0
        for column in zero_columns
            j += 1
            null_space_A[:,j] = view(AA_inf,:,column)
        end

        return null_space_A
    end

end

function GF2_inverse(A::BitMatrix;IAEF=false)

    # IAEF: invertible augmented echelon form

    M,N = size(A)

    if !IAEF

        if M != N
            throw(
                DimensionMismatch(
                    lazy"matrix is not square: dimensions are ($M,N)"
                )
            )
        end

        AA = [A; I]

        invertible = GF2_column_echelon_form!(AA,N)

    else
        MM,N = size(A)
        # confirm that is indeed an IAEF
        if N != MM÷2
            throw(
                ArgumentError(
                    "A is not in an invertible augmented echelon form"
                )
            )
        elseif !(istril(A[1:N,:] - I,-1))
            throw(
                ArgumentError(
                    "A is not in an invertible augmented echelon form"
                )
            )
        else
            invertible = true
            AA = copy(A)
        end
    end

    if !invertible
        throw(
            SingularException(1)
        )
    end

    GF2_reduce!(AA,N)

    return AA[N+1:end,:]

end

function GF2_column_echelon_form!(A::BitMatrix,N::Int64)

    full_rank_sub_matrix = true

    for j=1:N
        if A[j,j] != 1
            p = j+1
            while p <= N && A[j,p] != 1
                p +=1
            end
            if p <= N
                @. A[:,j] = A[:,j] ⊻ A[:,p]
                @. A[:,p] = A[:,p] ⊻ A[:,j]
            else
                full_rank_sub_matrix = false
            end
        end
        for i=j+1:N
            if A[j,i] == 1
                @. A[:,i] = A[:,i] ⊻ A[:,j]
            end
        end
    end

    return full_rank_sub_matrix

end

function GF2_reduce!(A::BitMatrix,M::Int64)

    for i=M:-1:1
        for j=i-1:-1:1
            if A[i,j] == 1
                @. A[:,j] = A[:,j] ⊻ A[:,i]
            end
        end
    end

end

function isGF2invertible(A::BitMatrix)

    M,N = size(A)
    AA = [A; I]

    if M != N
        invertible = false
    else
        invertible = GF2_column_echelon_form!(AA,N)
    end

    return invertible, AA

end

function find_GF2_invertible_matrix(M::Int64)

    A = bitrand(M,M)

    invertible, AA = isGF2invertible(A)
    while !invertible
        A = bitrand(M,M)
        invertible, AA = isGF2invertible(A)
    end

    A_inv = GF2_inverse(AA;IAEF=true)

    return A, A_inv

end

