################################################################################
# Allan Eduardo Feitosa
# 9 set 2024
# GF(2) Matrix Functions

using Random
using LinearAlgebra

########################### GF2 matrix multiplication ##########################

# C = A*B BitArray

function
    gf2_mat_mult(
        A::BitMatrix,
        B::BitMatrix
    )::BitMatrix

    mA, nA = size(A)
    mB, nB = size(B)

    if nA == mB
        C = gf2_mat_mult_core(Matrix(A),Matrix(B),mA,nA,nB)
    else
        throw(
            DimensionMismatch(
                lazy"A has dimensions ($mA,$nA) but B has dimensions ($mB,$nB)"
            )
        )
    end

    return BitMatrix(C)

end

# C = A*B Matrix{Bool}

function
    gf2_mat_mult(
        A::Matrix{Bool},
        B::Matrix{Bool}
    )::Matrix{Bool}

    mA, nA = size(A)
    mB, nB = size(B)

    if nA == mB
        C = gf2_mat_mult_core(A,B,mA,nA,nB)
    else
        throw(
            DimensionMismatch(
                lazy"A has dimensions ($mA,$nA) but B has dimensions ($mB,$nB)"
            )
        )
    end

    return C

end

# y = A*x BitArray

function 
    gf2_mat_mult(
        A::BitMatrix,
        x::BitVector
    )::BitVector

    mA, nA = size(A)
    if nA == length(x)
        y = gf2_mat_mult_core(Matrix(A),Vector(x),mA,nA)
    else
        throw(
            DimensionMismatch(
                lazy"second dimension of A, $nA, does not match length of x, $L"
            )
        )
    end

    return BitVector(y)

end

# y = A*x Matrix{Bool}

function 
    gf2_mat_mult(
        A::Matrix{Bool},
        x::Vector{Bool}
    )::Vector{Bool}

    mA, nA = size(A)
    if nA == length(x)
        y = gf2_mat_mult_core(A,x,mA,nA)
    else
        throw(
            DimensionMismatch(
                lazy"second dimension of A, $nA, does not match length of x, $L"
            )
        )
    end

    return y

end

# Core computation of C = A*B

function 
    gf2_mat_mult_core(
        A::Matrix{Bool},
        B::Matrix{Bool},
        mA::Integer,
        nA::Integer,
        nB::Integer
    )::Matrix{Bool}

    C = zeros(Bool,mA,nB)
    for i in 1:mA
        for j in 1:nB
            for k in 1:nA
                @inbounds C[i,j] ⊻= A[i,k] && B[k,j]
            end
        end
    end

    return C
end

# Core computation of y = A*x

function 
    gf2_mat_mult_core(
        A::Matrix{Bool},
        x::Vector{Bool},
        mA::Integer,
        nA::Integer,
    )::Vector{Bool}

    y = zeros(Bool,mA)
    for i in 1:mA
        for k in 1:nA
            @inbounds y[i] ⊻= A[i,k] && x[k]
        end
    end

    return y
end

################################# GF2 NULLSPACE ################################

function gf2_nullspace(A::BitMatrix)

    M,N = size(A)

    AA = [A;I]

    _ = gf2_column_echelon_form!(AA,N)


    AA_sup = view(AA,1:M,:)
    AA_inf = view(AA,M+1:M+N,:)

    # find the zero columns of AA_sup
    zero_columns = []

    for j = 1:N
        if iszero(view(AA_sup,:,j))
            append!(zero_columns, j)
        end
    end

    # The nullspace of A is the columns of AA_inf corresponding to the zero 
    # columns of AA_sup

    nullspace_A = falses(N,length(zero_columns))
    j = 0
    for column in zero_columns
        j += 1
        nullspace_A[:,j] = view(AA_inf,:,column)
    end

    return nullspace_A

end

############################# GF2 MATRIX INVERSION #############################


function gf2_inverse(A::BitMatrix;ACCEF=false)

    # ACCEF: augmented complete column echelon form

    # ⌈   1      0      0    ...    0      0      0  ⌉
    # |  x₂₁     1      0    ...    0      0      0  |
    # |  x₃₁    x₃₂     1    ...    0      0      0  |
    # |   ⋮       ⋮       ⋮    ⋱      ⋮      ⋮       ⋮  |
    # | xₙ₋₂₁  xₙ₋₂₂  xₙ₋₂₃  ...    1      0      0  |
    # | xₙ₋₁₁  xₙ₋₁₂  xₙ₋₁₃  ...  xₙₘ₋₃    1      0  |
    # |  xₙ₁    xₙ₂    xₙ₃   ...  xₙₘ₋₂  xₙₘ₋₁    1  |   ⌈ X ⌉
    # |----------------------------------------------| = |---|
    # |  y₁₁    y₁₂    y₁₃   ...  y₁ₘ₋₂  y₁ₘ₋₁   y₁ₘ |   ⌊ Y ⌋
    # |   ⋮       ⋮       ⋮    ⋱      ⋮      ⋮       ⋮  |
    # ⌊  yₘ₁    yₘ₂    yₘ₃   ...  yₘₘ₋₂  yₘₘ₋₁   yₘₘ ⌋

    #       ⌈ X ⌉
    # where |---| is the result of Gaussian elimination by applying a sequence
    #       ⌊ Y ⌋                           ⌈ Z ⌉
    # of elementary columns operations over |---|.
    #                                       ⌊ I ⌋ 

    M,N = size(A)

    if !(ACCEF)

        if M ≠ N
            throw(
                DimensionMismatch(
                    lazy"matrix is not square: dimensions are ($M,N)"
                )
            )
        end

        AA = [A; I]

        invertible = gf2_column_echelon_form!(AA,N)

    else
        if M ≠ 2*N
            throw(
                ArgumentError(
                    lazy"matrix dimensions are incompatible ($M ≠ $(2*N))"
                )
            )
        elseif !(istril(A[1:N,:] - I,-1))
            throw(
                ArgumentError(
                    lazy"matrix X in [X;Y] is not lower triangular with diag(X) = I"
                )
            )
        end
        println("""WARNING: Be sure A is the result of a column-wise Gaussian 
                elimination over an augmented matrix [B; I]. Otherwise result
                will almost certainly be wrong.""")
        invertible = true
        AA = A
    end

    if !(invertible)
        throw(
            SingularException(1)
        )
    end

    gf2_reduce!(AA,N)

    return AA[N+1:end,:]

end

function gf2_column_echelon_form!(AA::BitMatrix,N::Integer)

    full_rank_sub_matrix = true

    for j in 1:N-1
        if !(AA[j,j])
            p = j+1
            while p <= N && !(AA[j,p])
                p +=1
            end
            if p <= N
                @. AA[:,j] ⊻= AA[:,p]
                @. AA[:,p] ⊻= AA[:,j]
            else
                full_rank_sub_matrix = false
            end
        end
        for k in j+1:N
            if AA[j,k]
                @. AA[:,k] ⊻= AA[:,j]
            end
        end
    end
    if !(AA[N,N])
        full_rank_sub_matrix = false
    end

    return full_rank_sub_matrix

end

function gf2_reduce!(AA::BitMatrix,N::Integer)

    for j in N:-1:2
        for k in j-1:-1:1
            if AA[j,k]
                @. AA[:,k] ⊻= AA[:,j]
            end
        end
    end

end

function isgf2invertible(A::BitMatrix)

    M,N = size(A)
    AA = [A; I]

    if M ≠ N
        invertible = false
    else
        invertible = gf2_column_echelon_form!(AA,N)
    end

    return invertible, AA

end

function find_gf2_invertible_matrix(M::Integer)

    A = bitrand(M,M)

    invertible, AA = isgf2invertible(A)
    while !invertible
        A = bitrand(M,M)
        invertible, AA = isgf2invertible(A)
    end

    A_inv = gf2_inverse(AA;IAEF=true)

    return A, A_inv

end