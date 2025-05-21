################################################################################
# Allan Eduardo Feitosa
# 9 set 2024
# GF(2) Matrix Functions

using Random
using LinearAlgebra

import Base.*

########################### GF2 matrix multiplication ##########################
*(A::AbstractMatrix{Bool},B::AbstractMatrix{Bool}) = gf2_mat_mult(A,B)
*(A::AbstractMatrix{Bool},B::AbstractVector{Bool}) = gf2_mat_mult(A,B)

function
    gf2_mat_mult(
        A::AbstractMatrix{Bool},
        B::AbstractArray{Bool}
    )

    (mA, nA) =  (ndims(A) == 2) ? size(A) : (length(A),1)
    (mB, nB) =  (ndims(B) == 2) ? size(B) : (length(B),1)

    if nA == mB
        C = _gf2_mat_mult(A,B,mA,nA,nB)
    else
        throw(
            DimensionMismatch(
                lazy"A has dimensions ($mA,$nA) but B has dimensions ($mB,$nB)"
            )
        )
    end

    return C

end

function 
    _gf2_mat_mult(
        A::AbstractMatrix{Bool},
        B::AbstractMatrix{Bool},
        mA::Integer,
        nA::Integer,
        nB::Integer
    )
    
    C = similar(A,mA,nB)
    C .*= false
    for i in 1:mA
        for j in 1:nB
            for k in 1:nA
                @inbounds C[i,j] ⊻= A[i,k] && B[k,j]
            end
        end
    end

    return C

end

function 
    _gf2_mat_mult(
        A::AbstractMatrix{Bool},
        B::AbstractVector{Bool},
        mA::Integer,
        nA::Integer,
        nB::Integer
    )
    
    v = similar(A,mA)
    v .*= false
    for i in 1:mA
        for k in 1:nA
            @inbounds v[i] ⊻= A[i,k] && B[k]
        end
    end

    return v

end

################################# GF2 NULLSPACE ################################

function gf2_nullspace(A::AbstractMatrix{Bool})

    M,N = size(A)

    AA = [A;I]

    _ = gf2_column_echelon_form!(AA,N)


    @inbounds AA_sup = view(AA,1:M,:)
    @inbounds AA_inf = view(AA,M+1:M+N,:)

    # find the zero columns of AA_sup
    zero_columns = []

    for j = 1:N
        if @inbounds iszero(view(AA_sup,:,j))
            append!(zero_columns, j)
        end
    end

    # The nullspace of A is the columns of AA_inf corresponding to the zero 
    # columns of AA_sup

    nullspace_A = similar(A,N,length(zero_columns))
    j = 0
    for column in zero_columns
        j += 1
        @inbounds nullspace_A[:,j] = view(AA_inf,:,column)
    end

    return nullspace_A

end

############################# GF2 MATRIX INVERSION #############################


function gf2_inverse(A::AbstractMatrix{Bool};ACCEF=false)

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

    return @inbounds AA[N+1:end,:]

end

function gf2_row_echelon_form!(AA::AbstractMatrix{Bool})

    M,N = size(AA)
    AAt = zeros(Bool,N,M)
    AAt .= AA'
    gf2_column_echelon_form!(AAt,M)
    AA .= AAt'

end

function gf2_column_echelon_form!(AA::AbstractMatrix{Bool})

    _,N = size(AA)
    return gf2_column_echelon_form!(AA,N)

end

function gf2_column_echelon_form!(AA::AbstractMatrix{Bool},N::Integer)

    full_rank_sub_matrix = true

    for j in 1:N-1
        if !(@inbounds AA[j,j])
            p = j+1
            while p <= N && !(@inbounds AA[j,p])
                p +=1
            end
            if p <= N
                @. @inbounds AA[:,j] ⊻= AA[:,p]
                @. @inbounds AA[:,p] ⊻= AA[:,j]
            else
                full_rank_sub_matrix = false
            end
        end
        for k in j+1:N
            if @inbounds AA[j,k]
                @. @inbounds AA[:,k] ⊻= AA[:,j]
            end
        end
    end
    if !(@inbounds AA[N,N])
        full_rank_sub_matrix = false
    end

    return full_rank_sub_matrix

end



function gf2_reduce!(AA::AbstractMatrix{Bool})

    _,N = size(AA)
    gf2_reduce!(AA,N)

end

function gf2_reduce!(AA::AbstractMatrix{Bool},N::Integer)
    
    for j in N:-1:2
        for k in j-1:-1:1
            if @inbounds AA[j,k]
                @. @inbounds AA[:,k] ⊻= AA[:,j]
            end
        end
    end
end

function isgf2invertible(A::AbstractMatrix{Bool})

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

    A_inv = gf2_inverse(AA;ACCEF=true)

    return A, A_inv

end

function 
    gf2_solve_LU(
        A::Matrix{Bool},
        y::Vector{Bool}
    )

    M,N = size(A)
    if M ≠ N
        throw(
            DimensionMismatch(
                lazy"matrix is not square: dimensions are ($M,N)"
            )
        )
    end

    @inbounds begin
        if istril(A)
            x = zeros(Bool,M)
            for i in 1:M
                x[i] = y[i]
                for j=1:i-1
                    x[i] ⊻= A[i,j] && x[j]
                end
            end
        elseif istriu(A)
            x = zeros(Bool,M)
            for i in M:-1:1
                x[i] = y[i]
                for j=i+1:M
                    x[i] ⊻= A[i,j] && x[j]
                end
            end
        else
            throw(
                    ArgumentError(
                        lazy"Matrix must be triangular"
                    )
                )
        end
    end
    
    return x
end