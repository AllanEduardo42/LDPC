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

    if nB > 1
        C = similar(A,mA,nB)
    else
        C = similar(A,mA)
    end

    if nA == mB
        if nB > 1
            _gf2_mat_mult!(C,A,B,mA,nA,nB)
        else
            _gf2_mat_mult!(C,A,B,mA,nA)
        end
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
    gf2_mat_mult!(
        C::AbstractMatrix{Bool},
        A::AbstractMatrix{Bool},
        B::AbstractArray{Bool}
    )

    (mC, nC) =  (ndims(C) == 2) ? size(C) : (length(C),1)
    (mA, nA) =  (ndims(A) == 2) ? size(A) : (length(A),1)
    (mB, nB) =  (ndims(B) == 2) ? size(B) : (length(B),1)

    if nA == mB 
        if mC == mA && nC == nB
            if nB > 1
                _gf2_mat_mult!(C,A,B,mA,nA,nB)
            else
                _gf2_mat_mult!(C,A,B,mA,nA)
            end
        else
            throw(
                DimensionMismatch(
                    lazy"C has dimensions ($mC,$nC) but A*B has dimensions ($mA,$nB)"
                )
            )
        end
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
    _gf2_mat_mult!(
        C::AbstractMatrix{Bool},
        A::AbstractMatrix{Bool},
        B::AbstractMatrix{Bool},
        mA::Int,
        nA::Int,
        nB::Int
    )
    
    @inbounds for i in 1:mA
        for j in 1:nB
            result = false
            for k in 1:nA
                if A[i,k]
                    if B[k,j]
                         result ⊻= true
                    end
                end
            end
            C[i,j] = result
        end
    end

end

function 
    _gf2_mat_mult!(
        y::AbstractVector{Bool},
        A::AbstractMatrix{Bool},
        x::AbstractVector{Bool},
        mA::Int,
        nA::Int
    )
    
    @inbounds for i in 1:mA
        result = false
        for k in 1:nA
            if A[i,k]
                if x[k]
                    result ⊻= true
                end
            end
        end
        y[i] = result
    end
end

################################# GF2 NULLSPACE ################################

function gf2_nullspace(A::AbstractMatrix{Bool})

    M,N = size(A)

    AA = [A;I]

    _ = gf2_column_echelon_form!(AA,N)

    @inbounds begin
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

        nullspace_A = similar(A,N,length(zero_columns))
        j = 0
        for column in zero_columns
            j += 1
            nullspace_A[:,j] = view(AA_inf,:,column)
        end
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

function gf2_column_echelon_form!(AA::AbstractMatrix{Bool},N::Int)

    full_rank_sub_matrix = true

    @inbounds for j in 1:N-1
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



function gf2_reduce!(AA::AbstractMatrix{Bool})

    _,N = size(AA)
    gf2_reduce!(AA,N)

end

function gf2_reduce!(AA::AbstractMatrix{Bool},N::Int)
    
    @inbounds for j in N:-1:2
        for k in j-1:-1:1
            if AA[j,k]
                @. AA[:,k] ⊻= AA[:,j]
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

function find_gf2_invertible_matrix(M::Int)

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
        L::Matrix{Bool},
        U::Matrix{Bool},
        y::Vector{Bool}
    )

    M = is_LU_conformable(L,U,y)

    x = zeros(Bool,M)

    _gf2_solve_LU!(x,L,U,y,M)

    return x

end

function 
    gf2_solve_LU!(
        x::Vector{Bool},
        L::Matrix{Bool},
        U::Matrix{Bool},
        y::Vector{Bool}
    )

    M = is_LU_conformable(L,U,y)
    Lx = length(x)

    if Lx ≠ M
        throw(
            DimensionMismatch(
                lazy"x must have the dimension $M, got $Lx"
            )
        )
    end

    _gf2_solve_LU!(x,L,U,y,M)

end

function 
    is_LU_conformable(
        L::Matrix{Bool},
        U::Matrix{Bool},
        y::Vector{Bool}
    )

    ML,NL = size(L)
    MU,NU = size(U)
    L = length(y)

    if ML ≠ NL 
        throw(
            DimensionMismatch(
                lazy"matrix L is not square: dimensions are ($ML,$NL)"
            )
        )
    elseif !istril(L)
        throw(
            ArgumentError(
                lazy"Matrix L must be lower triangular triangular"
            )
        )
    elseif MU ≠ NU
        throw(
            DimensionMismatch(
                lazy"matrix U is not square: dimensions are ($MU,$NU)"
            )
        )
    elseif !istriu(U)
        throw(
            ArgumentError(
                lazy"Matrix U must be upper triangular"
            )
        )
    elseif ML ≠ MU
        throw(
            DimensionMismatch(
                lazy"matrices L and U must have same dimension: L has dimensions ($ML,$NL), U has dimensions ($MU,$NU)"
            )
        )
    elseif L ≠ ML
        throw(
            DimensionMismatch(
                lazy"y must have the dimension $ML, got $L"
            )
        )
    end

    return ML

end


function 
    _gf2_solve_LU!(
        x::Vector{Bool},
        L::Matrix{Bool},
        U::Matrix{Bool},
        y::Vector{Bool},
        M::Int
    )

    @inbounds begin
        _gf2_solve_L!(x,L,y,M)
        _gf2_solve_U!(x,U,M)
    end
end

function 
    _gf2_solve_L!(
        x::Vector{Bool},
        L::Matrix{Bool},
        y::Vector{Bool},
        M::Int
    )

    @inbounds begin
        for i in 1:M
            result = y[i]
            for j=1:i-1
                if L[i,j]
                    if x[j]
                        result ⊻= true
                    end
                end
            end
            x[i] = result
        end
    end
end

function 
    _gf2_solve_U!(
        x::Vector{Bool},
        U::Matrix{Bool},
        M::Int
    )

    @inbounds begin
        for i in M:-1:1
            result = x[i]
            for j=i+1:M
                if U[i,j]
                    if x[j]
                        result ⊻= true
                    end
                end
            end
            x[i] = result
        end
    end
end