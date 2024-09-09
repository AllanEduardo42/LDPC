################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Auxiliary functions

function find_indices_col(H::BitMatrix)

    N = size(H,2)
    indices_col = Vector{Vector{Int64}}(undef, N)
    for n=1:N
        indices_col[n] = findall(x -> x == 1, H[:,n])
    end

    return indices_col

end

function find_indices_row(H::BitMatrix)

    M = size(H,1)
    indices_row = Vector{Vector{Int64}}(undef, M)
    for m=1:M
        indices_row[m] = findall(x -> x == 1, H[m,:])
    end

    return indices_row

end

function normalize!(f::Matrix{Float64})

    N = size(f,1)
    
    for n=1:N
        @inbounds @fastmath α = f[n,1] + f[n,2]
        @inbounds @fastmath f[n,1] = f[n,1]/α
        @inbounds @fastmath f[n,2] = f[n,2]/α
    end

end

function 
    calc_syndrome!(
        syndrome::Vector{Bool},
        d::Vector{Bool},
        indices_row::Vector{Vector{Int64}}
    )

    syndrome .*= false
    M = length(indices_row)
    for m=1:M
        for n in indices_row[m]
            syndrome[m] ⊻= d[n]
        end
    end

end

function ϕ(x::Float64)::Float64
    @fastmath m = exp(x)-1
    @fastmath log(1 + 2/m)
end

function
    bitwise_mat_mult(
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
    bitwise_mat_mult(
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