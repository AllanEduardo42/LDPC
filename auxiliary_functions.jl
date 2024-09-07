################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Auxiliary functions

function findindices_M(H::BitMatrix)

    N = size(H,2)
    indices_m = Vector{Vector{Int64}}(undef, N)
    for n=1:N
        indices_m[n] = findall(x -> x == 1, H[:,n])
    end

    return indices_m

end

function findindices_N(H::BitMatrix)

    M = size(H,1)
    indices_n = Vector{Vector{Int64}}(undef, M)
    for m=1:M
        indices_n[m] = findall(x -> x == 1, H[m,:])
    end

    return indices_n

end

function normalize(f, N)

    for n=1:N
        α = f[n,1] + f[n,2]
        f[n,1] = f[n,1]/α
        f[n,2] = f[n,2]/α
    end

end

function calc_syndrome!(syndrome::Vector{Int64},
                        indices_n::Vector{Vector{Int64}},
                        d::Vector{Int64})

    syndrome .*= 0
    M = length(indices_n)
    for m=1:M
        for n in indices_n[m]
            syndrome[m] += d[n]
        end
    end

    syndrome .%= 2

end

function ϕ(x::Float64)::Float64
    @fastmath m = exp(x)-1
    @fastmath log(1 + 2/m)
end

function get_index(arg::Float64)::Int64
    z = unsafe_trunc(Int,arg)
    if z >= SIZE
        i = SIZE
    else
        i = z + 1
    end
    
    return i
end

function bitwise_mat_mult(A::BitMatrix,
                          B::BitMatrix)::BitMatrix

    Ax, Ay = size(A)
    Bx, By = size(B)

    if Ay == Bx
        C = falses(Ax,By)
        for i in 1:Ax
            for j in 1:By
                for k in 1:Ay
                    @inbounds C[i,j] ⊻= A[i,k] && B[k,j]
                end
            end
        end
    end

    return C

end