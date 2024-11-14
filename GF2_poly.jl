################################################################################
# Allan Eduardo Feitosa
# 28 Out 2024
# GF(2) polynomial functions

function
    gf2_poly(p::String)

    L = length(p)
    if L == 1
        if p[1] == '1'
            return [true]
        elseif p[1] == 'x'
            return [true, false]
        else
            return nothing
        end
    end
    first = true
    coeffs = Vector{Bool}()
    N = 0
    l = 1
    while l ≤ L
        if p[l] == 'x'
            if l < L
                l += 1
                if p[l] == '^'
                    l += 1
                    str = ""
                    while l ≤ L && p[l] ≠ ' ' && p[l] ≠ '+'
                        str *= p[l]
                        l += 1
                    end
                    c = parse(Int,str)
                else
                    c = 1
                end                
                if first
                    first = false
                    N = c + 1
                    coeffs = zeros(Bool,c+1)
                end
                coeffs[N - c] = true
            else
                coeffs[N-1] = true
            end
        elseif p[l] == '1'
            coeffs[N] = true
        end
        l += 1
    end
    return coeffs
end

function divide_poly(
    p::Union{Vector{Bool},Vector{Union{Bool,Missing}}},
    d::Vector{Bool}
    )

    if d[1] != 1
        throw(ArgumentError(
            lazy"d[1] must be one."
        ))
    end
    N = length(p)
    # p(x) = p[1]*x^(N-1) + p[2]*x^(N-2) + ... + p[N-2]*x^2 + p[N-1]*x + p[N]
    M = length(d)
    # d(x) = d[1]*x^(M-1) + d[2]*x^(M-2) + ... + d[M-2]*x^2 + d[M-1]*x + d[M]

    # p(x) = q(x)*d(x) + r(x)
    # Degree{q(x)} = Degree{p(x)} - Degree{d(x)} ==>
    # Length{q(x)} - 1 = Length{p(x)} - 1 - (Length{d(x)} - 1) ==>
    # Length{q(x)} = Length{p(x)} - Length{d(x)} + 1 ==>
    # Length{q(x)} = N - M + 1
    L = N-M+1
    if L < 1
        return [0], p
    end
    q = zeros(Bool,L)
    # q(x) = q[1]*x^(L-1) + q[2]*x^(L-2) + ... + q[L-2]*x^2 + q[L-1]*x + q[L]
    a = copy(p)
    for i = 1:L
        if a[i]
            q[i] = true
            for j = 1:M
                # a[i+j-1]*x^(N+1-i-j) + q[i]*x^(L-i) * d[j]*x^(M-j) =
                # a[i+j-1]*x^(N+1-i-j) + (q[i]*d[j]) * x^(L+M-i-j) =
                # (a[i+j-1] + q[i]*d[j])*x^(N+1-i-j), since N+1=L+M
                a[j+i-1] ⊻= d[j]
            end
        end
    end

    # Degree{r(x)} < Degree{d(x)} ==>
    # Length{r(x)} - 1 < Length{d(x)} - 1 ==>
    # Length{r(x)} < M ==>
    # Length{r(x)} ≤ M - 1 ==>
    # r(x) = a[N-(M-1)+1:end]

    return q, a[N-M+2:end]

end

function mul_poly(q::Vector{Bool},d::Vector{Bool})

    L = length(q)
    M = length(d)
    N = L + M - 1

    p = zeros(Bool,N)

    for i=1:L
        for j=1:M
            p[i+j-1] ⊻= q[i] && d[j]
        end
    end

    return p

end

function sum_poly(p::Vector{Bool},r::Vector{Bool})

    N = length(p)
    R = length(r)

    if N < R
        a = [zeros(Bool,R-N);p]
        b = r
    elseif R < N
        a = p
        b = [zeros(Bool,N-R);r]
    else
        a = p
        b = r
    end

    return a .⊻ b

end