################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Auxiliary functions

function find_indices_col(H::BitMatrix)

    N = size(H,2)
    indices_col = Vector{Vector{Int64}}(undef, N)
    for n in 1:N
        indices_col[n] = findall(x -> x == 1, H[:,n])
    end

    return indices_col

end

function find_indices_row(H::BitMatrix)

    M = size(H,1)
    indices_row = Vector{Vector{Int64}}(undef, M)
    for m in 1:M
        indices_row[m] = findall(x -> x == 1, H[m,:])
    end

    return indices_row

end

function normalize!(f::Matrix{Float64})

    N = size(f,1)
    
    for n in 1:N
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
    m = 0
    for indices in indices_row
        m += 1
        for n in indices
            @inbounds syndrome[m] ⊻= d[n]
        end
    end

end