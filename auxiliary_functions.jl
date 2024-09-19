################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Auxiliary functions

"""For each column j of matrix H, find the indices i where H[i,j] = 1, and 
return a vector "nodes2checks" where nodes2checks[j] is a vector containing the
indices i where H[i,j] = 1"""
function find_nodes2checks(H::BitMatrix)

    N = size(H,2)
    nodes2checks = Vector{Vector{Int}}()
    for n in 1:N
        push!(nodes2checks,findall(x -> x == true, H[:,n]))
    end

    return nodes2checks

end

"""For each row i of matrix H, find the indices j where H[i,j] = 1, and 
return a vector "checks2nodes" where checks2nodes[i] is a vector containing the
indices j where H[i,j] = 1"""
function find_checks2nodes(H::BitMatrix)

    M = size(H,1)
    checks2nodes = Vector{Vector{Int}}()
    for m in 1:M
        push!(checks2nodes,findall(x -> x == true, H[m,:]))
    end

    return checks2nodes

end

"""Normalize a binary probability distribution function"""
function normalize!(f::Matrix{<:AbstractFloat})

    N = size(f,1)
    
    for n in 1:N
        @inbounds @fastmath α = f[n,1] + f[n,2]
        @inbounds @fastmath f[n,1] = f[n,1]/α
        @inbounds @fastmath f[n,2] = f[n,2]/α
    end

end

"""Return the syndrome H*d, where H is the parity check matrix and d an estimate
of the transmited codeword."""
function 
    calc_syndrome!(
        syndrome::Vector{Bool},
        d::Vector{Bool},
        checks2nodes::Vector{Vector{T}} where {T<:Integer}
    )

    syndrome .*= false
    m = 0
    for indices in checks2nodes
        m += 1
        for n in indices
            @inbounds syndrome[m] ⊻= d[n]
        end
    end

end