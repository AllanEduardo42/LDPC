################################################################################
# Allan Eduardo Feitosa
# 26 mai 2025
# Auxiliary functions

# Functions to create lists of VNs and CNs
"""For each column j of matrix H, find the indices i where H[i,j] = 1, and 
return a vector "vn2cn" where vn2cn[j] is a vector containing the
indices i where H[i,j] = 1"""
function make_vn2cn_list(H::Matrix{Bool})::Vector{Vector{Int}}

    N = size(H,2)
    vn2cn = Vector{Vector{Int}}()
    @inbounds for n in 1:N
        push!(vn2cn,findall(isone, H[:,n]))
    end

    return vn2cn

end

"""For each row i of matrix H, find the indices j where H[i,j] = 1, and 
return a vector "cn2vn" where cn2vn[i] is a vector containing the
indices j where H[i,j] = 1"""
function make_cn2vn_list(H::Matrix{Bool})::Vector{Vector{Int}}

    M = size(H,1)
    cn2vn = Vector{Vector{Int}}()
    @inbounds for m in 1:M
        push!(cn2vn,findall(isone, H[m,:]))
    end

    return cn2vn

end

# Returns true if the algorithm uses Residual Decay (RD) strategy 
function is_RD_algorithm(algorithm::String)
    
    return algorithm == "RD-RBP"   || 
           algorithm == "List-RBP" ||
           algorithm == "C-RBP"    ||
           algorithm == "C&R-RBP"  ||
           algorithm == "C&DR-RBP"

end

# Function to find the girth (shortest cicle) in a bipartite graph using random
# walk

function find_girth(H,max)::Int

    M,N = size(H)

    count = -ones(Int,N)

    visited = zeros(Bool,N)

    girth = N

    node = rand(1:N)

    count[node] = 1

    visited[node] = true

    vn2cn = make_vn2cn_list(H)

    check = rand(vn2cn[node])

    cn2vn = make_cn2vn_list(H)

    for i = 1:max
        if check == -1
            visited .*= false
            node = rand(1:N)
        else
            nodes = filter(x->x!=node,cn2vn[check])
            if !isempty(nodes)
                node = rand(nodes)
            end
        end
        if visited[node]
            girth = min(girth,i - count[node])
        end
        count[node] = i
        visited[node] = true
        checks = filter(x->x!=check,vn2cn[node])
        if length(checks) == 0
            check = -1
        else
            check = rand(checks)
        end
    end

    return 2*girth

end

function generate_parity_matrix(
    H::Matrix{Bool},
    L::Matrix{Bool},
    U::Matrix{Bool}
)

    M,N = size(H)
    K = N - M
    P = zeros(Bool,M,K)
    H1 = H[:,1:K]

    @turbo for k in axes(P,2)
        P[:,k] = gf2_solve_LU(L,U,H1[:,k])
    end

    return P

end