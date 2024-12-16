################################################################################
# Allan Eduardo Feitosa
# 16 dez 2024
# Function to find the girth (shortest cicle) in a bipartite graph using random
# walk

include("auxiliary_functions.jl")

function find_girth(H,max)

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
            node = rand(nodes)
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