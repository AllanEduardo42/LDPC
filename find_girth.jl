include("auxiliary_functions.jl")

function find_girth(H,max)

    M,N = size(H)

    count = -ones(Int,N)

    visited = zeros(Bool,N)

    girth = N

    # println("i = 1")

    node = rand(1:N)

    # println("node = $node")

    count[node] = 1

    visited[node] = true

    vn2cn = make_vn2cn_list(H)

    check = rand(vn2cn[node])

    # println("check = $check")

    cn2vn = make_cn2vn_list(H)

    for i = 1:max
        # println()
        # println("i = $(i+1)")
        if check == -1
            visited .*= false
            node = rand(1:N)
        else
            nodes = filter(x->x!=node,cn2vn[check])
            node = rand(nodes)
        end
        # println("node = $node")
        if visited[node]
            girth = min(girth,i - count[node])
            # println("    girth = $(2*girth)")
        end
        count[node] = i
        visited[node] = true
        checks = filter(x->x!=check,vn2cn[node])
        if length(checks) == 0
            check = -1
        else
            check = rand(checks)
        end
        # println("check = $check")
    end

    return 2*girth

end