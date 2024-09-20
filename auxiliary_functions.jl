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

function 
    llr_init_q!(
        Lq::Matrix{<:AbstractFloat},
        ΔLf::Vector{<:AbstractFloat},
        nodes2checks::Vector{Vector{T}} where {T<:Integer}
    )
    node = 0
    for checks in nodes2checks
        node += 1
        for check in checks
            @inbounds Lq[check,node] = ΔLf[node]
        end
    end
end