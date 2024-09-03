# PEG Algorithm

function PEG(M, N, d)

    E = zeros(Int,M,N)
    check_degree = zeros(Int,M)
    girth = 0

    for i in 1:N
        for k in 1:d[i]
            if k == 1
                _,j = findmin(check_degree)
                E[j,i] = 1
                check_degree[j] += 1
            else
                checks_L0 = findall(isone,E[:,i])
                count = 0
                E, check_degree, girth = subgraph(E, M, N, check_degree, checks_L0, i, count, girth)
            end
        end
    end

    return E, girth

end

function PEG(M, N, λ, ρ)

    Nv, Nc = node_degrees(M, N, λ, ρ)

    E = zeros(Int,M,N)
    check_degree = zeros(Int,M)
    girth = 0

    for i in 1:N
        degree = Nv
        for k in 1:degree
            if k == 1
                _,j = findmin(check_degree)
                E[j,i] = 1
                check_degree[j] += 1
            else
                checks_L0 = findall(isone,E[:,i])
                count = 0
                E, check_degree, girth = subgraph(E, M, N, check_degree, checks_L0, i, count, girth)
            end
        end
    end

    return E, girth

end

function subgraph(E, M, N, count_degree, checks_L0, root, count, girth; parents=[root], all_checks=copy(checks_L0))    
    
    count += 1
    checks_L1 = Int.([])
    L0 = length(all_checks)
    L1 = 0
    L1_nodes = Int.([])
    for check in checks_L0
        # remove parent nodes
        row =  E[check,:]   
        for parent in parents
            row[parent] = 0
        end
        check_nodes = findall(isone,row)
        union!(L1_nodes, check_nodes)
        for node in check_nodes
            column = E[:,node]
            column[check] = 0
            node_checks = findall(isone,column)
            union!(checks_L1, node_checks)
        end
    end
    union!(all_checks,checks_L1)
    L1 = length(all_checks)
    if L1 == M  # all checks reached
        girth = 2*(count+1)
        _,j = findmin(count_degree[checks_L1])
        E[checks_L1[j],root] = 1
        count_degree[checks_L1[j]] += 1
    elseif L1 == L0 # stopped increasing
        compl = collect(1:M)
        sort!(all_checks)
        deleteat!(compl, all_checks)        
        _,j = findmin(count_degree[compl])
        E[compl[j],root] = 1
        count_degree[compl[j]] += 1
    else
        E, count_degree, girth = subgraph(E, M, N,count_degree, checks_L1, root, count, girth; parents = L1_nodes, all_checks=all_checks)
    end

    return E, count_degree, girth
end

function node_degrees(M, N, λ, ρ)

    # d = zeros(Int, N)
    Nv = zeros(Int, length(λ))
    Nc = zeros(Int, length(ρ))
    Λ = 0
    for i in eachindex(λ)
        Λ += λ[i]/i
    end
    R = 0
    for i in eachindex(ρ)
        R += ρ[i]/i
    end

    E = round(Int,(M/R + N/Λ)/2)

    for i in eachindex(λ)
        Nv[i] = round(Int, E*λ[i]/i)
    end

    for i in eachindex(ρ)
        Nc[i] = round(Int,E*ρ[i]/i)
    end

    return Nv, Nc

end