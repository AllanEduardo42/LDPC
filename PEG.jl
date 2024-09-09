################################################################################
# Allan Eduardo Feitosa
# 3 set 2024
# PEG Algorithm

using SparseArrays

function 
    PEG!(
        d::Vector{Int64},
        M::Int64
    )

    sort!(d)
    N = length(d)
    H = BitArray(undef,M,N)
    H *= false
    check_degrees = zeros(Int,M)
    girth = 2*N

    for i in 1:N
        for k in 1:d[i]
            if k == 1
                _,j = findmin(check_degrees)
                H[j,i] = 1
                check_degrees[j] += 1
            else
                L0_checks = findall(isone,H[:,i])
                level = 
                    subtree!(
                        H,
                        M,
                        i,
                        0,
                        check_degrees,
                        L0_checks,[i],
                        copy(L0_checks),
                        copy(L0_checks)
                    )
                ;
                if level > 0
                    girth = min(girth, 2*(level+1))
                end
            end
        end
    end

    return H, girth

end

function 
    subtree!(
        H::BitMatrix,
        M::Int64,
        root::Int64,
        level::Int64,
        check_degrees::Vector{Int64},
        L0_checks::Vector{Int64},
        parent_nodes::Vector{Int64},
        L0_check_set::Vector{Int64},
        L1_check_set::Vector{Int64}
    )    
    # L0: previous level in the subtree
    # L1: current level in the subtree

    STOPPED = false                     # flag: subtree stopped increasing
    level += 1                          # absolute level in the subtree
    L1_checks = Vector{Int64}(undef,0)  # checks in the current level
    L1_nodes = Vector{Int64}(undef,0)   # nodes in the current level
    # L1_check_set = copy(L0_check_set)
    for check in L0_checks
        row =  H[check,:]
        row[parent_nodes] .= 0              # remove parent nodes
        check_nodes = findall(isone,row)    # nodes linked to the current check
        append_sort_unique!(L1_nodes,check_nodes)
        for node in check_nodes
            column = H[:,node]
            column[check] = 0  
            # include checks linked to the current node
            append_sort_unique!(L1_checks,findall(isone,column))
        end
    end
    append_sort_unique!(L1_check_set,L1_checks)
    len0 = length(L0_check_set)
    len1 = length(L1_check_set)
    if len1 == len0                     # subtee stopped increasing 
        STOPPED = true
        level = 0
    end
    if STOPPED || len1 == M             # stopped increasing or all checks reached
        compl = collect(1:M)
        deleteat!(compl, L0_check_set)  # take the complement set
        _,j = findmin(check_degrees[compl])
        idx = compl[j]
        H[idx,root] = 1
        check_degrees[idx] += 1
    else
        append_sort_unique!(L0_check_set,L1_checks)
        level = 
            subtree!(
                H,
                M,
                root,
                level,
                check_degrees,
                L1_checks,
                L1_nodes,
                L0_check_set,
                L1_check_set
            )
        ;
    end

    return level
end

function 
    append_sort_unique!(
        x::Vector{Int64},
        y::Vector{Int64}
    )

    append!(x,y)
    sort!(x)
    unique!(x)

end

# function node_degrees(M, N, λ, ρ)

#     # d = zeros(Int, N)
#     Nv = zeros(Int, length(λ))
#     Nc = zeros(Int, length(ρ))
#     Λ = 0
#     for i in eachindex(λ)
#         Λ += λ[i]/i
#     end
#     R = 0
#     for i in eachindex(ρ)
#         R += ρ[i]/i
#     end

#     H = round(Int,(M/R + N/Λ)/2)

#     for i in eachindex(λ)
#         Nv[i] = round(Int, H*λ[i]/i)
#     end

#     for i in eachindex(ρ)
#         Nc[i] = round(Int,H*ρ[i]/i)
#     end

#     return Nv, Nc

# end