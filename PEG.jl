################################################################################
# Allan Eduardo Feitosa
# 3 set 2024
# PEG Algorithm

function 
    PEG(
        d::Vector{<:Integer},
        M::Integer,
        N::Integer
    )

    sort!(d)
    H = zeros(Bool,M,N)
    check_degrees = zeros(Int,M)
    girth = Inf

    @inbounds for n in 1:N
        for k in 1:d[n]
            if k == 1
                _,m = findmin(check_degrees)
                H[m,n] = 1
                check_degrees[m] += 1
            else
                L0_checks = findall(isone,H[:,n])
                level = subtree!(H,M,n,0,check_degrees,L0_checks,[n],
                            copy(L0_checks),copy(L0_checks))
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
        H::Matrix{Bool},
        M::Integer,
        root::Integer,
        level::Integer,
        check_degrees::Vector{<:Integer},
        L0_checks::Vector{<:Integer},
        parent_nodes::Vector{<:Integer},
        L0_check_set::Vector{<:Integer},
        L1_check_set::Vector{<:Integer}
    )    
    # L0: previous level in the subtree
    # L1: current level in the subtree

    @inbounds begin

        STOPPED = false                     # flag: subtree stopped increasing
        level += 1                          # absolute level in the subtree
        L1_checks = Vector{Int}(undef,0)    # checks in the current level
        L1_nodes = Vector{Int}(undef,0)     # nodes in the current level

        for check in L0_checks
            row = H[check,:]
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
            _,m = findmin(check_degrees[compl])
            idx = compl[m]
            H[idx,root] = 1
            check_degrees[idx] += 1
        else
            append_sort_unique!(L0_check_set,L1_checks)
            level = 
                subtree!(H,M,root,level,check_degrees,L1_checks,L1_nodes,L0_check_set,
                    L1_check_set)
        end
    end

    return level
end

function 
    append_sort_unique!(
        x::Vector{<:Integer},
        y::Vector{<:Integer}
    )

    append!(x,y)
    sort!(x)
    unique!(x)

end