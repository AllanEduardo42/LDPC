################################################################################
# Allan Eduardo Feitosa
# 3 set 2024
# PEG Algorithm

using Polynomials

function 
    PEG(
        lambda::Vector{Float64},
        rho::Vector{Float64},
        M::Int,
        N::Int
    )::Tuple{Matrix{Bool},Int}

    @inbounds begin
        BitNode_poly = lambda[end:-1:1]
        # ChkNode_poly = rho[end:-1:1]
        Nv = zeros(Int,length(BitNode_poly))
        # Nc = zeros(Int,length(ChkNode_poly))
        d = Vector{Int}()

        # Computing the Integral of lamda and rho polynomial distribution
        Lambda = integrate(Polynomial(BitNode_poly))
        # Rou = integrate(Polynomial(ChkNode_poly))
        temp1 = N/(Lambda(1) - Lambda(0))
        # temp2 = M/(Rou(1) - Rou(0))

        for j in eachindex(BitNode_poly)
            Nv[j] = round(Int,temp1*BitNode_poly[j]/j)
        end
        
        # for j in eachindex(ChkNode_poly)
        #     Nc[j] = round(Int,temp2*ChkNode_poly[j]/j)
        # end
        
        if sum(Nv) ≠ N
            Nv[end] = N - sum(Nv[1:end-1])
        end
        
        # if sum(Nc) ≠ M
        #     Nc[end] = M - sum(Nc[1:end-1])
        # end

        for kk in eachindex(Nv)
            append!(d,kk*ones(Int,Nv[kk]))
        end
        H = zeros(Bool,M,N)
        check_degrees = zeros(Int,M)
        girth = typemax(Int)

        for n in 1:N
            for k in 1:d[n]
                if k == 1
                    _,m = findmin(check_degrees)
                    H[m,n] = true
                    check_degrees[m] += 1
                else
                    L0_checks = findall(isone,H[:,n])
                    level = subtree!(H,M,N,n,0,check_degrees,L0_checks,[n],
                                copy(L0_checks),copy(L0_checks))
                    if level > 0
                        girth = min(girth, 2*(level+1))
                    end
                end
            end
        end
    end

    return H, girth

end

function 
    subtree!(
        H::Matrix{Bool},
        M::Int,
        N::Int,
        root::Int,
        level::Int,
        check_degrees::Vector{Int},
        L0_checks::Vector{Int},
        parent_nodes::Vector{Int},
        L0_check_set::Vector{Int},
        L1_check_set::Vector{Int}
    )    
    # L0: previous level in the subtree
    # L1: current level in the subtree

    @inbounds begin

        STOPPED = false                     # flag: subtree stopped increasing
        level += 1                          # absolute level in the subtree
        L1_checks = Vector{Int}(undef,0)    # checks in the current level
        L1_nodes = Vector{Int}(undef,0)     # nodes in the current level

        row = zeros(Bool,N)
        column = zeros(Bool,M)

        for check in L0_checks
            for i in 1:N
                row[i] = H[check,i]
            end
            row[parent_nodes] .= 0              # remove parent nodes
            check_nodes = findall(row)    # nodes linked to the current check
            append_sort_unique!(L1_nodes,check_nodes)
            for node in check_nodes
                for i in 1:M
                    column[i] = H[i,node]
                end
                column[check] = 0  
                # include checks linked to the current node
                append_sort_unique!(L1_checks,findall(column))
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
            H[idx,root] = true
            check_degrees[idx] += 1
        else
            append_sort_unique!(L0_check_set,L1_checks)
            level = 
                subtree!(H,M,N,root,level,check_degrees,L1_checks,L1_nodes,L0_check_set,
                    L1_check_set)
        end
    end

    return level
end

function 
    append_sort_unique!(
        x::Vector{Int},
        y::Vector{Int}
    )

    append!(x,y)
    sort!(x)
    unique!(x)

end