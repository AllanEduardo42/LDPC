################################################################################
# Allan Eduardo Feitosa
# 05 dez 2024
# Function to find the maximum residue

# RBP
function 
    update_residue!(
        addressinv::Matrix{<:Integer},
        residues::Vector{<:AbstractFloat},
        ::Integer,
        ::Integer,
        l::Integer,
        x::AbstractFloat,
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        ::Vector{<:AbstractFloat},
        ::Vector{<:Integer},
        ::Vector{<:Integer},
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Integer,
        ::Integer,
        count_size::Integer,
        ::Nothing
    )

    @fastmath @inbounds x *= Factors[l]
    @inbounds residues[addressinv[l]] = x

    return count_size

end

# Local-RBP
function 
    update_residue!(
        ::Nothing,
        ::Nothing,
        m::Integer,
        n::Integer,
        ::Integer,
        x::AbstractFloat,
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer},
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Integer,
        ::Integer,
        count_size::Integer,
        ::Nothing
    )

    @fastmath @inbounds x *= Factors[l]
    @fastmath @inbounds if x > listres[1]
        listres[1] = x
        listm[1] = m
        listn[1] = n
    end

    return count_size

end

# List-RBP
function 
    update_residue!(
        ::Nothing,
        ::Nothing,
        m::Integer,
        n::Integer,
        l::Integer,
        x::AbstractFloat,
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        listsize::Integer,
        listsize2::Integer,
        count_size::Integer,
        inlist::Matrix{<:Integer}
    )

    @fastmath @inbounds if inlist[l]  # if residue(m,n) is in the list
        inlist[l] = false   # remove from the list
        count_size += 1
        pos = 0
        for i = 1:listsize
            if listm[i] == m
                if listn[i] == n
                    pos = i
                    break
                end
            end
        end
        if pos == 0
            throw(error("($m,$n) is on the list, but it's not registered."))
        end
        for j in pos:listsize
            listres[j] = listres[j+1]
            listm[j] = listm[j+1]
            listn[j] = listn[j+1]
        end
    end

    @fastmath @inbounds if x != 0.0
        x *= Factors[l]
        if listsize2 == 0
            update_list!(inlist,listres,listm,listn,x,m,n,listsize)
        elseif count_size == 1
            if x > listres2[1]
                listres2[1] = x
                listm2[1] = m
                listn2[1] = n
            end
        else
            update_list!(nothing,listres2,listm2,listn2,x,m,n,count_size)
        end
    end

    return count_size
    
end