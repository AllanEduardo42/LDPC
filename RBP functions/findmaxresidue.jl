################################################################################
# Allan Eduardo Feitosa
# 05 dez 2024
# Function to find the maximum residue

# RBP
function 
    findmaxresidue!(
        addressinv::Matrix{<:Integer},
        residues::Vector{<:AbstractFloat},
        m::Integer,
        n::Integer,
        l::Integer,
        x::AbstractFloat,
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer},
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Integer,
        ::Integer,
        ::Nothing
    )

    @inbounds residues[addressinv[l]] = x

    @fastmath @inbounds if x > listres[1]
        listres[1] = x
        listm[1] = m
        listn[1] = n
    end

end

# Local-RBP
function 
    findmaxresidue!(
        ::Nothing,
        ::Nothing,
        m::Integer,
        n::Integer,
        ::Integer,
        x::AbstractFloat,
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer},
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Integer,
        ::Integer,
        ::Nothing
    )

    @fastmath @inbounds if x > listres[1]
        listres[1] = x
        listm[1] = m
        listn[1] = n
    end

end

# List-RBP
function 
    findmaxresidue!(
        ::Nothing,
        ::Nothing,
        m::Integer,
        n::Integer,
        l::Integer,
        x::AbstractFloat,
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        listsize::Integer,
        listsize2::Integer,
        inlist::Matrix{<:Integer}
    )

    @inbounds if inlist[l]  # if residue(m,n) is in the list
        inlist[l] = false   # remove from the list
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

    if listsize2 == 0
        update_list!(inlist,listres,listm,listn,x,m,n,listsize)
    else
        update_list!(nothing,listres2,listm2,listn2,x,m,n,listsize2)
    end
    
end