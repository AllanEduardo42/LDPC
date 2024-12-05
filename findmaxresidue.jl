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
        listres1::Vector{<:AbstractFloat},
        listm1::Vector{<:Integer},
        listn1::Vector{<:Integer},
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Integer,
        ::Integer,
        ::Nothing
    )

    @inbounds residues[addressinv[l]] = x

    @fastmath @inbounds if x > listres1[1]
        listres1[1] = x
        listm1[1] = m
        listn1[1] = n
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
        listres1::Vector{<:AbstractFloat},
        listm1::Vector{<:Integer},
        listn1::Vector{<:Integer},
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Integer,
        ::Integer,
        ::Nothing
    )

    @fastmath @inbounds if x > listres1[1]
        listres1[1] = x
        listm1[1] = m
        listn1[1] = n
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
        listres1::Vector{<:AbstractFloat},
        listm1::Vector{<:Integer},
        listn1::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        listsize1::Integer,
        listsize2::Integer,
        inlist::Matrix{<:Integer}
    )

    @inbounds if inlist[l]  # if residue(m,n) is in the list
        inlist[l] = false   # remove from the list
        pos = 0
        for i = 1:listsize1
            if listm1[i] == m
                if listn1[i] == n
                    pos = i
                    break
                end
            end
        end
        if pos == 0
            throw(error("($m,$n) is on the list, but it's not registered."))
        end
        for j in pos:listsize1
            listres1[j] = listres1[j+1]
            listm1[j] = listm1[j+1]
            listn1[j] = listn1[j+1]
        end
    end

    if listsize2 == 0
        update_list!(inlist,listres1,listm1,listn1,x,m,n,listsize1)
    else
        update_list!(nothing,listres2,listm2,listn2,x,m,n,listsize2)
    end
    
end