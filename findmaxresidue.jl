# RBP and Random-RBP
function 
    findmaxresidue!(
        Residues::Matrix{<:AbstractFloat},
        maxcoords::Vector{<:Integer},
        maxresidue::AbstractFloat,
        m::Integer,
        n::Integer,
        x::AbstractFloat,
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Integer,
        ::Integer,
        ::Nothing
    )

    @inbounds Residues[m,n] = x

    if @fastmath x > maxresidue
        maxresidue = x
        @inbounds maxcoords[1] = m
        @inbounds maxcoords[2] = n
    end

    return maxresidue
end

# Local-RBP
function 
    findmaxresidue!(
        ::Nothing,
        maxcoords::Vector{<:Integer},
        maxresidue::AbstractFloat,
        m::Integer,
        n::Integer,
        x::AbstractFloat,
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Integer,
        ::Integer,
        ::Nothing
    )
    if @fastmath x > maxresidue
        maxresidue = x
        @inbounds maxcoords[1] = m
        @inbounds maxcoords[2] = n
    end

    return maxresidue
end

# List-RBP
function 
    findmaxresidue!(
        ::Nothing,
        maxcoords::Vector{<:Integer},
        maxresidue::AbstractFloat,
        m::Integer,
        n::Integer,
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

    @inbounds if inlist[m,n] # if residue(m,n) is in the list
        inlist[m,n] = false   # remove from the list
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
        @inbounds maxcoords[1], maxcoords[2] = listm1[1], listn1[1]
    else
        update_list!(nothing,listres2,listm2,listn2,x,m,n,listsize2)
    end

    @inbounds return listres1[1]
    
end