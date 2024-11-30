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
        listadd1::Matrix{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listadd2::Union{Matrix{<:Integer},Nothing},
        listsize1::Integer,
        listsize2::Integer,
        inlist::Matrix{<:Integer}
    )

    @inbounds if inlist[m,n] # if residue(m,n) is in the list
        inlist[m,n] = false   # remove from the list
        pos = 0
        for i = 1:listsize1
            if listadd1[1,i] == m
                if listadd1[2,i] == n
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
            mm = listadd1[1,j+1]
            nn = listadd1[2,j+1]
            listadd1[1,j] = mm
            listadd1[2,j] = nn
        end
    end

    if listsize2 == 0
        update_list!(inlist,listres1,listadd1,x,m,n,listsize1)
        @inbounds maxcoords[1], maxcoords[2] = listadd1[1], listadd1[2]
    else
        update_list!(nothing,listres2,listadd2,x,m,n,listsize2)
    end

    @inbounds return listres1[1]
    
end