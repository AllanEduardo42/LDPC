function add_to_list!(
        inlist::Union{Matrix{Bool},Nothing},
        listres::Vector{<:AbstractFloat},
        coords::Matrix{<:Integer},
        residue::AbstractFloat,
        m::Integer,
        n::Integer,
        listsize::Integer
    )

    @fastmath @inbounds if residue > listres[listsize]

        if listsize > 1
            if residue ≥ listres[1]
                i = 1
            else
                d = listsize >> 1
                i = d
                while d > 1
                    d >>= 1
                    if residue ≥ listres[i]
                        i -= d
                    else
                        i += d
                    end
                end
                if residue < listres[i]
                    i += 1
                end
            end

            update_inlist!(inlist,coords,m,n)

            for j=listsize:-1:i+1
                listres[j] = listres[j-1]
                coords[1,j] = coords[1,j-1]
                coords[2,j] = coords[2,j-1]
            end
        else
            i = 1
        end

        coords[1,i] = m
        coords[2,i] = n
        listres[i] = residue        
    end
end

function update_inlist!(
    inlist::Matrix{Bool},
    coords::Matrix{<:Integer},
    m::Integer,
    n::Integer,
)
    mlast = coords[1,end-1]
    nlast = coords[2,end-1]
    if mlast ≠ 0
        inlist[mlast,nlast] = false
    end
    inlist[m,n] = true

end

function update_inlist!(
    ::Nothing,
    ::Matrix{<:Integer},
    ::Integer,
    ::Integer
)

end
