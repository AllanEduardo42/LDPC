function update_list!(
        inlist::Union{Matrix{Bool},Nothing},
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer},
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

            update_inlist!(inlist,listm,listn,m,n)

            for j=listsize:-1:i+1
                listres[j] = listres[j-1]
                listm[j] = listm[j-1]
                listn[j] = listn[j-1]
            end
        else
            i = 1
        end

        listm[i] = m
        listn[i] = n
        listres[i] = residue        
    end
end

function update_inlist!(
    inlist::Matrix{Bool},
    listm::Vector{<:Integer},
    listn::Vector{<:Integer},
    m::Integer,
    n::Integer
)
    mm = listm[end-1]
    if mm ≠ 0
        nn = listn[end-1]
        inlist[mm,nn] = false
    end
    inlist[m,n] = true

end

function update_inlist!(
    ::Nothing,
    ::Vector{<:Integer},
    ::Vector{<:Integer},
    ::Integer,
    ::Integer
)

end
