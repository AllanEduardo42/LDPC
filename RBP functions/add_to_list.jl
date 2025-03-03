function add_to_list!(
        inlist::Union{Matrix{Bool},Nothing},
        listres::Vector{<:AbstractFloat},
        indices_res::Vector{<:Integer},
        residue::AbstractFloat,
        li::Integer,
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

            update_inlist!(inlist,indices_res,li)

            for j=listsize:-1:i+1
                listres[j] = listres[j-1]
                indices_res[j] = indices_res[j-1]
            end
        else
            i = 1
        end

        indices_res[i] = li
        listres[i] = residue        
    end
end

function update_inlist!(
    inlist::Matrix{Bool},
    indices_res::Vector{<:Integer},
    li::Integer
)
    last = indices_res[end-1]
    if last ≠ 0
        inlist[last] = false
    end
    inlist[li] = true

end

function update_inlist!(
    ::Nothing,
    ::Vector{<:Integer},
    ::Integer
)

end
