function update_list!(
        inlist::Union{Matrix{Bool},Nothing},
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer},
        x::AbstractFloat,
        m::Integer,
        n::Integer,
        listsize::Integer
    )

    @fastmath @inbounds if x > listres[listsize]
        if x ≥ listres[1]
            i = 1
        else
            d = listsize >> 1
            i = d
            while d > 1
                d >>= 1
                if x ≥ listres[i]
                    i -= d
                else
                    i += d
                end
            end
            if x < listres[i]
                i += 1
            end
        end
        if inlist !== nothing
            mm = listm[end-1]
            if mm ≠ 0
                nn = listn[end-1]
                inlist[mm,nn] = false
            end
        end
        for j=listsize:-1:i+1
            listres[j] = listres[j-1]
            listm[j] = listm[j-1]
            listn[j] = listn[j-1]
        end
        listm[i] = m
        listn[i] = n
        listres[i] = x
        if inlist !== nothing
            inlist[m,n] = true
        end
    end
end

