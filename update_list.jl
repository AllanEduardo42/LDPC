function update_list!(
        inlist::Matrix{Bool},
        listres::Vector{<:AbstractFloat},
        listadd::Matrix{<:Integer},
        listaddinv::Matrix{<:Integer},
        x::AbstractFloat,
        m::Integer,
        n::Integer,
        listsize::Integer
    )
    @fastmath @inbounds if x > listres[listsize]
        if x > listres[1]
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
        mm = listadd[1,end-1]
        nn = listadd[2,end-1]
        if mm ≠ 0
            inlist[mm,nn] = false
            listaddinv[mm,nn] = 0
        end
        for j=listsize:-1:i+1
            listres[j] = listres[j-1]
            mm = listadd[1,j-1]
            nn = listadd[2,j-1]
            listadd[1,j] = mm
            listadd[2,j] = nn
            if listadd[1,j] != 0
                listaddinv[mm,nn] = j
            end
        end
        listadd[1,i] = m
        listadd[2,i] = n
        listaddinv[m,n] = i
        listres[i] = x
        inlist[m,n] = true
    end

end