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
        ::Nothing,
        ::Nothing,
        listaddinv1::Union{Matrix{<:Integer},Nothing},
        listsize1::Integer,
        listsize2::Integer,
        inlist1::Matrix{<:Integer}
    )

    # println("m = $m, n = $n, x = $x, inlist = $(inlist1[m,n])")

    if @inbounds inlist1[m,n] # if residue(m,n) is in the list
        @inbounds inlist1[m,n] = false   # remove from the list
        @inbounds pos = listaddinv1[m,n]
        if pos == 0
            throw(error("($m,$n) is on the list, but it's not registered."))
        end
        @inbounds listaddinv1[m,n] = 0
        @inbounds for j in pos:listsize1
            listres1[j] = listres1[j+1]
            listadd1[1,j] = listadd1[1,j+1]
            listadd1[2,j] = listadd1[2,j+1]
            if listadd1[1,j] != 0
                listaddinv1[listadd1[1,j],listadd1[2,j]] = j
            end
        end
    end

    @fastmath @inbounds for i in 1:listsize1
        y = listres1[i]
        if x > y
            mm = listadd1[1,end-1]
            nn = listadd1[2,end-1]
            if mm ≠ 0
                inlist1[mm,nn] = false
                listaddinv1[mm,nn] = 0
            end
            for j=listsize1:-1:i+1
                listres1[j] = listres1[j-1]
                listadd1[1,j] = listadd1[1,j-1]
                listadd1[2,j] = listadd1[2,j-1]
                if listadd1[1,j] != 0
                    listaddinv1[listadd1[1,j],listadd1[2,j]] = j
                end
            end
            listadd1[1,i] = m
            listadd1[2,i] = n
            listaddinv1[m,n] = i
            listres1[i] = x
            inlist1[m,n] = true
            break
        end
    end

    @inbounds maxcoords[1] = listadd1[1,1]
    @inbounds maxcoords[2] = listadd1[2,1]

    @inbounds return listres1[1]
end

# 2-List-RBP
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
        listres2::Vector{<:AbstractFloat},
        listadd2::Matrix{<:Integer},
        listaddinv1::Union{Matrix{<:Integer},Nothing},
        listsize1::Integer,
        listsize2::Integer,
        inlist1::Matrix{<:Integer}
    )

    # println("m = $m, n = $n, x = $x, inlist = $(inlist1[m,n])")

    if @inbounds inlist1[m,n] # if residue(m,n) is in the list
        # println()
        # println("### inlist at m = $m and n = $n ###")
        # println()
        @inbounds inlist1[m,n] = false   # remove from the list
        @inbounds pos = listaddinv1[m,n]
        if pos == 0
            throw(error("($m,$n) is on the list, but it's not registered."))
        end
        @inbounds listaddinv1[m,n] = 0
        @inbounds for j in pos:listsize1
            listres1[j] = listres1[j+1]
            listadd1[1,j] = listadd1[1,j+1]
            listadd1[2,j] = listadd1[2,j+1]
            if listadd1[1,j] != 0
                listaddinv1[listadd1[1,j],listadd1[2,j]] = j
            end
        end
    end

    # if @inbounds inlist2[m,n] # if residue(m,n) is in the list
    #     @inbounds inlist2[m,n] = false   # remove from the list
    #     @inbounds pos = listaddinv2[m,n]
    #     if pos == 0
    #         throw(error("($m,$n) is on the list, but it's not registered."))
    #     end
    #     @inbounds listaddinv2[m,n] = 0
    #     @inbounds for j in pos:listsize1
    #         listres2[j] = listres2[j+1]
    #         listadd2[1,j] = listadd2[1,j+1]
    #         listadd2[2,j] = listadd2[2,j+1]
    #         if listadd2[1,j] != 0
    #             listaddinv2[listadd2[1,j],listadd2[2,j]] = j
    #         end
    #     end
    # end

    @fastmath @inbounds for i in 1:listsize1
        y = listres2[i]
        if x > y
            # mm = listadd2[1,end-1]
            # nn = listadd2[2,end-1]
            # if mm ≠ 0
            #     inlist2[mm,nn] = false
            #     listaddinv2[mm,nn] = 0
            # end
            for j=listsize1:-1:i+1
                listres2[j] = listres2[j-1]
                listadd2[1,j] = listadd2[1,j-1]
                listadd2[2,j] = listadd2[2,j-1]
                # if listadd2[1,j] != 0
                #     listaddinv2[listadd2[1,j],listadd2[2,j]] = j
                # end
            end
            listadd2[1,i] = m
            listadd2[2,i] = n
            # listaddinv2[m,n] = i
            listres2[i] = x
            # inlist2[m,n] = true
            break
        end
    end

    # display(listres2)
    # display(listadd2)

    @inbounds return listres1[1]
end