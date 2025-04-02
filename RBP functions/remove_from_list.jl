################################################################################
# Allan Eduardo Feitosa
# 23 Feb 2024
# remove edge from list

function 
    remove_from_list!(
        cnmax::Integer,
        vnmax::Integer,
        listsize::Integer,
        listres::Vector{<:AbstractFloat},
        coords::Matrix{<:Integer},
        inlist::Matrix{Bool},
        pos::Integer
    )

    @inbounds inlist[cnmax,vnmax] = false

    # update the list
    @inbounds for i in pos:listsize
        listres[i] = listres[i+1]
        coords[1,i] = coords[1,i+1]
        coords[2,i] = coords[2,i+1]
    end    
end