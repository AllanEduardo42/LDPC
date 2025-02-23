################################################################################
# Allan Eduardo Feitosa
# 23 Feb 2024
# remove edge from list

function 
    remove_from_list!(
        lmax::Integer,
        listsize::Integer,
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer},
        inlist::Matrix{Bool},
        index::Integer
    )

    @inbounds inlist[lmax] = false

    # update the list
    @inbounds for i in index:listsize
        listres[i] = listres[i+1]
        listm[i] = listm[i+1]
        listn[i] = listn[i+1]
    end    
end