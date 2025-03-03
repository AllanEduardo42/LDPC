################################################################################
# Allan Eduardo Feitosa
# 23 Feb 2024
# remove edge from list

function 
    remove_from_list!(
        li::Integer,
        listsize::Integer,
        listres::Vector{<:AbstractFloat},
        indices_res::Vector{<:Integer},
        inlist::Matrix{Bool},
        pos::Integer
    )

    @inbounds inlist[li] = false

    # update the list
    @inbounds for i in pos:listsize
        listres[i] = listres[i+1]
        indices_res[i] = indices_res[i+1]
    end    
end