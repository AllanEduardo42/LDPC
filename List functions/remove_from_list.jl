################################################################################
# Allan Eduardo Feitosa
# 5 Apr 2024

function 
    remove_from_list!(
        li::Int,
        listsize::Int,
        list::Vector{Float64},
        coords::Matrix{Int},
        inlist::Matrix{Bool},
        pos::Int
    )

    @inbounds inlist[li] = false

    # update the list
    @inbounds for i in pos:listsize
        list[i] = list[i+1]
        coords[1,i] = coords[1,i+1]
        coords[2,i] = coords[2,i+1]
        coords[3,i] = coords[3,i+1]
    end    
end