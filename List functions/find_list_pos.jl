################################################################################
# Allan Eduardo Feitosa
# 25 Ago 2025

function 
    find_list_pos(
        li::Int,
        listsize::Int,
        coords::Matrix{<:Int},
        ci::Int,
        vj::Int
    )

    pos = 0
    @inbounds for i = 1:listsize
        if coords[3,i] == li
            pos = i
            break
        end
    end
    if pos == 0
        throw(error("($ci,$vj) is registered as being on the list, but it's not."))
    end

    return pos

end