################################################################################
# Allan Eduardo Feitosa
# 2 Mar 2025
# Update the local list of residues for the List-RBP algorithm

include("add_to_list.jl")
include("remove_from_list.jl")

# List-RBP
function
    update_local_list!(
        list::Vector{Float64},
        coords::Matrix{Int},
        local_list::Vector{Float64},
        local_coords::Matrix{Int},
        listsizes::Vector{Int},
        inlist::Matrix{Bool},
        li::Int,
        ci::Int,
        vj::Int,
        residue::Float64    
    )
    
    @inbounds begin
        if inlist[li]  # if residue(ci,vj) is in the list
            # display("($ci,$vj) is on the list")
            pos = 0
            for i = 1:listsizes[1]
                if coords[3,i] == li
                    pos = i
                    break
                end
            end
            # if pos == 0
            #     throw(error("($(coords[1,i]),$(coords[1,i])) is registered as being on the list, but it's not."))
            # end
            remove_from_list!(li,listsizes[1],list,coords,inlist,pos)
        end                

        add_to_list!(nothing,local_list,local_coords,residue,li,ci,vj,listsizes[2])

    end

end

