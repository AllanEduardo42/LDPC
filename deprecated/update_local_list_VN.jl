################################################################################
# Allan Eduardo Feitosa
# 2 Mar 2025
# Update the local list of alpha for the List-RBP algorithm

include("add_residue_VN.jl")

# List-RBP
function
    update_local_list_VN!(
        alpha::Vector{Float64},
        coords::Vector{Int},
        localresidues::Vector{Float64},
        localcoords::Vector{Int},
        listsizes::Vector{Int},
        inlist::Vector{Bool},
        vj::Int,
        residue::Float64    
    )
    
    @inbounds begin
        if inlist[vj]  # if residue(ci,vj) is in the list
            # display("$vj is on the list")
            pos = 0
            for i = 1:listsizes[1]
                if coords[i] == vj
                    pos = i
                    break
                end
            end
            # if pos == 0
            #     throw(error("($(coords[1,i]),$(coords[1,i])) is registered as being on the list, but it's not."))
            # end
            # remove from list
            inlist[vj] = false
            remove_residue_VN!(vj,listsizes[1],alpha,coords,inlist,pos)
        end                

        add_residue_VN!(nothing,localresidues,localcoords,residue,vj,listsizes[2])

    end

end

