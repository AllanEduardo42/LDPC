################################################################################
# Allan Eduardo Feitosa
# 2 Mar 2025
# Update the local list of residues for the List-RBP algorithm

include("add_residue.jl")

# List-RBP
function
    update_local_list!(
        residues::Vector{Float64},
        coords::Matrix{Int},
        localresidues::Vector{Float64},
        localcoords::Matrix{Int},
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
            # remove from list
            inlist[li] = false
            remove_residue!(li,listsizes[1],residues,coords,inlist,pos)
        end                

        add_residue!(nothing,localresidues,localcoords,residue,li,ci,vj,listsizes[2])

    end

end

# RBP
function
    update_local_list!(
        residues::Vector{Float64},
        ::Matrix{Int},
        ::Nothing,
        ::Nothing,
        ::Vector{Int},
        indices::Matrix{Int},
        li::Int,
        ::Int,
        ::Int,
        residue::Float64    
    )

    @inbounds residues[indices[li]] = residue

end

