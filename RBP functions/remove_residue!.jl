################################################################################
# Allan Eduardo Feitosa
# 5 Apr 2024
# remove max residue

# List-RBP
function 
    remove_residue!(
        lmax::Int,
        listsize::Int,
        residues::Vector{Float64},
        coords::Matrix{Int},
        inlist::Matrix{Bool},
        pos::Int
    )

    @inbounds inlist[lmax] = false

    # update the list
    @inbounds for i in pos:listsize
        residues[i] = residues[i+1]
        coords[1,i] = coords[1,i+1]
        coords[2,i] = coords[2,i+1]
        coords[3,i] = coords[3,i+1]
    end    
end