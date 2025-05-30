################################################################################
# Allan Eduardo Feitosa
# 5 Apr 2024
# remove max residue

# List-RBP
function 
    remove_residue_VN!(
        vjmax::Int,
        listsize::Int,
        alpha::Vector{Float64},
        coords::Vector{Int},
        inlist::Vector{Bool},
        pos::Int
    )

    @inbounds inlist[vjmax] = false

    # update the list
    @inbounds for i in pos:listsize
        alpha[i] = alpha[i+1]
        coords[i] = coords[i+1]
    end    
end