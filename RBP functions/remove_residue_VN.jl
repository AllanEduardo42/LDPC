################################################################################
# Allan Eduardo Feitosa
# 5 Apr 2024
# remove max residue

# List-RBP
function 
    remove_residue_VN!(
        vjmax::Integer,
        listsize::Integer,
        alpha::Vector{<:AbstractFloat},
        coords::Vector{<:Integer},
        inlist::Vector{Bool},
        pos::Integer
    )

    @inbounds inlist[vjmax] = false

    # update the list
    @inbounds for i in pos:listsize
        alpha[i] = alpha[i+1]
        coords[i] = coords[i+1]
    end    
end