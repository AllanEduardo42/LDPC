################################################################################
# Allan Eduardo Feitosa
# 5 Apr 2024
# remove max residue

# List-RBP
function 
    remove_residue!(
        lmax::Integer,
        listsize::Integer,
        residues::Vector{<:AbstractFloat},
        coords::Matrix{<:Integer},
        inlist::Matrix{Bool},
        pos::Integer
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

#RBP
function 
    remove_residue!(
        ::Integer,
        ::Integer,
        residues::Vector{<:AbstractFloat},
        ::Matrix{<:Integer},
        ::Matrix{<:Integer},
        max_edge::Integer
    )

    @inbounds residues[max_edge] = 0.0
end