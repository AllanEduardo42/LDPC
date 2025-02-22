################################################################################
# Allan Eduardo Feitosa
# 16 dez 2024
# Set the maximum residue to zero or remove it from the list of residues

#RBP
function 
    set_zero_or_remove!(
        addressinv::Matrix{<:Integer},
        residues::Vector{<:AbstractFloat},
        lmax::Integer,
        ::Integer,
        listres::Vector{<:AbstractFloat},
        ::Vector{<:Integer},
        ::Vector{<:Integer},
        ::Nothing,
        index::Integer
    )

    @inbounds listres[index] = 0
    @inbounds residues[addressinv[lmax]] = 0.0

end

# Local-RBP
function 
    set_zero_or_remove!(
        ::Nothing,
        ::Nothing,
        ::Integer,
        ::Integer,
        listres::Vector{<:AbstractFloat},
        ::Vector{<:Integer},
        ::Vector{<:Integer},
        ::Nothing
    )
    
    @inbounds listres[1] = 0

end

# List-RBP
function 
    set_zero_or_remove!(
        ::Nothing,
        ::Nothing,
        lmax::Integer,
        listsize::Integer,
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer},
        inlist::Matrix{Bool},
        index::Integer
    )

    @inbounds inlist[lmax] = false

    # update the list
    @inbounds for i in index:listsize
        listres[i] = listres[i+1]
        listm[i] = listm[i+1]
        listn[i] = listn[i+1]
    end    
end