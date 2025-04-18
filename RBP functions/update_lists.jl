################################################################################
# Allan Eduardo Feitosa
# 2 Mar 2025
# Update the list of residues for the List-RBP algorithm

include("add_residue.jl")

# List-RBP
function
    update_local_list!(
        residues::Vector{<:AbstractFloat},
        coords::Matrix{<:Integer},
        localresidues::Vector{<:AbstractFloat},
        localcoords::Matrix{<:Integer},
        listsizes::Vector{<:Integer},
        inlist::Matrix{<:Integer},
        li::Integer,
        m::Integer,
        n::Integer,
        residue::AbstractFloat    
    )
    
    @inbounds begin
        if inlist[li]  # if residue(m,n) is in the list
            # display("($m,$n) is on the list")
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

        add_residue!(nothing,localresidues,localcoords,residue,li,m,n,listsizes[2])

    end

end

# RBP
function
    update_local_list!(
        residues::Vector{<:AbstractFloat},
        ::Matrix{<:Integer},
        ::Nothing,
        ::Nothing,
        ::Vector{<:Integer},
        indices::Matrix{<:Integer},
        li::Integer,
        m::Integer,
        n::Integer,
        residue::AbstractFloat    
    )

    @inbounds residues[indices[li]] = residue

end

function update_global_list!(
    residues::Vector{<:AbstractFloat},
    coords::Matrix{<:Integer},
    localresidues::Vector{<:AbstractFloat},
    localcoords::Matrix{<:Integer},
    listsizes::Vector{<:Integer},
    inlist::Matrix{<:Integer}
)
    
    @inbounds begin
        pos = listsizes[1]-1
        for i=listsizes[2]:-1:1
            if localcoords[3,i] != 0
                pos = listsizes[1] - i + 1
                break
            end
        end
        for i = pos:listsizes[1]-1
            li = coords[3,i]
            if li ≠ 0
                inlist[li] = false
                residues[i] = 0.0
            else
                break
            end
        end
        for i in 1:listsizes[2]
            m = localcoords[1,i]
            n = localcoords[2,i]
            li = localcoords[3,i]
            add_residue!(inlist,residues,coords,localresidues[i],li,m,n,listsizes[1])
        end
        # clear list 2
        localresidues .*= 0.0
        localcoords .*= 0
    end
end

function update_global_list!(
    ::Vector{<:AbstractFloat},
    ::Matrix{<:Integer},
    ::Nothing,
    ::Nothing,
    ::Vector{<:Integer},
    ::Matrix{<:Integer}
)

end