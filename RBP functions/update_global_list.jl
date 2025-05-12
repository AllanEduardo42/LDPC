################################################################################
# Allan Eduardo Feitosa
# 2 Mar 2025
# Update the global list of residues for the List-RBP algorithm

include("add_residue.jl")

#List-RBP
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

# RBP
function update_global_list!(
    ::Vector{<:AbstractFloat},
    ::Matrix{<:Integer},
    ::Nothing,
    ::Nothing,
    ::Vector{<:Integer},
    ::Matrix{<:Integer}
)

end