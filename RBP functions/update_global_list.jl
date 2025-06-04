################################################################################
# Allan Eduardo Feitosa
# 2 Mar 2025
# Update the global list of residues for the List-RBP algorithm

include("add_residue.jl")

#List-RBP
function update_global_list!(
    residues::Vector{Float64},
    coords::Matrix{Int},
    localresidues::Vector{Float64},
    localcoords::Matrix{Int},
    listsizes::Vector{Int},
    inlist::Matrix{Bool}
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