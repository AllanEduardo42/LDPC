################################################################################
# Allan Eduardo Feitosa

include("add_residue_VN.jl")

#List-RBP
function update_global_list_VN!(
    alpha::Vector{Float64},
    coords::Vector{Int},
    localresidues::Vector{Float64},
    localcoords::Vector{Int},
    listsizes::Vector{Int},
    inlist::Vector{Bool}
)
    
    @inbounds begin
        pos = listsizes[1]-1
        for i=listsizes[2]:-1:1
            if localcoords[i] != 0
                pos = listsizes[1] - i + 1
                break
            end
        end
        for i = pos:listsizes[1]-1
            li = coords[i]
            if li ≠ 0
                inlist[li] = false
                alpha[i] = 0.0
            else
                break
            end
        end
        for i in 1:listsizes[2]
            vj = localcoords[i]
            add_residue_VN!(inlist,alpha,coords,localresidues[i],vj,listsizes[1])
        end
        # clear list 2
        localresidues .*= 0.0
        localcoords .*= 0
    end
end