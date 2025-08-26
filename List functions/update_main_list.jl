################################################################################
# Allan Eduardo Feitosa
# 2 Mar 2025
# Update the global list of residues for the List-RBP algorithm

include("add_to_list.jl")

#List-RBP
function update_main_list!(
    list::Vector{Float64},
    coords::Matrix{Int},
    local_list::Vector{Float64},
    local_coords::Matrix{Int},
    listsize::Int,
    listsize2::Int,
    inlist::Matrix{Bool}
)
    
    @inbounds begin
        pos = listsize-1
        for i=listsize2:-1:1
            if local_coords[3,i] != 0
                pos = listsize - i + 1
                break
            end
        end
        for i = pos:listsize-1
            li = coords[3,i]
            if li ≠ 0
                inlist[li] = false
                list[i] = 0.0
            else
                break
            end
        end
        for i in 1:listsize2
            m = local_coords[1,i]
            n = local_coords[2,i]
            li = local_coords[3,i]
            add_to_list!(inlist,list,coords,local_list[i],li,m,n,listsize)
        end
        # clear list 2
        local_list .*= 0.0
        local_coords .*= 0
    end
end