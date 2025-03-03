################################################################################
# Allan Eduardo Feitosa
# 2 Mar 2025
# Update the list of residues for the List-RBP algorithm

include("add_to_list.jl")

function
    update_list2!(
        listres1::Vector{<:AbstractFloat},
        indices_res1::Vector{<:Integer},
        listres2::Vector{<:AbstractFloat},
        indices_res2::Vector{<:Integer},
        listsizes::Vector{<:Integer},
        new_listsize2::Integer,
        inlist::Matrix{<:Integer},
        li::Integer,
        residue::AbstractFloat    
    )
    
    if inlist[li]  # if residue(m,n) is in the list
        # display("($m,$n) is on the list")
        if listsizes[2] == 1
            new_listsize2 += 1
        end
        pos = 0
        for i = 1:listsizes[1]
            if indices_res1[i] == li
                pos = i
                break
            end
        end
        if pos == 0
            ci = CartesianIndex(inlist)[li]
            throw(error("($(ci[1]),$(ci[2])) is registered as being on the list, but it's not."))
        end
        # remove from list
        inlist[li] = false
        remove_from_list!(li,listsizes[1],listres1,indices_res1,inlist,pos)
    end                
    if residue != 0.0
        add_to_list!(nothing,listres2,indices_res2,residue,li,new_listsize2)
    end

    return new_listsize2
end

function update_list1!(
    listres1::Vector{<:AbstractFloat},
    indices_res1::Vector{<:Integer},
    listres2::Vector{<:AbstractFloat},
    indices_res2::Vector{<:Integer},
    listsizes::Vector{<:Integer},
    new_listsize2::Integer,
    inlist::Matrix{<:Integer},
    listsize1m1::Integer,
    difflistsizes::Integer 
)
    if listsizes[2] ≠ 1
        for i = listsize1m1:-1:difflistsizes
            li = indices_res1[i]
            if li ≠ 0
                inlist[li] = false
                listres1[i] = 0.0
            end
        end
    end
    for i in 1:new_listsize2
        add_to_list!(inlist,listres1,indices_res1,listres2[i],indices_res2[i],
            listsizes[1])
    end
end

