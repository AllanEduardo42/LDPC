################################################################################
# Allan Eduardo Feitosa
# 2 Mar 2025
# Update the list of residues for the List-RBP algorithm

include("add_to_list.jl")

function
    update_list2!(
        listres1::Vector{<:AbstractFloat},
        coords1::Matrix{<:Integer},
        listres2::Vector{<:AbstractFloat},
        coords2::Matrix{<:Integer},
        listsizes::Vector{<:Integer},
        new_listsize2::Integer,
        inlist::Matrix{<:Integer},
        m::Integer,
        n::Integer,
        residue::AbstractFloat    
    )
    
    if inlist[m,n]  # if residue(m,n) is in the list
        # display("($m,$n) is on the list")
        if listsizes[2] == 1
            new_listsize2 += 1
        end
        pos = 0
        for i = 1:listsizes[1]
            if coords1[1,i] == m && coords1[2,i] == n
                pos = i
                break
            end
        end
        if pos == 0
            throw(error("($(coords1[1,i]),$(coords1[1,i])) is registered as being on the list, but it's not."))
        end
        # remove from list
        inlist[m,n] = false
        remove_from_list!(m,n,listsizes[1],listres1,coords1,inlist,pos)
    end                
    if residue != 0.0
        add_to_list!(nothing,listres2,coords2,residue,m,n,new_listsize2)
    end

    return new_listsize2
end

function update_list1!(
    listres1::Vector{<:AbstractFloat},
    coords1::Matrix{<:Integer},
    listres2::Vector{<:AbstractFloat},
    coords2::Matrix{<:Integer},
    listsizes::Vector{<:Integer},
    new_listsize2::Integer,
    inlist::Matrix{<:Integer},
    listsize1m1::Integer,
    difflistsizes::Integer 
)
    if listsizes[2] ≠ 1
        for i = listsize1m1:-1:difflistsizes
            m = coords1[1,i]
            n = coords1[2,i]
            if m ≠ 0
                inlist[m,n] = false
                listres1[i] = 0.0
            end
        end
    end
    for i in 1:new_listsize2
        add_to_list!(inlist,listres1,coords1,listres2[i],coords2[1,i],coords2[2,i],
            listsizes[1])
    end
end

