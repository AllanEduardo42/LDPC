################################################################################
# Allan Eduardo Feitosa
# 25 Feb 2025
# Calculate the residues for the RBP algorithm

include("../update_Lr.jl")
include("_calc_residue.jl")

# RBP
function
    calc_residues!(
        addressinv::Union{Matrix{<:Integer},Nothing},
        residues::Union{Vector{<:AbstractFloat},Nothing},
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        Ms::Matrix{<:AbstractFloat},
        Lr::Union{Matrix{<:AbstractFloat},Nothing},                      
        Lq::Matrix{<:AbstractFloat},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        cnmax::Integer,
        vnmax::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        checks::Union{Vector{<:Integer},Base.OneTo{Int64}}
    )

    for m in checks
        if m ≠ cnmax    
            # calculate the new check to node messages
            update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)
            # calculate the residues
            @fastmath @inbounds for n in cn2vn[m]
                if n ≠ vnmax
                    residue, index = _calc_residue(Ms,Lr,Factors,Lrn,Lq,m,n)
                    residues[addressinv[index]] = residue
                end
            end
        end
    end
end

# Local-RBP
function
    calc_residues!(
        max_residue::Vector{<:AbstractFloat},
        max_coords::Vector{<:Integer},
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        Ms::Matrix{<:AbstractFloat},
        Lr::Union{Matrix{<:AbstractFloat},Nothing},                      
        Lq::Matrix{<:AbstractFloat},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        cnmax::Integer,
        vnmax::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        checks::Union{Vector{<:Integer},Base.OneTo{Int64}}
    )

    for m in checks
        if m ≠ cnmax    
            # calculate the new check to node messages
            update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)
            # calculate the residues
            @fastmath @inbounds for n in cn2vn[m]
                if n ≠ vnmax
                    residue, _ = _calc_residue(Ms,Lr,Factors,Lrn,Lq,m,n)
                    if residue > max_residue[1]
                        max_residue[1], max_residue[2] = residue, max_residue[1]
                        max_coords[3] = max_coords[1]
                        max_coords[4] = max_coords[2]
                        max_coords[1] = m
                        max_coords[2] = n
                    elseif residue > max_residue[2]
                        max_residue[2] = residue
                        max_coords[3] = m
                        max_coords[4] = n
                    end   
                end
            end
        end
    end

    # update list
    if max_residue[1] < max_residue[3]
        max_coords[1], max_coords[5] = max_coords[5], max_coords[1]
        max_coords[2], max_coords[6] = max_coords[6], max_coords[2]
        max_residue[1], max_residue[3] = max_residue[3], max_residue[1]            
    else
        max_residue[3] = max_residue[2]
        max_coords[5] = max_coords[3]
        max_coords[6] = max_coords[4]
    end
end

# List-RBP
function
    calc_residues!(
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        Ms::Matrix{<:AbstractFloat},
        Lr::Union{Matrix{<:AbstractFloat},Nothing},                      
        Lq::Matrix{<:AbstractFloat},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        cnmax::Integer,
        vnmax::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        checks::Union{Vector{<:Integer},Base.OneTo{Int64}},
        listres1::Vector{<:AbstractFloat},
        listm1::Vector{<:Integer},
        listn1::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        listsizes::Vector{<:Integer},
        inlist::Union{Matrix{<:Integer},Nothing}   
    )
    new_listsize2 = listsizes[2]
    for m in checks
        if m ≠ cnmax    
            # calculate the new check to node messages
            update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)
            # calculate the residues
            @fastmath @inbounds for n in cn2vn[m]
                if n ≠ vnmax
                    residue, index = _calc_residue(Ms,Lr,Factors,Lrn,Lq,m,n)
                    if inlist[index]  # if residue(m,n) is in the list
                        # display("($m,$n) is on the list")
                        if listsizes[2] == 1
                            new_listsize2 += 1
                        end
                        pos = 0
                        for i = 1:listsizes[1]
                            if listm1[i] == m
                                if listn1[i] == n
                                    pos = i
                                    break
                                end
                            end
                        end
                        if pos == 0
                            throw(error("($m,$n) is registered as being on the list, but it's not."))
                        end
                        # remove from list
                        inlist[index] = false
                        remove_from_list!(index,listsizes[1],listres1,listm1,
                            listn1,inlist,pos)
                    end                
                    if residue != 0.0
                        if listsizes[2] == 0
                            update_list!(inlist,listres1,listm1,listn1,residue,
                                m,n,listsizes[1])
                        else
                            update_list!(nothing,listres2,listm2,listn2,residue,
                                m,n,new_listsize2)
                        end
                    end
                end
            end
        end
    end

    # update list 1 
    if listsizes[2] ≠ 0
        count = 0
        k = listsizes[1]
        while count < listsizes[2]-1
            count += 1
            k -= 1
            m = listm1[k]
            if m != 0
                listres1[k] = 0.0
                n = listn1[k]
                inlist[m,n] = false
            end
        end
        for k in new_listsize2:-1:1
            update_list!(inlist,listres1,listm1,listn1,listres2[k],listm2[k],
                         listn2[k],listsizes[1])
        end
        listres2 .*= 0.0
        listm2 .*= 0
        listn2 .*= 0
    end
end

