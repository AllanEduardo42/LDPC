################################################################################
# Allan Eduardo Feitosa
# 24 set 2024
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
        vnmax::Integer,
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )
    # calculate the new check to node messages
    update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)

    # calculate the residues
    @fastmath @inbounds for n in cn2vn[m]
        if n ≠ vnmax
            index = LinearIndices(Ms)[m,n]
            residue = _calc_residue(Ms,Lr,index,Lrn,Lq)
            residue *= Factors[index]
            residues[addressinv[index]] = residue
        end
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
        vnmax::Integer,
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        listres1::Vector{<:AbstractFloat},
        listm1::Vector{<:Integer},
        listn1::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        listsize1::Integer,
        listsize2::Integer,
        new_listsize2::Integer,
        inlist::Union{Matrix{<:Integer},Nothing}   
    )
    # calculate the new check to node messages
    update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)

    # calculate the residues
    @fastmath @inbounds for n in cn2vn[m]
        if n ≠ vnmax

            index = LinearIndices(Ms)[m,n]
            residue = _calc_residue(Ms,Lr,index,Lrn,Lq)
            residue *= Factors[index]

            @fastmath @inbounds if inlist[index]  # if residue(m,n) is in the list
                # display("($m,$n) is on the list")
                inlist[index] = false   # remove from the list
                if listsize2 == 1
                    new_listsize2 += 1
                end
                pos = 0
                for i = 1:listsize1
                    if listm1[i] == m
                        if listn1[i] == n
                            pos = i
                            break
                        end
                    end
                end
                if pos == 0
                    throw(error("($m,$n) is on the list, but it's not registered."))
                end
                # remove from list
                remove_from_list!(index,listsize1,listres1,listm1,listn1,inlist,pos)
            end
        
            @fastmath @inbounds if residue != 0.0
                if listsize2 == 0
                    update_list!(inlist,listres1,listm1,listn1,residue,m,n,listsize1)
                elseif new_listsize2 == 1
                    if residue > listres2[1]
                        listres2[1] = residue
                        listm2[1] = m
                        listn2[1] = n
                    end
                else
                    update_list!(nothing,listres2,listm2,listn2,residue,m,n,new_listsize2)
                end
            end
        end
    end

    return new_listsize2

end

